import pandas as pd
import numpy as np
import subprocess
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from numpy.random import seed
from scipy.stats import mannwhitneyu
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor

def search_target(query):
    """
    Search for targets using the ChEMBL API.
    
    Parameters:
        query (str): Target keyword or name.
        
    Returns:
        pd.DataFrame: Dataframe containing target information.
    """
    target = new_client.target
    target_query = target.search(query)
    return pd.DataFrame.from_dict(target_query)


def fetch_bioactivity_data(target_chembl_id, standard_type="IC50"):
    """
    Fetch bioactivity data for a given target and filter by standard_type.
    
    Parameters:
        target_chembl_id (str): ChEMBL ID of the target.
        standard_type (str): Bioactivity measure type (default is "IC50").
        
    Returns:
        pd.DataFrame: Bioactivity data.
    """
    activity = new_client.activity
    res = activity.filter(target_chembl_id=target_chembl_id).filter(standard_type=standard_type)
    return pd.DataFrame.from_dict(res)


def preprocess_bioactivity_data(df):
    """
    Preprocess the bioactivity data.
    
    Steps:
        - Remove rows with missing values.
        - Drop duplicates based on 'canonical_smiles'.
        - Select relevant columns.
        - Convert 'standard_value' to numeric.
    
    Parameters:
        df (pd.DataFrame): Raw bioactivity data.
        
    Returns:
        pd.DataFrame: Cleaned bioactivity data.
    """
    df_clean = df.dropna(subset=['canonical_smiles', 'standard_value'])
    df_clean = df_clean.drop_duplicates(['canonical_smiles'])
    selection = ['molecule_chembl_id', 'canonical_smiles', 'standard_value']
    df_clean = df_clean[selection]
    df_clean['standard_value'] = pd.to_numeric(df_clean['standard_value'], errors='coerce')
    df_clean.reset_index(drop=True, inplace=True)
    return df_clean


def classify_bioactivity(df, threshold_active=1000, threshold_inactive=10000):
    """
    Classify bioactivity into active, inactive, and intermediate categories.
    
    Parameters:
        df (pd.DataFrame): Dataframe containing 'standard_value'.
        threshold_active (float): Upper limit for active compounds.
        threshold_inactive (float): Lower limit for inactive compounds.
        
    Returns:
        pd.Series: Bioactivity classification.
    """
    bioactivity_threshold = [
        "inactive" if float(value) >= threshold_inactive
        else "active" if float(value) <= threshold_active
        else "intermediate"
        for value in df.standard_value
    ]
    return pd.Series(bioactivity_threshold, name='class')


def clean_canonical_smiles(df):
    """
    Clean and standardize canonical SMILES by retaining the longest fragment.
    
    Parameters:
        df (pd.DataFrame): Dataframe containing 'canonical_smiles'.
        
    Returns:
        pd.Series: Cleaned SMILES.
    """
    smiles = [max(str(i).split('.'), key=len) for i in df.canonical_smiles.tolist()]
    return pd.Series(smiles, name='canonical_smiles')


def lipinski(smiles):
    """
    # Inspired by: https://codeocean.com/explore/capsules?query=tag:data-curation
    Compute Lipinski descriptors for a list of SMILES.
    
    Parameters:
        smiles (list): List of SMILES strings.
        
    Returns:
        pd.DataFrame: Lipinski descriptors (MW, LogP, NumHDonors, NumHAcceptors).
    """
    moldata = [Chem.MolFromSmiles(s) for s in smiles]
    descriptors = [
        [
            Descriptors.MolWt(mol),
            Descriptors.MolLogP(mol),
            Lipinski.NumHDonors(mol),
            Lipinski.NumHAcceptors(mol)
        ]
        for mol in moldata
    ]
    column_names = ["MW", "LogP", "NumHDonors", "NumHAcceptors"]
    return pd.DataFrame(descriptors, columns=column_names)


def pIC50(df):
    """
    Convert IC50 values (nM) to pIC50 (-log10 of molar values).
    
    Parameters:
        df (pd.DataFrame): Dataframe containing 'standard_value_norm'.
        
    Returns:
        pd.DataFrame: Dataframe with a new column 'pIC50'.
    """
    df['pIC50'] = [-np.log10(i * 10**-9) for i in df['standard_value_norm']]
    return df.drop(columns=['standard_value_norm'])


def norm_value(df, max_value=100000000):
    """
    Normalize 'standard_value' by capping at max_value.
    
    Parameters:
        df (pd.DataFrame): Dataframe containing 'standard_value'.
        max_value (float): Maximum value for normalization.
        
    Returns:
        pd.DataFrame: Dataframe with normalized values in 'standard_value_norm'.
    """
    df['standard_value_norm'] = [min(i, max_value) for i in df['standard_value']]
    return df.drop(columns=['standard_value'])


def mannwhitney(descriptor, data):
    """
    # https://machinelearningmastery.com/nonparametric-statistical-significance-tests-in-python/

    Perform Mann-Whitney U Test on descriptor values for active and inactive classes.
    
    Parameters:
        descriptor (str): Descriptor column name.
        data (pd.DataFrame): Dataframe containing descriptors and 'class'.
        
    Returns:
        pd.DataFrame: Results of the test.
    """
    active = data[data['class'] == 'active'][descriptor]
    inactive = data[data['class'] == 'inactive'][descriptor]
    stat, p = mannwhitneyu(active, inactive)
    alpha = 0.05
    interpretation = "Different distribution (reject H0)" if p <= alpha else "Same distribution (fail to reject H0)"
    results = pd.DataFrame({
        'Descriptor': descriptor,
        'Statistics': stat,
        'p': p,
        'alpha': alpha,
        'Interpretation': interpretation
    }, index=[0])
    results.to_csv(f'data/mannwhitneyu_{descriptor}.csv', index=False)
    return results


def run_bash_command(command):
    """
    Execute a bash command and return its output.
    
    Parameters:
        command (str): Bash command to execute.
        
    Returns:
        str: Output or error message.
    """
    print(f"Running command: {command}")
    
    try:
        result = subprocess.run(command, shell=True, capture_output=True, text=True, check=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        return e.stderr


def remove_low_variance(df):
    """
    Remove low-variance features from the input data.
    
    Parameters:
        df (pd.DataFrame): Input dataframe with features and target.
        
    Returns:
        tuple: Filtered features (X) and target (Y).
    """
    X = df.drop('pIC50', axis=1)
    Y = df['pIC50']
    selector = VarianceThreshold(threshold=(.8 * (1 - .8)))
    X = selector.fit_transform(X)
    return X, Y


def define_input_output(df):
    """
    Split data into training and testing sets after removing low-variance features.
    
    Parameters:
        df (pd.DataFrame): Input dataframe with features and target.
        
    Returns:
        tuple: Training and testing sets for features and target.
    """
    X, Y = remove_low_variance(df)
    return train_test_split(X, Y, test_size=0.2)


def build_regression_model(X_train, X_test, Y_train, Y_test):
    """
    Train a Random Forest regression model and evaluate its performance.
    
    Parameters:
        X_train, X_test (array): Training and testing feature sets.
        Y_train, Y_test (array): Training and testing target sets.
        
    Returns:
        tuple: R-squared value and predicted target values.
    """
    model = RandomForestRegressor(n_estimators=100)
    model.fit(X_train, Y_train)
    r2 = model.score(X_test, Y_test)
    Y_pred = model.predict(X_test)
    return r2, Y_pred
