import pandas as pd
from chembl_webresource_client.new_client import new_client

def search_target(query):
    """
    Search for targets using ChEMBL API.
    """
    target = new_client.target
    target_query = target.search(query)
    return pd.DataFrame.from_dict(target_query)

def fetch_bioactivity_data(target_chembl_id, standard_type="IC50"):
    """
    Fetch bioactivity data for a given target and filter by standard_type.
    """
    activity = new_client.activity
    res = activity.filter(target_chembl_id=target_chembl_id).filter(standard_type=standard_type)
    return pd.DataFrame.from_dict(res)

def preprocess_bioactivity_data(df):
    """
    Preprocess the bioactivity data:
    - Remove rows with missing values
    - Drop duplicates based on 'canonical_smiles'
    - Select relevant columns
    """
    df_clean = df.dropna(subset=['canonical_smiles', 'standard_value'])
    df_clean = df_clean.drop_duplicates(['canonical_smiles'])
    selection = ['molecule_chembl_id', 'canonical_smiles', 'standard_value']
    return df_clean[selection]

def classify_bioactivity(df, threshold_active=1000, threshold_inactive=10000):
    """
    Classify bioactivity based on IC50 values:
    - Active: IC50 <= threshold_active
    - Inactive: IC50 >= threshold_inactive
    - Intermediate: Otherwise
    """
    bioactivity_threshold = []
    for value in df['standard_value']:
        if float(value) >= threshold_inactive:
            bioactivity_threshold.append("inactive")
        elif float(value) <= threshold_active:
            bioactivity_threshold.append("active")
        else:
            bioactivity_threshold.append("intermediate")
    return pd.Series(bioactivity_threshold, name='class')
