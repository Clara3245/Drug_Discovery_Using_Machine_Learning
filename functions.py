import pandas as pd
import numpy as np
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski


# def search_target(query, filter):
#     """
#     Search for targets using ChEMBL API.
#     """
#     target = new_client.target
#     target_query = target.search(query)
#     if filter == 'y':
#         filtered_target_query = target_query.filter(organism='Homo sapiens')
#     else:
#         filtered_target_query = target_query

#     return pd.DataFrame.from_dict(filtered_target_query)

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

def convert_data(df):
    """
    Convert the result to a DataFrame
    """
    return pd.DataFrame.from_dict(df) 

def preprocess_bioactivity_data(df):
    """
    Preprocess the bioactivity data:
    - Remove rows with missing values
    - Drop duplicates based on 'canonical_smiles'
    - Select relevant columns
    - Convert 'standard_value' to numeric
    """
    df_clean = df.dropna(subset=['canonical_smiles', 'standard_value'])
    df_clean = df_clean.drop_duplicates(['canonical_smiles'])
    selection = ['molecule_chembl_id', 'canonical_smiles', 'standard_value']
    df_clean = df_clean[selection]
    
    # Ensure 'standard_value' is numeric
    df_clean['standard_value'] = pd.to_numeric(df_clean['standard_value'], errors='coerce')
    df_clean.reset_index(drop=True, inplace=True)
    return df_clean


def classify_bioactivity(df, threshold_active=1000, threshold_inactive=10000):
    """
    Classify bioactivity based on IC50 values:
    - Active: IC50 <= threshold_active
    - Inactive: IC50 >= threshold_inactive
    - Intermediate: Otherwise
    """
    bioactivity_threshold = []
    for value in df.standard_value:
        if float(value) >= threshold_inactive:
            bioactivity_threshold.append("inactive")
        elif float(value) <= threshold_active:
            bioactivity_threshold.append("active")
        else:
            bioactivity_threshold.append("intermediate")

    return pd.Series(bioactivity_threshold, name='class')

def clean_canonical_smiles(df):
    smiles = []
    for i in df.canonical_smiles.tolist():
        cpd = str(i).split('.')
        cpd_longest = max(cpd, key = len)
        smiles.append(cpd_longest)

    return pd.Series(smiles, name = 'canonical_smiles')

def lipinski(smiles, verbose=False):

    moldata= []
    for elem in smiles:
        mol=Chem.MolFromSmiles(elem) 
        moldata.append(mol)
       
    baseData= np.arange(1,1)
    i=0  
    for mol in moldata:        
       
        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)
           
        row = np.array([desc_MolWt,
                        desc_MolLogP,
                        desc_NumHDonors,
                        desc_NumHAcceptors])   
    
        if(i==0):
            baseData=row
        else:
            baseData=np.vstack([baseData, row])
        i=i+1      
    
    columnNames=["MW","LogP","NumHDonors","NumHAcceptors"]   
    descriptors = pd.DataFrame(data=baseData,columns=columnNames)
    
    return descriptors

def pIC50(input):
    pIC50 = []

    for i in input['standard_value_norm']:
        molar = i*(10**-9) # Converts nM to M
        pIC50.append(-np.log10(molar))

    input['pIC50'] = pIC50
    x = input.drop('standard_value_norm', 1)
        
    return x

def norm_value(input, norm_value = 100000000):
    norm = []

    for i in input['standard_value']:
        if i > norm_value:
          i = norm_value
        norm.append(i)

    input['standard_value_norm'] = norm
    x = input.drop('standard_value', axis = 1)
        
    return x
