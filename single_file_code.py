# Import necessary libraries
import pandas as pd
import numpy as np
import sys
# import seaborn as sns
# import matplotlib.pyplot as plt
from chembl_webresource_client.new_client import new_client
# from rdkit import Chem
# from rdkit.Chem import Descriptors, Lipinski

# Target search for pancreas
target = new_client.target
target_query = target.search('pancreas')  # Search for 'pancreas' in ChEMBL
targets = pd.DataFrame.from_dict(target_query)  # Convert the result to a DataFrame

# Select and retrieve bioactivity data for Human pancreas (7th entry)
selected_target = targets.target_chembl_id[6]  # Select the 7th target (index 6)
activity = new_client.activity  # Initialize activity resource
res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")  # Filter activity data for IC50
df = pd.DataFrame.from_dict(res)  # Convert the result to a DataFrame

# Save raw bioactivity data to CSV
df.to_csv('pancreaticlipase_01_bioactivity_data_raw.csv', index=False)

# Drop rows with missing values in 'canonical_smiles' or 'standard_value'
df2 = df.dropna(subset=['canonical_smiles', 'standard_value'])

# Remove duplicate rows based on 'canonical_smiles'
df2_nr = df2.drop_duplicates(['canonical_smiles'])

# Select relevant columns for further analysis
selection = ['molecule_chembl_id', 'canonical_smiles', 'standard_value']
df3 = df2_nr[selection]  # Create a new DataFrame with selected columns

# Save the preprocessed bioactivity data to CSV
df3.to_csv('pancreaticlipase_02_bioactivity_data_preprocessed.csv', index=False)

# Load the preprocessed data from CSV
df4 = pd.read_csv('pancreaticlipase_02_bioactivity_data_preprocessed.csv')

# Classify bioactivity based on IC50 values
bioactivity_threshold = []  # List to store bioactivity classifications
for i in df4.standard_value:
    if float(i) >= 10000:  # IC50 >= 10000 is considered inactive
        bioactivity_threshold.append("inactive")
    elif float(i) <= 1000:  # IC50 <= 1000 is considered active
        bioactivity_threshold.append("active")
    else:  # IC50 between 1000 and 10000 is considered intermediate
        bioactivity_threshold.append("intermediate")

# Create a new column 'class' for the bioactivity classification
bioactivity_class = pd.Series(bioactivity_threshold, name='class')
df5 = pd.concat([df4, bioactivity_class], axis=1)  # Add the new column to the DataFrame

# Save the curated bioactivity data to CSV
df5.to_csv('pancreaticlipase_03_bioactivity_data_curated.csv', index=False)
