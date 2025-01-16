import os
import pandas as pd
# from functions import search_target, fetch_bioactivity_data, convert_data, preprocess_bioactivity_data, classify_bioactivity
from functions import *

'''
If you're working on a new project and it has similar requirements, 
you can reuse the same environment by activating it (run in the terminal: "conda activate bioinformatics").
Check environment packages with "conda list"
'''

def main():
    ########### PART 1: DATA CONCISED ###########
    # Create a data folder if it doesn't exist
    data_folder = "data"
    os.makedirs(data_folder, exist_ok=True)

    # Step 1: Search for pancreas-related targets
    search_label = input('Enter search label: \n')
    targets = search_target(search_label)
    targets.to_csv(os.path.join(data_folder, (search_label + '_targets.csv')), index=False)

    # Step 2: Select the target and fetch bioactivity data
    target_index = 3
    selected_target = targets.target_chembl_id[target_index]
    bioactivity_data = fetch_bioactivity_data(selected_target)
    selected_target_name = targets.pref_name[target_index]
    formatted_target_name = selected_target_name.replace(" ", "_")

    # Step 3: Convert the bioactivity data
    data_frame = convert_data(bioactivity_data)
    data_frame_path = os.path.join(data_folder, formatted_target_name + '_01_bioactivity_data_raw.csv')
    data_frame.to_csv(data_frame_path, index=False)

    # Step 4: Preprocess the bioactivity data
    preprocessed_data = preprocess_bioactivity_data(data_frame)
    preprocessed_data_path = os.path.join(data_folder, formatted_target_name + '_02_bioactivity_data_preprocessed.csv')
    preprocessed_data.to_csv(preprocessed_data_path, index=False)

    # Step 5: Classify bioactivity based on IC50 values
    bioactivity_classes = classify_bioactivity(preprocessed_data)
    curated_data = pd.concat([preprocessed_data, bioactivity_classes], axis=1)
    curated_data_path = os.path.join(data_folder, formatted_target_name + '_03_bioactivity_data_curated.csv')
    curated_data.to_csv(curated_data_path, index=False)

    # Print summary Part 1
    print(targets)
    print(f"Raw bioactivity data saved to: {data_frame_path}")
    print(data_frame)
    print(f"Preprocessed bioactivity data saved to: {preprocessed_data_path}")
    print(preprocessed_data)
    print(f"Curated bioactivity data with classifications saved to: {curated_data_path}")
    print(curated_data)

    # ########### PART 2: EXPLORATORY DATA ANALYSIS ###########
    # # Step 1: Clean canonical smiles column
    # df_no_smiles = curated_data.drop(columns='canonical_smiles')
    # smiles = clean_canonical_smiles(curated_data)
    # df_clean_smiles = pd.concat([df_no_smiles,smiles], axis=1)

    # # Step 2: Calculate Lipinski descriptors
    # df_lipinski = lipinski(df_clean_smiles.canonical_smiles)

    # # Step 3: Combine data frames
    # df_combined = pd.concat([curated_data,df_lipinski], axis=1)

    # # Step 4: Convert IC50 to pIC50
    # df_norm = norm_value(df_combined)
    # df_final = pIC50(df_norm)

    # # Step 5: Save data
    # df_final.to_csv('_04_bioactivity_data_3class_pIC50.csv')



if __name__ == "__main__":
    main()
