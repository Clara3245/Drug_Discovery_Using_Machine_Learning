import os
import pandas as pd
from functions import search_target, fetch_bioactivity_data, convert_data, preprocess_bioactivity_data, classify_bioactivity

'''
If you're working on a new project and it has similar requirements, 
you can reuse the same environment by activating it (run in the terminal: "conda activate bioinformatics").
'''

def main():
    ########### PART 1: DATA CONCISED ###########
    # Create a data folder if it doesn't exist
    data_folder = "data"
    os.makedirs(data_folder, exist_ok=True)

    # Step 1: Search for pancreas-related targets
    targets = search_target('pancreas')
    targets.to_csv(os.path.join(data_folder, 'pancreas_targets.csv'), index=False)

    # Step 2: Select the 7th target (index 6) and fetch bioactivity data
    selected_target = targets.target_chembl_id[6]
    bioactivity_data = fetch_bioactivity_data(selected_target)
    raw_data_path = os.path.join(data_folder, 'pancreaticlipase_data_raw.csv')
    bioactivity_data.to_csv(raw_data_path, index=False)

    # Step 3: Convert the bioactivity data
    data_frame = convert_data(bioactivity_data)
    data_frame_path = os.path.join(data_folder, 'pancreaticlipase_01_bioactivity_data_raw.csv')
    data_frame.to_csv(data_frame_path, index=False)

    # Step 4: Preprocess the bioactivity data
    preprocessed_data = preprocess_bioactivity_data(data_frame)
    preprocessed_data_path = os.path.join(data_folder, 'pancreaticlipase_02_bioactivity_data_preprocessed.csv')
    preprocessed_data.to_csv(preprocessed_data_path, index=False)
    print('preprocessed_data.dtypes')
    print(preprocessed_data.dtypes)   
    print(preprocessed_data)  

    # Step 5: Classify bioactivity based on IC50 values
    bioactivity_classes = classify_bioactivity(preprocessed_data)
    curated_data = pd.concat([preprocessed_data, bioactivity_classes], axis=1)
    curated_data_path = os.path.join(data_folder, 'pancreaticlipase_03_bioactivity_data_curated.csv')
    curated_data.to_csv(curated_data_path, index=False)

    # Print summary Part 1
    print(f"Raw bioactivity data saved to: {data_frame_path}")
    print(data_frame)
    print(f"Preprocessed bioactivity data saved to: {preprocessed_data_path}")
    print(preprocessed_data)
    print(f"Curated bioactivity data with classifications saved to: {curated_data_path}")
    print(curated_data)

    ########### PART 2: EXPLORATORY DATA ANALYSIS ###########
    # df_no_smiles = curated_data.drop(columns='canonical_smiles')

if __name__ == "__main__":
    main()
