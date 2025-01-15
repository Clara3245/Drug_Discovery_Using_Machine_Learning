import os
import pandas as pd
from functions import search_target, fetch_bioactivity_data, preprocess_bioactivity_data, classify_bioactivity

'''
If you're working on a new project and it has similar requirements, 
you can reuse the same environment by activating it (run in the terminal: "conda activate bioinformatics").
'''

def main():
    # Create a data folder if it doesn't exist
    data_folder = "data"
    os.makedirs(data_folder, exist_ok=True)

    # Step 1: Search for pancreas-related targets
    targets = search_target('pancreas')
    targets.to_csv(os.path.join(data_folder, 'pancreas_targets.csv'), index=False)

    # Step 2: Select the 7th target (index 6) and fetch bioactivity data
    selected_target = targets.target_chembl_id[6]
    bioactivity_data = fetch_bioactivity_data(selected_target)
    raw_data_path = os.path.join(data_folder, 'pancreaticlipase_01_bioactivity_data_raw.csv')
    bioactivity_data.to_csv(raw_data_path, index=False)

    # Step 3: Preprocess the bioactivity data
    preprocessed_data = preprocess_bioactivity_data(bioactivity_data)
    preprocessed_data_path = os.path.join(data_folder, 'pancreaticlipase_02_bioactivity_data_preprocessed.csv')
    preprocessed_data.to_csv(preprocessed_data_path, index=False)

    # Step 4: Classify bioactivity based on IC50 values
    bioactivity_classes = classify_bioactivity(preprocessed_data)
    curated_data = pd.concat([preprocessed_data, bioactivity_classes], axis=1)
    curated_data_path = os.path.join(data_folder, 'pancreaticlipase_03_bioactivity_data_curated.csv')
    curated_data.to_csv(curated_data_path, index=False)

    # Print summary
    print(f"Raw bioactivity data saved to: {raw_data_path}")
    print(bioactivity_data.head(10))
    print(f"Preprocessed bioactivity data saved to: {preprocessed_data_path}")
    print(preprocessed_data)
    print(f"Curated bioactivity data with classifications saved to: {curated_data_path}")
    print(curated_data)

if __name__ == "__main__":
    main()
