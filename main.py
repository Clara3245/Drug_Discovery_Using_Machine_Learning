import os
import pandas as pd
# from functions import search_target, fetch_bioactivity_data, convert_data, preprocess_bioactivity_data, classify_bioactivity
from functions import *
from plots import *

'''
If you're working on a new project and it has similar requirements, 
you can reuse the same environment by activating it (run in the terminal: "conda activate bioinformatics").
Check environment packages with "conda list"
'''

def main(
    search_label="pancreas",
    target_index=3
):
    ########### PART 1: DATA CONCISED ###########
    # Create a data folder if it doesn't exist
    data_folder = "data"
    os.makedirs(data_folder, exist_ok=True)

    # Step 1: Search for pancreas-related targets
    targets = search_target(search_label)
    targets_data_path = os.path.join(data_folder, f"{search_label}_targets.csv")
    targets.to_csv(targets_data_path, index=False)
    # targets.to_csv(os.path.join(data_folder, (search_label + '_targets.csv')), index=False)

    # Step 2: Select the target and fetch bioactivity data
    selected_target_name = targets.pref_name[target_index]
    formatted_target_name = selected_target_name.replace(" ", "_")
    
    selected_target = targets.target_chembl_id[target_index]
    bioactivity_data = fetch_bioactivity_data(selected_target)
    bioactivity_data_path = os.path.join(data_folder, f"{formatted_target_name}_01_bioactivity_data_raw.csv")
    bioactivity_data.to_csv(bioactivity_data_path, index=False)

    # Step 3: Preprocess the bioactivity data
    preprocessed_data = preprocess_bioactivity_data(bioactivity_data)
    preprocessed_data_path = os.path.join(data_folder, f'{formatted_target_name}_02_bioactivity_data_preprocessed.csv')
    preprocessed_data.to_csv(preprocessed_data_path, index=False)

    # Step 5: Classify bioactivity based on IC50 values
    bioactivity_classes = classify_bioactivity(preprocessed_data)
    curated_data = pd.concat([preprocessed_data, bioactivity_classes], axis=1)
    curated_data_path = os.path.join(data_folder, f'{formatted_target_name}_03_bioactivity_data_curated.csv')
    curated_data.to_csv(curated_data_path, index=False)

    ########### PART 2: EXPLORATORY DATA ANALYSIS ###########
    # Step 1: Clean canonical smiles column
    df_no_smiles = curated_data.drop(columns='canonical_smiles')
    smiles = clean_canonical_smiles(curated_data)
    df_clean_smiles = pd.concat([df_no_smiles,smiles], axis=1)

    # Step 2: Calculate Lipinski descriptors
    df_lipinski = lipinski(df_clean_smiles.canonical_smiles)

    # Step 3: Combine data frames
    df_combined = pd.concat([curated_data,df_lipinski], axis=1)

    # Step 4: Convert IC50 to pIC50
    df_norm = norm_value(df_combined)
    df_final = pIC50(df_norm)
    df_final_path = os.path.join(data_folder, f'{formatted_target_name}_04_bioactivity_data_3class_pIC50.csv')
    df_final.to_csv(df_final_path, index=False)

    # Step 5: Remove intermidiate
    df_2class = df_final[df_final['class'] != 'intermediate']
    df_final_path = os.path.join(data_folder, f'{formatted_target_name}_05_bioactivity_data_2class_pIC50.csv')
    df_2class.to_csv(df_final_path, index=False)

    # Step 6: Plot results
    frequency_plot('class', df_2class)
    scatter_plot('MW','LogP',df_2class,'class','pIC50')

    box_plot('class','pIC50',df_2class)
    box_plot('class','MW',df_2class)
    box_plot('class','LogP',df_2class)
    box_plot('class','NumHDonors',df_2class)
    box_plot('class','NumHAcceptors',df_2class)

    # Step 8: Statistical analysis
    mannwhitney('pIC50', df_2class)
    mannwhitney('MW',df_2class)
    mannwhitney('LogP',df_2class)
    mannwhitney('NumHDonors',df_2class)
    mannwhitney('NumHAcceptors',df_2class) 

    ########### PART 3: DESCRIPTOR DATA PREPARATION ###########
    # Step 1: Load PaDEL-Descriptor
    # run_bash_command("! wget https://github.com/dataprofessor/bioinformatics/raw/master/padel.zip")
    # run_bash_command("! unzip padel.zip")

    selection = ['canonical_smiles','molecule_chembl_id']
    df3_selection = df_final[selection] # df_final = pd.read_csv(f'{data_folder}/{formatted_target_name}_04_bioactivity_data_3class_pIC50.csv')
    df3_selection.to_csv('molecule.smi', sep='\t', index=False, header=False)

    run_bash_command("! bash padel.sh")

    df3_X = pd.read_csv('data/descriptors_output.csv')
    df3_X = df3_X.drop(columns=['Name'])
    df3_Y = df_final['pIC50']
    dataset3 = pd.concat([df3_X,df3_Y], axis=1)
    dataset3.to_csv(f'{data_folder}/{formatted_target_name}_06_bioactivity_data_3class_pIC50_pubchem_fp.csv', index=False)

    print("main ran successfully!")

if __name__ == "__main__":
    main()
