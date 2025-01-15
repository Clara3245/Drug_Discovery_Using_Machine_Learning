# Import necessary libraries
import pandas as pd
from chembl_webresource_client.new_client import new_client

# Target search for coronavirus
target = new_client.target
target_query = target.search('pancreas')
targets = pd.DataFrame.from_dict(target_query)
print('\nTarget search for pancreas')
print(targets.head(10))

selected_target = targets.target_chembl_id[6]
print('\nSelected target')
print(selected_target)

activity = new_client.activity
# res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
res = activity.filter(target_chembl_id=selected_target)


df = pd.DataFrame.from_dict(res)
print('\nBioactivity data for Human Pancreatic lipase (CHEMBL1812) that are reported as pChEMBL values')
print(df)

df.to_csv('pancreas_07_bioactivity_data_raw.csv', index=False)