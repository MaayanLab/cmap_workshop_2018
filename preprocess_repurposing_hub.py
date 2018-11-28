from __future__ import print_function
import pandas as pd

# Load data downloaded from RepurposingHub
drug_moa_df = pd.read_csv('data/Repurposing_Hub_export_2018.txt', sep='\t')
print(drug_moa_df.shape)
drug_moa_df['pert_ids'] = drug_moa_df['Id']\
    .map(lambda x: ','.join(set(['-'.join(item.split('-')[0:2]) for item in x.split(', ')])))

drug_moa_df['pert_id_count'] = drug_moa_df['pert_ids']\
    .map(lambda x: len(x.split(',')))

drug_moa_df.set_index('Name', inplace=True, verify_integrity=True)

# A dict from pert_id to name
d_pert_name = {}
for name, row in drug_moa_df.iterrows():
    for pert_id in row['pert_ids'].split(','):
        d_pert_name[pert_id] = name
len(d_pert_name)

df_pert_name = pd.DataFrame({'pert_id': list(d_pert_name.keys()), 'Name': list(d_pert_name.values())})\
    .set_index('pert_id')

drug_moa_df = df_pert_name.merge(drug_moa_df[['MOA', 'Target', 'Indication', 'Phase']], 
                                 left_on='Name', 
                                 right_index=True
                                )

print(drug_moa_df.shape)

drug_moa_df.to_csv('data/parsed_Repurposing_Hub.csv')

