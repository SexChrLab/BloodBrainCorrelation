import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Subset tissue tpm files to contain samples that are shared between 2 tissues')
parser.add_argument('--tissue_1_tpm',required=True,help='Input the path to tissue 1 tpm file after renaming. For example: Blood_tpm_sampleID.tsv')
parser.add_argument('--tissue_2_tpm',required=True,help='Input the path to tissue 2 tpm file after renaming. For example: Putamen_tpm_sampleID.tsv')
parser.add_argument('--tissue_1_tpm_shared',required=True,help='Input the output file for tissue 1 that contains shared samples')
parser.add_argument('--tissue_2_tpm_shared',required=True,help='Input the output file for tissue 2 that contains shared samples')

args = parser.parse_args()

# Find shared sampleIDs
tissue_1_id = set()
with open(args.tissue_1_tpm, 'r') as f:
    items = f.readline().rstrip('\n').split('\t')
    for tissue_id in items[2:]:
        tissue_1_id.add(tissue_id)

tissue_2_id = set()
with open(args.tissue_2_tpm, 'r') as f:
    items = f.readline().rstrip('\n').split('\t')
    for tissue_id in items[2:]:
        tissue_2_id.add(tissue_id)

shared_id = tissue_1_id.intersection(tissue_2_id)
print ('Number of shared samples is ', str(len(shared_id)))

# Subset each tissue tpm to contain samples in shared
tissue_1_tpm = pd.read_csv(args.tissue_1_tpm, sep='\t', chunksize=1000)
tissue_1_chunk_list = []
for data_chunk in tissue_1_tpm:
    data_chunk_subset = data_chunk[['Name', 'Description']]
    for sample in shared_id:
        sample_subset = data_chunk[[sample]]
        data_chunk_subset = data_chunk_subset.join(sample_subset)
    tissue_1_chunk_list.append(data_chunk_subset)

tissue_1_tpm_subset = pd.concat(tissue_1_chunk_list) #concatenate all the chunks together
print ('Number of columns of tissue 1 after subsetting for shared samples is ', str(tissue_1_tpm_subset.shape[1]))
tissue_1_tpm_subset.to_csv(args.tissue_1_tpm_shared, sep='\t', index=False) #save to a file

tissue_2_tpm = pd.read_csv(args.tissue_2_tpm, sep='\t', chunksize=1000)
tissue_2_chunk_list = []
for data_chunk in tissue_2_tpm:
    data_chunk_subset = data_chunk[['Name', 'Description']]
    for sample in shared_id:
        sample_subset = data_chunk[[sample]]
        data_chunk_subset = data_chunk_subset.join(sample_subset)
    tissue_2_chunk_list.append(data_chunk_subset)

tissue_2_tpm_subset = pd.concat(tissue_2_chunk_list) #concatenate all the chunks together
print ('Number of columns of tissue 2 after subsetting for shared samples is ', str(tissue_2_tpm_subset.shape[1]))
tissue_2_tpm_subset.to_csv(args.tissue_2_tpm_shared, sep='\t', index=False) #save to a file