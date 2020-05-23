import argparse

parser = argparse.ArgumentParser(description='Subset genes with certain tpm threshold')
parser.add_argument('--tpm_threshold',required=True,help='Input the value for the tpm threshold. For example: 1')
parser.add_argument('--tissue_1_tpm',required=True,help='Input the path to tissue 1 tpm. For example: Blood_tpm_sampleID_shared_Putamen.tsv')
parser.add_argument('--tissue_2_tpm',required=True,help='Input the path to tissue 2 tpm. For example: Putamen_tpm_sampleID_shared_Blood.tsv')
parser.add_argument('--tissue_1_out',required=True,help='Input the path to tissue 1 output that contains only genes where tpm across all samples is greater than threshold.')
parser.add_argument('--tissue_2_out',required=True,help='Input the path to tissue 2 output that contains only genes where tpm across all samples is greater than threshold.')

args = parser.parse_args()
threshold = float(args.tpm_threshold)


def CountFrequency(my_list):
    # Creating an empty dictionary
    freq = {}
    for item in my_list:
        if (item in freq):
            freq[item] += 1
        else:
            freq[item] = 1

    for key, value in freq.items():
        if value > 1:
            print (key, value)

tissue_1_genes = set()
with open(args.tissue_1_tpm, 'r') as f:
    for line in f:
        if not line.startswith('Name'):
            tpm_greater_threshold = []
            items = line.rstrip('\n').split('\t')
            for tpm in items[2:]:
                if float(tpm) >= threshold:
                    tpm_greater_threshold.append('True')
                else:
                    tpm_greater_threshold.append('False')
            if all(i=='True' for i in tpm_greater_threshold):
                tissue_1_genes.add(items[1])

tissue_2_genes = set()
with open(args.tissue_2_tpm, 'r') as f:
    for line in f:
        if not line.startswith('Name'):
            tpm_greater_threshold = []
            items = line.rstrip('\n').split('\t')
            for tpm in items[2:]:
                if float(tpm) >= threshold:
                    tpm_greater_threshold.append('True')
                else:
                    tpm_greater_threshold.append('False')
            if all(i=='True' for i in tpm_greater_threshold):
                tissue_2_genes.add(items[1])

genes_with_tpm_threshold = tissue_1_genes.intersection(tissue_2_genes)
print ('Number of genes with tpm greater than threshold of ', str(threshold), ' is: ' , str(len(genes_with_tpm_threshold)))

genes_list = []
genes_set = set()

tissue_1_out = open(args.tissue_1_out, 'w')
with open(args.tissue_1_tpm, 'r') as f:
    for line in f:
        if line.startswith('Name'):
            print (line.rstrip('\n'), file=tissue_1_out)
        else:
            if line.rstrip('\n').split('\t')[1] in genes_with_tpm_threshold:
                genes_list.append(line.rstrip('\n').split('\t')[1])
                genes_set.add(line.rstrip('\n').split('\t')[1])
                print (line.rstrip('\n'), file=tissue_1_out)
print (len(genes_list))
print (len(genes_set))

print (CountFrequency(genes_list))

tissue_2_out = open(args.tissue_2_out, 'w')
with open(args.tissue_2_tpm, 'r') as f:
    for line in f:
        if line.startswith('Name'):
            print (line.rstrip('\n'), file=tissue_2_out)
        else:
            if line.rstrip('\n').split('\t')[1] in genes_with_tpm_threshold:
                print (line.rstrip('\n'), file=tissue_2_out)