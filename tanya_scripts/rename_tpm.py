import argparse

parser = argparse.ArgumentParser(description='Rename tissue ID in TPM file to sample ID')
parser.add_argument('--infile',required=True,help='Input the path to the infile. This file is the file after running subset_tpm.py. Example is Blood_tpm.tsv')
parser.add_argument('--outfile',required=True,help='Input the path to the output file where the tissue ID is renamed to sample ID')

args = parser.parse_args()

new_file = open(args.outfile, 'w')
with open(args.infile, 'r') as f:
    for line in f:
        if line.startswith('Name'):
            items = line.rstrip('\n').split('\t')
            new_line = [items[0], items[1]]
            for tissue_id in items[2:]:
                i = tissue_id.split('-')
                sample_id = i[0] + '-' + i[1]
                new_line.append(sample_id)
            print ('\t'.join(new_line), file=new_file)
        else:
            print (line.rstrip('\n'), file=new_file)