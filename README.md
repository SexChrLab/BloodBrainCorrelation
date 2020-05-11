# BloodBrainCorrelation
Testing for correlations between expression in the blood and brain

## Tanya Phung's script

### 1. Rename tissue ID in TPM file to sample ID
- Rename tissue ID in TPM file to sample ID in order to select the individuals that have TPM values when comparing 2 tissues
- Use the script `rename_tpm.py`:
  ```
  python rename_tpm.py -h
  usage: rename_tpm.py [-h] --infile INFILE --outfile OUTFILE

  Rename tissue ID in TPM file to sample ID

  optional arguments:
  -h, --help         show this help message and exit
  --infile INFILE    Input the path to the infile. This file is the file after
                     running subset_tpm.py. Example is Blood_tpm.tsv
  --outfile OUTFILE  Input the path to the output file where the tissue ID is
                     renamed to sample ID
  ```
  - Example:
  ```
  python rename_tpm.py --infile Blood_tpm.tsv --outfile Blood_tpm_sampleID.tsv
  python rename_tpm.py --infile Putamen_tpm.tsv --outfile Putamen_tpm_sampleID.tsv
  ```

### 2. Subset the tissue tpm files to contain samples that are shared between 2 tissues
- Use the script `subset_tissue_tpm_for_shared_samples.py`:
  ```
  python subset_tissue_tpm_for_shared_samples.py -h
  usage: subset_tissue_tpm_for_shared_samples.py [-h] --tissue_1_tpm
                                               TISSUE_1_TPM --tissue_2_tpm
                                               TISSUE_2_TPM
                                               --tissue_1_tpm_shared
                                               TISSUE_1_TPM_SHARED
                                               --tissue_2_tpm_shared
                                               TISSUE_2_TPM_SHARED

  Subset tissue tpm files to contain samples that are shared between 2 tissues

  optional arguments:
  -h, --help            show this help message and exit
  --tissue_1_tpm TISSUE_1_TPM
                        Input the path to tissue 1 tpm file after renaming.
                        For example: Blood_tpm_sampleID.tsv
  --tissue_2_tpm TISSUE_2_TPM
                        Input the path to tissue 2 tpm file after renaming.
                        For example: Putamen_tpm_sampleID.tsv
  --tissue_1_tpm_shared TISSUE_1_TPM_SHARED
                        Input the output file for tissue 1 that contains
                        shared samples
  --tissue_2_tpm_shared TISSUE_2_TPM_SHARED
                        Input the output file for tissue 2 that contains
                        shared samples
  ```
  - Example:
  ```
  python subset_tissue_tpm_for_shared_samples.py --tissue_1_tpm Blood_tpm_sampleID.tsv --tissue_2_tpm Putamen_tpm_sampleID.tsv --tissue_1_tpm_shared Blood_tpm_sampleID_shared_Blood.tsv --tissue_2_tpm_shared Putamen_tpm_sampleID_shared_Putamen.csv
  ```
    - In addition to 2 new output files `Blood_tpm_sampleID_shared_Blood.tsv` and `Putamen_tpm_sampleID_shared_Putamen.csv`, this script prints out:
    ```
    Number of shared samples is  160
    Number of columns of tissue 1 after subsetting for shared samples is  162
    Number of columns of tissue 2 after subsetting for shared samples is  162
    ```
      - So you know that the number of shared samples match. 
      
 ### 3. Subset genes with certain tpm threshold
 - Use the script `subset_genes_with_tpm_threshold.py`:
  ```
  python subset_genes_with_tpm_threshold.py -h
  usage: subset_genes_with_tpm_threshold.py [-h] --tpm_threshold TPM_THRESHOLD
                                          --tissue_1_tpm TISSUE_1_TPM
                                          --tissue_2_tpm TISSUE_2_TPM
                                          --tissue_1_out TISSUE_1_OUT
                                          --tissue_2_out TISSUE_2_OUT

  Subset genes with certain tpm threshold

  optional arguments:
  -h, --help            show this help message and exit
  --tpm_threshold TPM_THRESHOLD
                        Input the value for the tpm threshold. For example: 1
  --tissue_1_tpm TISSUE_1_TPM
                        Input the path to tissue 1 tpm. For example:
                        Blood_tpm_sampleID_shared_Putamen.tsv
  --tissue_2_tpm TISSUE_2_TPM
                        Input the path to tissue 2 tpm. For example:
                        Putamen_tpm_sampleID_shared_Blood.tsv
  --tissue_1_out TISSUE_1_OUT
                        Input the path to tissue 1 output that contains only
                        genes where tpm across all samples is greater than
                        threshold.
  --tissue_2_out TISSUE_2_OUT
                        Input the path to tissue 2 output that contains only
                        genes where tpm across all samples is greater than
                        threshold.
  ```
  - Example:
   ```
   python subset_genes_with_tpm_threshold.py --tpm_threshold 1 --tissue_1_tpm Blood_tpm_sampleID_shared_Blood.tsv --tissue_2_tpm Putamen_tpm_sampleID_shared_Putamen.tsv --tissue_1_out Blood_tpm_sampleID_shared_Blood_genes_tpm_greater_than_1.tsv --tissue_2_out Putamen_tpm_sampleID_shared_Putamen_genes_tpm_greater_than_1.tsv
   ```
   - This script also prints out: 
   ```
   Number of genes with tpm greater than threshold of  1.0  is:  3243
   ```
