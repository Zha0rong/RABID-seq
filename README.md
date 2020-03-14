# Rabid-seq bioinformatics workflow

## Step 1: Find Rabies barcodes in fastq files
### **Rabie.Barcode.processing.py**: processes all fastq files in folder, generating the following tables: 

Output 1: Cell_statistics.tsv

rabies barcode  | # with only correct 3' handle | # with only correct 5' handle | # with both handles | # with correct structure
--------------- | ------------------------------| ------------------------------| --------------------| --------------

Currently we decide that any read that contains one handle (5' or 3') and the correct barcode structure is kept

Output 2: cell.rabies.counts.tsv

cell barcode  | rabies barcode  | counts
------------- | ----------------| -----

## Step 2: Error correct barcodes
### **run.error.correction.py**

Output: cell.rabies.corrected_counts.tsv

cell barcode  | rabies barcode  | counts
------------- | ----------------| -----

## Step 3: Visualize barcode distributions
### **viz.corrected.barcodes.R cell.rabies.corrected_counts.tsv

Output: viz.corrected.pdf

## Step 4: Rarefaction curve
### **run.rarefaction.sh**
#### This script subsamples the fastq files and determines the number error corrected barcodes. Plots are generated to determine if the library has been sequenced to sufficient depth


