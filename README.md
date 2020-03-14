# Rabid-seq bioinformatics workflow

## Step 1: Find Rabies barcodes in fastq files

**Rabie.Barcode.processing.py**: processes all fastq files in folder, generating the following tables: 

Output 1: Cell_statistics.tsv

rabies barcode  | # with only correct 3' handle | # with only correct 5' handle | # with both handles | # with correct structure
--------------- | ------------------------------| ------------------------------| --------------------| --------------

Currently we decide that any read that contains one handle (5' or 3') and the correct barcode structure is kept

Output 2: cell.rabies.count.tsv

cell barcode  | rabies barcode  | count
------------- | ----------------| -----

## Step 2: Error correct barcodes

**run.error.correction.py**

Output: Cell_statistics_corrected.tsv

cell barcode  | # with error correction
------------- | -----------------------

## Step 3: Visualize barcode distributions

**viz.corrected.barcodes.R**

Output: viz.corrected.pdf

## Step 4: Generate rabies count table 

**viz.corrected.barcodes.R**

cell barcode | rabies barcode | # with error correction
------------ | ---------------| ---------------|
