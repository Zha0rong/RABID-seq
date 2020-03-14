# Rabid-seq bioinformatics workflow

## Step 1: Find Rabies barcodes in fastq files

**Rabie.Barcode.processing.py**: processes all fastq files in folder, generating the following tables: 

Output: Cell_statistics.tsv

barcode  | # with only correct 3' handle | # with only correct 5' handle | # with both handles | # with correct structure
-------- | --------------------------| --------------------------| --------------------------| --------------------------


**error.correct.barcodes.py**

## Step 2: Quality control 

### viz.corrected.barcodes.R

Output: viz.corrected.pdf
