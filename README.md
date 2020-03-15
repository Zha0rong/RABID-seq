# Rabid-seq bioinformatics workflow

## Step 0: Convert paired end reads to a single fastq file
### Use inDrops script to convert a PE file with R1 cell barcode and R2 Rabies barcode, into a SE file containing cell barcode and UMI in the header and rabies barcode as the read.

```diff
-To do: The steps from raw data to SE fastq need to be outlined
```

## Step 1: Extract Rabies barcodes in fastq files
### **Rabie.Barcode.processing.py**: processes all fastq files in folder, generating the following tables: 

```diff
-To do: usage for this script
```

```python
python Rabie.Barcode.processing.py [] [] []

```

Output 1: Cell_statistics.tsv

rabies barcode  | # with only correct 3' handle | # with only correct 5' handle | # with both handles | # with correct structure
--------------- | ------------------------------| ------------------------------| --------------------| --------------

Currently we decide that any read that contains one handle (5' or 3') and the correct barcode structure is kept

Output 2: cell.rabies.counts.tsv

cell barcode  | rabies barcode  | counts
------------- | ----------------| -----

## Step 2: Error correct barcodes
### **run.error.correction.py**

```diff
-To do: usage for this script
```

```python
python run.error.correction.py [] [] []

```

Output: cell.rabies.corrected_counts.tsv

cell barcode  | rabies barcode  | counts
------------- | ----------------| -----

## Step 3: Visualize barcode distributions
### **viz.corrected.barcodes.R**

```diff
-Iain to modify this script to incorporate cell data once available
-To do: usage for this script
```

```R
Rscript viz.corrected.barcodes.R 

```

Barcode_loss.pdf shows the number of correct barcodes over the extraction and correction process

Rabies_histogram.pdf shows # of rabies barcodes per cell

Cell_histogram.pdf shows # of cells per rabies barcodes

Rabies_rannked.pdf shows the ranked abundance of rabies barcode counts per cell

Cell_ranked.pdf shows the ranked abdundance of the number of cells per rabies barcode

## Step 4: Rarefaction curve
### **run.rarefaction.sh**
#### This script subsamples the fastq files and determines the number error corrected barcodes. Plots are generated to determine if the library has been sequenced to sufficient depth

```shell
# IDEA: This is a shell script that subsamples fastq and pipes to the full pipeline above. No need to rewrite the above processing steps
run.rarefaction.sh

```



