# RABID-seq bioinformatics pipeline

## Step 1: Convert paired end reads to a single fastq file

#### Use inDrops script to convert a PE file with R1 cell barcode and R2 Rabies barcode, into a SE file containing cell barcode and UMI in the header and rabies barcode as the read.

An example project.yaml file is provided in input/project.yaml 

    python indrops.py project.yaml filter 
      
#### The read structure in this single fastq file will be:

   Read_ID Cellbarcode:Rabiesbarcode:UMI
   Sequence
   +
   Quality of the Sequence

## Step 2: Extract, filter, and correct Rabies barcodes
#### Use the extract_rabies_barcodes.py to extract Rabies barcodes, perform structural filtering, and barcode error correction. Specify the mouse id using -m
	python extract_rabies_barcodes.py -m 1
	
#### The structure of the output table is: 
Mouse_SampleID_Cellbarcode | Mouse_Rabiesbarcode | UMI counts
------------ | ------------- | --------------
1_S501_A | 1_Rabiesbarcode | 29

## Step 3: Perform rarefaction analysis of rabies barcodes
#### Use the run_rarefaction.py script to perform a rarefaction analysis of the Rabies barcodes

	python run_rarefaction.py -r1 read1.fastq -r2 read2.fastq
	
#### The structure of the output table is: 
Read depth | Unique Rabies barocdes
------------ | ------------- 
1000 | 5000
10000 | 50000
100000 | 500000
.|.
.|.
.|.

## Step 4: Rarefaction curve
### python function to generate data for rarefaction curve
#### This step must only be run after step 2 has been run, since we need to use the clustering results in this process.
This function will use Random Number Generator to subsample the fastq and check the number of unique rabie barcodes at different subsampling size.
### Output at this stage: S505.rarefactioncurve.tsv A tsv file that has the following format:
Number of reads	Unique sequences
S505.rarefactioncurve.statistics.tsv A tsv file that shows how many reads has the handles and pass the filters.



## Step 5: Visualize barcode distributions
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

```shell
# IDEA: This is a shell script that subsamples fastq and pipes to the full pipeline above. No need to rewrite the above processing steps
run.rarefaction.sh

```



