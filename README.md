# RABID-seq bioinformatics pipeline

## Step 1: Convert paired end reads to a single fastq file

#### Use inDrops script to convert a PE file with R1 cell barcode and R2 Rabies barcode, into a SE file containing cell barcode and UMI in the header and rabies barcode as the read.

An example project.yaml file is provided in input/project.yaml 

    python indrops.py project.yaml filter 
      
#### The read structure in this single fastq file will be:

Read_ID CellBarcode_1:RABIESBarcode_2:UMI
Sequence
+
Quality of the Sequence

## Step 2: Extract, filter and correct Rabies barcodes
#### Use the extract_rabies_barcodes.py to extract Rabies barcodes, perform structural filtering, and barcode error correction. Specify the mouse id using -m
	python extract_rabies_barcodes.py -m 1
	
#### The structure of the output table is: 
1_Index_Cellbarocde 	\t 	1_Rabiesbarcode \t Counts

## Step 2: Use the starcode to cluster the reads from sample_filtered.fastq
The hamming distance we used is 1.
The output format is like this: 
Cluster center sequence (correct barcode sequence) tab Number of reads in this cluster tab members sequences of cluster (the sequences that are only 1 hamming distance away from the cluster center).
### Output at this stage: text file called sample.clustering.results

## Step 3: Use the sample_filtered.fastq and sample.clustering.results to build the count matrix for the rabie barcode
### At this stage, what we do first is to read in the sample.clustering.results: 
	remove any cluster that only has 1 read.
### Use the filtered.fastq, I correct the rabie barcode and count the number of reads belong to each barcode in each cell. (More specific mechanism in the code annotation)

### Output at this stage: 
Count Matrix: S505.connection.csv this is a count matrix similar to the transcriptome output.
Count Table: S505.connection.ideal.format.csv this is the output that remove the 0s in the count matrix, the format is like this
Cell	Rabie Barcode	Count
Quantification statistics: S505_Cell_quantification_statistics.tsv this statistics shows how many reads per cell, and how many reads need to be corrected.

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



