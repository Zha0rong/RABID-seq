# RABID-seq bioinformatics pipeline

## Step 1: Convert paired end reads to a single fastq file

#### Use inDrops script to convert a PE file with R1 cell barcode and R2 Rabies barcode, into a SE file containing cell barcode and UMI in the header and rabies barcode as the read.

An example project.yaml file is provided in input/project.yaml 

    python indrops.py project.yaml filter 
      
#### The output (rabies_se.fastq) will be a single fastq file with the following structure
    
    Read_ID Cellbarcode:UMI  
    Rabiesbarcode  
    +  
    Quality of the Sequence  

## Step 2: Extract, filter, and correct Rabies barcodes
#### Use the extract_rabies_barcodes.py to extract Rabies barcodes, perform structural filtering, and barcode error correction. Specify the mouse id using -m
	
	python extract_rabies_barcodes.py --in rabies_se.fastq -id S505

	# -in STR :                              Path to the single end file produced by "indrops.py project.yaml filter"
	# -id STR :                           	 A unique sample id. If left blank, sample index is used
	
#### The structure of the output table is: 
Mouse_SampleID_Cellbarcode | Mouse_Rabiesbarcode | UMI counts
------------ | ------------- | --------------
AGACGAGGAGATGGCT	 | ATGTATGTATCTTGCCGTATACATGCAG | 29

## Step 3: Perform rarefaction analysis of rabies barcodes
#### Use the run_rarefaction.py script to perform a rarefaction analysis of the Rabies barcodes

	python run_rarefaction.py -1 read1.fastq -2 read2.fastq
	
	# -1 STR :                             Path to unprocessed read1.fastq (same as used by indrops.py project.yaml filter)
	# -2 STR :                             Path to unprocessed read2.fastq (same as used by indrops.py project.yaml filter)
	
#### The structure of the output table is: 
Read depth | Unique Rabies barocdes
------------ | ------------- 
1000 | 5000
10000 | 50000
100000 | 500000
.|.
.|.
.|.


## Step 6: Combine data by mouse
#### Use the combine_samples.R script to combine samples from each mouse
	
	Rscript combine_samples.R in table.csv

table.csv must contain the following information

File | Mouse | Sample ID
-----|------ | ---------



## Step 5: Generate the igraph network 
#### Use the generate_network.R script to read filter Rabies barcodes and generate a network representation of the data from the output of Step 2 

Example input data files can be found in the input/ folder

	Rscript generate_network.R [--in table.csv] 
		[--out network.RData]
		[--meta metadata.csv]
		[--RNA rna_counts.csv] 
		[--cell_color cell_colors.csv]
		[--cluster_color cluster_colors.csv]

	# --in STR :                              Path to csv produced in step 2
	# --out STR :                             Path to R object containing the network 
	# --meta STR :                            Path to csv containing cell metadata
	# --RNA STR:                              Path to transcriptome data
	# --cell_color STR:                       Path to cell color file
	# --cluster_color STR:                    Path to cluster color file

## Step 6: Visualize the igraph network 
#### Generate a graph-based representation of the network with celltypes as vertex colors

Example input data files can be found in the input/ folder

	Rscript viz_network.R --in network.RData
	
	# --in STR :                              path the R object generated by generate_network.R

