# RABID-seq bioinformatics pipeline

## Step 1: Convert paired end reads to a single fastq file

#### Use inDrops script to convert a PE file with R1 cell barcode and R2 Rabies barcode, into a SE file containing cell barcode and UMI in the header and rabies barcode as the read.

An example project.yaml file is provided in input/project.yaml 

    python indrops.py project.yaml filter 
      
#### The read structure in this single fastq file will be:
   Read_ID Cellbarcode:UMI  
   Rabiesbarcode  
   +  
   Quality of the Sequence  

## Step 2: Extract, filter, and correct Rabies barcodes
#### Use the extract_rabies_barcodes.py to extract Rabies barcodes, perform structural filtering, and barcode error correction. Specify the mouse id using -m
	python extract_rabies_barcodes.py -m 1
	
#### The structure of the output table is: 
Mouse_SampleID_Cellbarcode | Mouse_Rabiesbarcode | UMI counts
------------ | ------------- | --------------
1_S505_AGACGAGGAGATGGCT	 | 1_ATGTATGTATCTTGCCGTATACATGCAG | 29

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

## Step 4: Generate the igraph network 
#### Use the generate_network.R script to read filter Rabies barcodes and generate a network representation of the data from the output of Step 2 

Example input data files can be found in the input/ folder

	Rscript generate_network.R [--in table.csv] 
		[--out network.RData]
		[--meta metadata.csv]
		[--RNA rna_counts.csv] 
		[--cell_color cell_colors.csv]
		[--cluster_color cluster_colors.csv]

	# --in STR :                              path to csv produced in step 2
	# --out STR :                             path to R object containing the network 
	# --meta STR :                            path to csv containing cell metadata
	# --RNA STR:                              path to transcriptome data
	# --cell_color STR:                       path to cell color file
	# --cluster_color STR:                    path to cluster color file

