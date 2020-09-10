# RABID-seq bioinformatics pipeline

## Scripts used in this pipeline
    RabidSeq.py
    generate_network.R
    viz_network.R
    

## Example input files are found in input/

~~project.yaml~~
~~read1.fastq~~ 
~~read2.fastq~~
    
~~sample_sheet.csv~~
~~table.csv~~
    
~~network.RData~~
~~metadata.csv~~
~~rna_counts.csv~~
~~cell_colors.csv~~
~~cluster_colors.csv~~


## Starting analyzing the Rabid-Seq data.

### Figuring out what you have in hands.
The first thing one needs to know before using the pipeline is what input he/she is having in hands. It is likely that you have one of these types of input:

1. 3 fastq files for 1 sample. The inDrop sequencing data include 3 fastq files: 

	1 file which includes the Cellbarcode 1 (8 bp), 
	1 file which includes Cellbarcode2 + UMI (14 bp),
	1 file which includes RNA read (varied length).

This type of data will be referred as 'Raw' in later description.

Caution: Before you decide that your data type is 'Raw', make sure that the names of files do not start with 'Undetermined'. If the names start with undetermined that means the data is not demultiplexed, and you may want to use the [inDrop pipeline](https://github.com/indrops/indrops).

2. 1 fastq file for 1 sample. This means that you have an output file from [inDrop pipeline](https://github.com/indrops/indrops). The inDrop pipeline extracts the cellbarcodes and UMIs from the Raw data type for you and include them in the Read ID in fastq file. The fastq file should look like this:

    Read_ID Cellbarcode1:Cellbarcode2:UMI  
    Rabiesbarcode  
    +  
    Quality of the Sequence  

This type of data will be referred as 'Filtered' in later description.
	
### Function: quantifying Rabid-seq data.

Just like any other type of single cell rna-seq data, we need to do quantify read for each cell before proceeding to any kind of analysis.

For Raw data type (3 fastq files):

    python RabidSeq --quantify_from_inDrop_raw_fastq_files [options] -R1 Cellbarcode1.fastq.gz (8bp Cellbarcode 1) -R2 Cellbarcode2andUMI.fastq.gz (8bp Cellbarcode2 and 6bp UMI) -R3 Read.fastq.gz -o outputdirectory/ -n outputname 
    Explanation:
       -R1 the name of fastq file that includes the 8bp Cellbarcode 1 (Make sure to include full directory if not in the same directory as script).
       -R2 the name of fastq file that includes the 8bp Cellbarcode 2 and 6bp UMI.
       -R3 the name of fastq file that includes the Read.
       -o the output directory, which will store the output files of the results.
       -n the name of the output. The name will be added to the results files as a prefix.
       [option]
       -l the levenshtein distance. The distance is used to correct the Rabid barcode sequencing error. The default distance is 1.

For Filtered data type (1 fastq file):

    python RabidSeq --quantify_from_inDrop_demultiplexed_fastq_files -R3 Read.fastq.gz -o outputdirectory/ -n outputname 
    Explanation:
       -R3 the name of fastq file that includes the Read. The cell barcode and umi information are already included in the fastq file
       -o the output directory, which will store the output files of the results.
       -n the name of the output. The name will be added to the results files as a prefix.
       [option]
       -l the levenshtein distance. The distance is used to correct the Rabid barcode sequencing error. The default distance is 1.








~~#### Use inDrops script to convert a PE file with R1 cell barcode and R2 Rabies barcode into a SE file containing cell barcode and UMI in the header and rabies barcode as the read.

~~An example project.yaml file is provided in input/project.yaml 

~~python indrops.py project.yaml filter 
      
~~#### The output (rabies_se.fastq) will be a single fastq file with the following structure
    
~~Read_ID Cellbarcode:UMI  
~~Rabiesbarcode  
~~+  
~~Quality of the Sequence  

~~## Step 2: Extract, filter, and correct Rabies barcodes by mouse
~~#### Use the extract_rabies_barcodes.py to extract Rabies barcodes, perform structural filtering, and barcode error correction. This script combines data by mouse by appending ~~the mouse identification to avoid barcode collisions with the cell barcode and allow rabies connections within mice. 
	
~~python extract_rabies_barcodes.py sample_sheet.csv

s~~ample_sheet.csv must contain the following information

~~File | Mouse | Sample ID
~-----|------ | ---------
~~rabies_se1.fastq | 1 | S503
~~rabies_se2.fastq | 1 | S504
~~rabies_se3.fastq | 3 | S504

~~#### The structure of the output (table.csv) is: 

~~Mouse_SampleID_Cellbarcode | Mouse_Rabiesbarcode            | UMI counts
~-------------------------- | ------------------------------ | --------------
~1_S503_AGACGAGGAGATGGCT	   | 1_ATGTATGTATCTTGCCGTATACATGCAG | 29

~~## Step 3: Perform rarefaction analysis of rabies barcodes
~#### Use the run_rarefaction.py script to perform a rarefaction analysis of the Rabies barcodes
~python run_rarefaction.py -1 read1.fastq -2 read2.fastq
	
~~# -1 STR :                             Path to unprocessed read1.fastq (same as used by indrops.py project.yaml filter)
~# -2 STR :                             Path to unprocessed read2.fastq (same as used by indrops.py project.yaml filter)
	
~~#### The structure of the output table is: 
~Read depth | Unique Rabies barcodes
------------ | ------------- 
1000 | 5000
10000 | 50000
100000 | 500000


## Step 4: Generate the igraph network 
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

## Step 5: Visualize the igraph network 
#### Generate a graph-based representation of the network with celltypes as vertex colors

Example input data files can be found in the input/ folder

	Rscript viz_network.R --in network.RData
	
	# --in STR :                              path the R object generated by generate_network.R

