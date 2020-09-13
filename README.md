# RABID-seq bioinformatics pipeline

## Scripts used in this pipeline
	RabidSeq.py
	generate_network.R
	viz_network.R
    
## Example input files are found in input/

	rna_counts.csv
	metadata.csv

	network.RData
	cell_colors.csv
	cluster_colors.csv

## Analyze Rabid-Seq data.

### Determine input file type
There are two types of input file formats that can be used in this pipeline. 

#### Case 1: After sequencing inDrops v3 libraries, there may be 3 fastq files per sample.  

	1 file containing Cellbarcode 1 (8 bp), 
	1 file containing Cellbarcode2 + UMI (14 bp),
	1 file containing RNA read (varied length).

This type of data will be referred as 'Raw' in later description. 

**Caution**: Before you decide that your data type is 'Raw', make sure that the names of files do not start with 'Undetermined'. If the names start with undetermined that means the data is not demultiplexed, and you may want to use the [inDrop pipeline](https://github.com/indrops/indrops).

#### Case 2: After sequencing, there may only be 1 fastq file per 1 sample. This means that you have an output file from [inDrop pipeline](https://github.com/indrops/indrops). The inDrop pipeline extracts the cellbarcodes and UMIs from the Raw data type for you and include them in the Read ID in fastq file. The fastq file should look like this:

    Read_ID Cellbarcode1:Cellbarcode2:UMI  
    Rabiesbarcode  
    +  
    Quality of the Sequence  

This type of data will be referred as 'Filtered' in later description.
	
### Step 1: Quantify Rabid-seq data.

Perform UMI counting on Rabies barcodes


#### Case 1: Raw data
For Raw data type (3 fastq files):

    python RabidSeq --quantify_from_inDrop_raw_fastq_files [options] 
   			[-R1 Cellbarcode1.fastq.gz]
    			[-R2 Cellbarcode2andUMI.fastq.gz] 
    			[-R3 Read.fastq.gz]
    			[-o outputdirectory/]
    			[-n outputname]
    Required
       -R1 the name of fastq file that includes the 8bp Cellbarcode 1 (Make sure to include full directory if not in the same directory as script).
       -R2 the name of fastq file that includes the 8bp Cellbarcode 2 and 6bp UMI.
       -R3 the name of fastq file that includes the Read.
       -o the output directory, which will store the output files of the results.
       -n the name of the output. The name will be added to the results files as a prefix.
       
    [option]
       -l the levenshtein distance. The distance is used to correct the Rabid barcode sequencing error. The default distance is 1.

#### Case 2: Raw data 
For Filtered data type (1 fastq file):

    python RabidSeq --quantify_from_inDrop_demultiplexed_fastq_files -R3 Read.fastq.gz -o outputdirectory/ -n outputname 
    Explanation:
       -R3 the name of fastq file that includes the Read. The cell barcode and umi information are already included in the fastq file
       -o the output directory, which will store the output files of the results.
       -n the name of the output. The name will be added to the results files as a prefix.
       [option]
       -l the levenshtein distance. The distance is used to correct the Rabid barcode sequencing error. The default distance is 1.



#### If there are more than one sample per 'sample'...

In our publication for each mouse we have prepared more than 1 sample, and our data shows that the Rabid-seq allows us to find connection between cells in different samples.

If there are more than 1 sample per 'sample', which means that you have multiple samples and you are interested in how the cells in these samples are connected to each other, you need to use an coma-delimited csv sheet which included the name of the overall sample, the names of individual samples and the input files.

We have prepared two templates here for this task, one for Raw data type and one for Filtered data type.

Suppose that we have sample A, B and C and they are from one large 'sample'. If you want to determine the connection of cells between different samples. Here is the sample sheet that you need to have as input. One template is prepared in the repository.

|Read1|Read2|Read3|Overall_Sample_Name|Individual_Sample_Name|
|:---:|:---:|:---:|:---:|:---:|
|A_R1.fastq.gz|A_R2.fastq.gz|A_R3.fastq.gz|Overall_sample1|A|
|B_R1.fastq.gz|B_R2.fastq.gz|B_R3.fastq.gz|Overall_sample1|B|
|C_R1.fastq.gz|C_R2.fastq.gz|C_R3.fastq.gz|Overall_sample1|C|

If you have filtered data type for the samples, here is the sample sheet that you need to have as input. One template is prepared in the repository.

|Read3|Overall_Sample_Name|Individual_Sample_Name|
|:---:|:---:|:---:|
|A_R3.fastq.gz|Overall_sample1|A|
|B_R3.fastq.gz|Overall_sample1|B|
|C_R3.fastq.gz|Overall_sample1|C|

#### Output of Quanfitying step

There will be 4 output files from the Quantifying step. You only need the samplename.table.csv in the next step to generate igraph network. 

##### outputname_filtered.fastq

This will be an intermediate file generated during the pipeline. Feel free to delete it.

#### outputname_statistics.tsv

This will be the cell statistics during the filtering step. The format will look like the table below.

|Cellname*|number of reads|number of reads with 5end handle|number of reads with 3end handle|number of reads with both handle|number of reads pass the structure filter|
|:---:|:---:|:---:|:---:|:---:|:---:|
|CTGTGACCAGCGCCTT|310871|13985|7695|280917|268786|
|...|...|...|...|...|...|

#### samplename.table.csv (This is required for the network generation in step 2)

This will be the quantification output. The format will look like the table below.
|Cellname*|Rabie|Counts|
|:---:|:---:|:---:|
|...|...|...|

#### samplename.clustering.results

This is the output from starcode, which is used to do error correction for the Rabid sequence. 


*If you are using the multiple samples option, the Cellname here will be Overall_sample_name_Individual_sample_name_CTGTGACCAGCGCCTT. The Overall sample name and individual sample name is added to the cell name as prefix in order to avoid cell name collision between samples.

### Step 2: Generate the igraph network 
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

### Step 3: Visualize the igraph network 
#### Generate a graph-based representation of the network with celltypes as vertex colors

Example input data files can be found in the input/ folder

	Rscript viz_network.R --in network.RData
	
	# --in STR :                              path the R object generated by generate_network.R

