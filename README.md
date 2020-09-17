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

#### Case 1: After sequencing inDrops v3 libraries, there may be 3 fastq files per sample. This type of data will be referred as 'Raw' in later description. 

	1 file containing Cellbarcode 1 (8 bp), 
	1 file containing Cellbarcode2 + UMI (14 bp),
	1 file containing RNA read (varied length).


**Caution**: Before you decide that your data type is 'Raw', make sure that the names of files do not start with 'Undetermined'. If the names start with undetermined that means the data is not demultiplexed, and you may want to use the [inDrop pipeline](https://github.com/indrops/indrops).

#### Case 2: After sequencing, there may only be 1 fastq file per 1 sample. This type of data will be referred as 'Filtered' in later description.
In this case, you have an output file from [inDrop pipeline](https://github.com/indrops/indrops). The inDrop pipeline extracts the cellbarcodes and UMIs from the Raw data type for you and include them in the Read ID in fastq file. The fastq file should look like this:

    Read_ID Cellbarcode1:Cellbarcode2:UMI  
    Rabiesbarcode  
    +  
    Quality of the Sequence  

	
### Step 1: Prepare your sample sheet 

#### Combining samples that correspond to a single mouse
To ensure unique cell barcoding, cells from a single mouse are divided into samples for sequencing and must be combined at the time of analysis. This allows all cell-cell connections to be determined within each mouse. 

The sample sheet is a **coma-delimited csv** file with the following information: 

	1. Overall_Sample_Name = Sample names of the biological unit (mice)
	2. Individual_Sample_Name = Sample names of individual library samples (sequencing libraries) 
	3. Read<1,2,3> = Names of the input files (sequencing files)

A different samples sheet format is needed depending on the type of input files you have. Please refer to the two templates we have provided for **Case 1: Raw data and Case 2: Filtered data.**. Suppose that you have **sample A, B and C** from **Mouse 1**. The sample sheet will look like:

#### Case 1: Raw data
|Read1|Read2|Read3|Overall_Sample_Name|Individual_Sample_Name|
|:---:|:---:|:---:|:---:|:---:|
|A_R1.fastq.gz|A_R2.fastq.gz|A_R3.fastq.gz|Overall_sample1|A|
|B_R1.fastq.gz|B_R2.fastq.gz|B_R3.fastq.gz|Overall_sample1|B|
|C_R1.fastq.gz|C_R2.fastq.gz|C_R3.fastq.gz|Overall_sample1|C|

#### Case 2: Filtered data
|Read3|Overall_Sample_Name|Individual_Sample_Name|
|:---:|:---:|:---:|
|A_R3.fastq.gz|Overall_sample1|A|
|B_R3.fastq.gz|Overall_sample1|B|
|C_R3.fastq.gz|Overall_sample1|C|


### Step 2: Quantify Rabid-seq data.
Now that the sample sheets are prepared, perform UMI counting on Rabies barcodes. 

#### Case 1: Raw data
For Raw data type (3 fastq files):

    python RabidSeq --quantify_from_inDrop_raw_fastq_files [option]
    			-R1 Cellbarcode1.fastq.gz
    			-R2 Cellbarcode2andUMI.fastq.gz
    			-R3 Read.fastq.gz
    			-o outputdirectory/
    			-n outputname
    Required
    #   -R1 the name of fastq file that includes the 8bp Cellbarcode 1 (Make sure to include full directory if not in the same directory as script).
    #   -R2 the name of fastq file that includes the 8bp Cellbarcode 2 and 6bp UMI.
    #   -R3 the name of fastq file that includes the Read.
    #   -o the output directory, which will store the output files of the results.
    #   -n the name of the output. The name will be added to the results files as a prefix.
       
    [option]
    #   -l the levenshtein distance. The distance is used to correct the Rabid barcode sequencing error. The default distance is 1.

#### Case 2: Filtered data 
For Filtered data type (1 fastq file):

    python RabidSeq --quantify_from_inDrop_demultiplexed_fastq_files [option]
    			-R3 Read.fastq.gz
    			-o outputdirectory/ 
    			-n outputname 
    
    Required
     #  -R3 the name of fastq file that includes the Read. The cell barcode and umi information are already included in the fastq file
     #  -o the output directory, which will store the output files of the results.
     #  -n the name of the output. The name will be added to the results files as a prefix.
     
    [option]
       -l the levenshtein distance. The distance is used to correct the Rabid barcode sequencing error. The default distance is 1.

#### Case 3: Combining samples into one large 'sample'

    python RabidSeq --quantify_from_multiple_samples [option]
    			-s samplesheet.csv
    			-o outputdirectory/ 
    			-n outputname 
    
    Required
     #  -s the name of sample sheet introduced in step 1.
     #  -o the output directory, which will store the output files of the results.
     #  -n the name of the output. The name will be added to the results files as a prefix.
     
    [option]
       -l the levenshtein distance. The distance is used to correct the Rabid barcode sequencing error. The default distance is 1.

#### Output of Quanfitying step

There will be 4 output files from the Quantifying step. You only need the samplename.table.csv in the next step to generate igraph network. 

**outputname_filtered.fastq** - intermediate file generated during the pipeline. Feel free to delete it.

	
**outputname_statistics.tsv** -  statistics from filtering step. The format is: 

|Cellname*|number of reads|number of reads with 5end handle|number of reads with 3end handle|number of reads with both handle|number of reads pass the structure filter|
|:---:|:---:|:---:|:---:|:---:|:---:|
|Overall_sample_Individual_sample_CTGTGACCAGCGCCTT|310871|13985|7695|280917|268786|
|...|...|...|...|...|...|


**samplename.table.csv** - UMI counts (This is required for the network generation in step 2). This will be the quantification output. The format is: 
	
|Cellname*|Rabie|Counts|
|:---:|:---:|:---:|
|Overall_sample_Individual_sample_CTGTGACCAGCGCCTT|...|...|

**samplename.clustering.results** - Barcode clustering (Starcode) output. This is used to perform error correction on Rabies barcode sequences. 

*NOTE: Cellname here will be Overall_sample_name_Individual_sample_name_CTGTGACCAGCGCCTT. The Overall sample name (mouse) and individual sample name (sequencing library) is added to the cell name as prefix in order to avoid cell barcode collisions between mice.
