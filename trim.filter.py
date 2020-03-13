#!/usr/bin/python
import os

namelist=os.listdir('./')
namelist=[x for x in namelist if '_STAR_Reads.fastq.gz' in x]
namelist=[x.replace('_STAR_Reads.fastq.gz','') for x in namelist]



for i in range(len(namelist)):
    name=namelist[i]
    #Trim the 5 end handle the *_firstroundtrimming.report.txt is the report of how many reads have handle GATGTCCACGAGGTCTCTGCTAGC
    os.system('cutadapt -g GATGTCCACGAGGTCTCTGCTAGC -e 0.3 %s_STAR_Reads.fastq.gz -o %s.handle.trimmed.fastq.gz > %s_firstroundtrimming.report.txt'%(name,name,name))
    #Trim the 3 end handle the *_secondroundtrimming.report.txt is the report of how many reads have handle GGCGCGCCCGTACGCTGCAGGTCGAC
    os.system('cutadapt -a GGCGCGCCCGTACGCTGCAGGTCGAC -e 0.3 %s.handle.trimmed.fastq.gz -o %s.final.fastq.gz> %s_secondroundtrimming.report.txt'%(name,name,name))
    #Use Regex to remove everything that does not fit in the structure
    os.system('zgrep --no-group-separator -A 2 -B 1 \'[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]AT[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]AT[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]\' %s.final.fastq.gz > %s.filtered.fastq'%(name,name))
    #Use STARCODE to find the correct barcode use clustering
    os.system('starcode --print-clusters --dist 1 -i %s.filtered.fastq -o %s.clustering.results'%(name,name))
    os.remove('%s.handle.trimmed.fastq.gz'%(name))
    os.remove('%s.final.fastq.gz'%(name))