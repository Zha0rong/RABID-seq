#!/usr/bin/python
import os
import csv
import subprocess
import sys
from collections import defaultdict, OrderedDict
from itertools import product, combinations
import gzip
import numpy
import pandas as pd
import re

#This function parsed the fastq file into the following format: [Sequence ID, [Sequence],[Quality]]
def ParseFastq(pathstofastqs):
    if pathstofastqs[0].endswith('.gz'):
        processes=[subprocess.Popen(['zcat',(fastq)],stdout=subprocess.PIPE) for fastq in pathstofastqs]
        totalreads = [r.stdout for r in processes]
    elif pathstofastqs[0].endswith('.bz2'):
        processes=[subprocess.Popen(['bzcat',(fastq)],stdout=subprocess.PIPE) for fastq in pathstofastqs]
        totalreads = [r.stdout for r in processes]
    elif pathstofastqs[0].endswith('.fastq'):
        processes=[subprocess.Popen(['cat',(fastq)],stdout=subprocess.PIPE) for fastq in pathstofastqs]
        totalreads = [r.stdout for r in processes]
    else:
        sys.exit('The format of the file %s is not recognized.'%(str(pathtofastq)))
    while True:
        names=[next(read).decode() for read in totalreads]
        Sequence=[next(read).decode() for read in totalreads]
        Blank=[next(read).decode() for read in totalreads]
        qualityscore= [next(read).decode() for read in totalreads]
        assert all(name==names[0] for name in names)
        if names:
            yield [names[0], Sequence, qualityscore]
        else:
            break
    for read in totalreads:
        read.close()

#This function writes read to fastq
def write_fastq(file,ID,seq,quality_score):
    file.write('%s\n'%ID)
    file.write('%s\n' % seq)
    file.write('+\n')
    file.write('%s\n' % quality_score)

class Rabid_Barcode:
    cell_information={}
    #Cell_information dictionary is a dictionary that stores read statistics for each cell, the value for each cell name is a list of 5 integer
    #Total number of reads in cell
	#number of reads,number of reads with 5end handle,number of reads with 3end handle
    #number of reads with only five end handle
    #number of reads with only three end handle
    #number of reads that pass the structure pattern filtering
    whitelist={}
    fastq=''
    samplename=''
    outputfile=''
    fiveendhandle='GCTAGC'
    threeendhandle='GGCGCGCC'
    pattern='[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]AT[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]AT[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]'
    #The pattern is formatted in regex
    def __init__(self,pathtofastq,samplename):
        self.fastq=pathtofastq
        self.samplename=samplename
        self.outputfile=samplename
    def processing(self):
        outputfile=open('%s_filtered.fastq'%(self.outputfile),'wt')#open a file for the filtering reads
        for read in ParseFastq([self.fastq]):
            readid=read[0].strip('\n').split(' ')[0] #Get the sequence ID
            cellname=read[0].strip('\n').split(' ')[1].split(':')#Extract Cell name from Sequence ID
            cellname=cellname[0].strip('\n')+cellname[1]#Extract Cell name from Sequence ID
            umi=read[0].split(' ')[1].split(':')#Extract UMI from Sequence ID
            Rabieread=read[1][0].strip('\n')#Extract The read sequence from Sequence ID
            QualityScore=read[2][0].strip('\n')#Extract The read quality from Sequence ID
            if cellname in self.cell_information.keys():
                self.cell_information[cellname][0]+=1 #If the cell is already in the dictionary, add 1 to the total read in cell statistics of this cell
                if self.fiveendhandle in Rabieread and self.threeendhandle in Rabieread: 
                    self.cell_information[cellname][3]+=1 #If the read has both handle, add 1 to the read with both handle in cell statistics of this cell
                    realRabieread=Rabieread[Rabieread.find(self.fiveendhandle)+len(self.fiveendhandle):Rabieread.find(self.threeendhandle)] #Trim the handle
                    realQualityScore=QualityScore[Rabieread.find(self.fiveendhandle)+len(self.fiveendhandle):Rabieread.find(self.threeendhandle)] #Trim the handle
                    if re.match(self.pattern,realRabieread):
                        self.cell_information[cellname][4]+=1 #Use Regex to check if the read is correct.
                        write_fastq(outputfile,read[0].strip('\n'),realRabieread,realQualityScore) #If correct, write read to the fastq
                else:
                    if self.fiveendhandle in Rabieread: #If only five end handle in the read,
                        self.cell_information[cellname][1]+=1 #add 1 to the read with only five end handle in cell statistics of this cell
                        realRabieread=Rabieread[Rabieread.find(self.fiveendhandle)+len(self.fiveendhandle):Rabieread.find(self.fiveendhandle)+len(self.fiveendhandle)+28] #Get the 28 base pair downstream of the five end handle
                        realQualityScore=QualityScore[Rabieread.find(self.fiveendhandle)+len(self.fiveendhandle):Rabieread.find(self.fiveendhandle)+len(self.fiveendhandle)+28] #Get corresponding quality score
                        if re.match(self.pattern,realRabieread):
                            self.cell_information[cellname][4]+=1 #Use Regex to check if the read is correct.
                            write_fastq(outputfile,read[0].strip('\n'),realRabieread,realQualityScore)
                    elif self.threeendhandle in Rabieread: #If only three end handle in the read,
                        self.cell_information[cellname][2]+=1 #add 1 to the read with only three end handle in cell statistics of this cell
                        realRabieread=Rabieread[Rabieread.find(self.threeendhandle)-28:Rabieread.find(self.threeendhandle)] #Get the 28 base pair upstream of the five end handle
                        realQualityScore=QualityScore[Rabieread.find(self.threeendhandle)-28:Rabieread.find(self.threeendhandle)] #Get corresponding quality score
                        if re.match(self.pattern,realRabieread): #Use Regex to check if the read is correct.
                            self.cell_information[cellname][4]+=1
                            write_fastq(outputfile,read[0].strip('\n'),realRabieread,realQualityScore)
            else:
                self.cell_information[cellname]=[0,0,0,0,0] #If cell not in dictionary, add in cell information
                self.cell_information[cellname][0]+=1 #add 1 to the total read in cell statistics of this cell
                if self.fiveendhandle in Rabieread and self.threeendhandle in Rabieread:
                    self.cell_information[cellname][3]+=1
                    realRabieread=Rabieread[Rabieread.find(self.fiveendhandle)+len(self.fiveendhandle):Rabieread.find(self.threeendhandle)]
                    realQualityScore=QualityScore[Rabieread.find(self.fiveendhandle)+len(self.fiveendhandle):Rabieread.find(self.threeendhandle)]
                    if re.match(self.pattern,realRabieread):
                        self.cell_information[cellname][4]+=1
                        write_fastq(outputfile,read[0].strip('\n'),realRabieread,realQualityScore)
                else:
                    if self.fiveendhandle in Rabieread:
                        self.cell_information[cellname][1]+=1
                        realRabieread=Rabieread[Rabieread.find(self.fiveendhandle)+len(self.fiveendhandle):Rabieread.find(self.fiveendhandle)+len(self.fiveendhandle)+28]
                        realQualityScore=QualityScore[Rabieread.find(self.fiveendhandle)+len(self.fiveendhandle):Rabieread.find(self.fiveendhandle)+len(self.fiveendhandle)+28]
                        if re.match(self.pattern,realRabieread):
                            self.cell_information[cellname][4]+=1
                            write_fastq(outputfile,read[0].strip('\n'),realRabieread,realQualityScore)
                    elif self.threeendhandle in Rabieread:
                        self.cell_information[cellname][2]+=1
                        realRabieread=Rabieread[Rabieread.find(self.threeendhandle)-28:Rabieread.find(self.threeendhandle)]
                        realQualityScore=QualityScore[Rabieread.find(self.threeendhandle)-28:Rabieread.find(self.threeendhandle)]
                        if re.match(self.pattern,realRabieread):
                            self.cell_information[cellname][4]+=1
                            write_fastq(outputfile,read[0].strip('\n'),realRabieread,realQualityScore)
        outputfile.close()
        #write out the cell statistics for each cell
        with open('%s_Cell_statistics.tsv'%(self.outputfile), "w", newline="") as csvfile:
            writer=csv.writer(csvfile,delimiter='\t')
            writer.writerow(['Cellname','number of reads','number of reads with 5end handle','number of reads with 3end handle','number of reads with both handle','number of reads pass the structure filter'])#number of reads, number of reads with both handle, number of reads pass the structure filter
            for cell in self.cell_information.keys():
                writer.writerow([cell,self.cell_information[cell][0],self.cell_information[cell][1],self.cell_information[cell][2],self.cell_information[cell][3],self.cell_information[cell][4]])
        csvfile.close()
        #Use starcode to cluster reads with 1 hamming distance, output results as samplename.clustering.results
        os.system('starcode --print-clusters --dist 1 -i %s_filtered.fastq -o %s.clustering.results'%(self.outputfile,self.outputfile))
        #Read in the clustering results
        starcode=open('%s.clustering.results'%(self.samplename),'r')
        starcode.line=starcode.readlines()
        connection={}
        #For now there are many reads that only has 1 copy and are more than 1 base pair different from every other reads. Now we filter those reads out.
        #In future we come up with another way of fixing the Rabie sequence, which will solve this problem.
        for line in starcode.line:
            #The structure of output of starcode is like this:
            #Cluster center \(correct Rabie Barcide\), tab, number of reads for this cluster, tab, members of this cluster
            connectionitem=line.split(sep='\t')[0]
            time=int(line.split(sep='\t')[1])
            cluster=line.split(sep='\t')[2]
            if time > 1: #Only choose the clusters with more than 1 reads.
                for member in cluster.split(sep=','):
                    connection[member.strip('\n')]=connectionitem.strip('\n')
                    #Here we build a dictionary with structure as member sequence:cluster center sequence, which allow O(1) time efficiecy for assigning reads to cluster.
        quantification={} #Cell statistics for the quantification process
        umi={} #A dictionary with structure as reads:[list of umi], for each read if the umi is already in the list, deem the read as duplicate
        valid_cells=[cell for cell in list(self.cell_information.keys()) if self.cell_information[cell][4]>0] # We only do the quantification for cells that has at least 1 read
        for cell in valid_cells:
            quantification[cell]=[0,0,0]#Total number of reads, number of reads that do not need correction, number of reads need correction
        starcode.close()
        value=connection.values() #Number of clusters in total
        value=list(dict.fromkeys(value))
        umi_stats=connection.keys() #Number of reads in total
        umi_stats=list(dict.fromkeys(umi_stats))
        for values in umi_stats:
            umi[values]=[]
        file_matrix=pd.DataFrame(0,index=list(value),columns=list(valid_cells)) #Use Panda dataframe to build a matrix, with row as Rabie read clusters and column as cell, similar to count matrices
        for read in ParseFastq(pathstofastqs=['%s_filtered.fastq'%(self.outputfile)]):
            cellname=read[0].split(' ')[1].split(':')
            cellname=cellname[0]+cellname[1]
            reads=read[1][0].strip('\n')
            umis=read[0].split(' ')[1].split(':')[2]
            if cellname in valid_cells: #If a cell has more than 1 read
                if reads in connection.keys(): #If the reads belongs to any one of the cluster
                    quantification[cellname][0]+=1 #Add 1 to the total number of reads to cell statstics
                    if connection[reads]==reads: #If the read is the cluster center sequence
                        quantification[cellname][1]+=1 #Add 1 to the number of reads that do not need to be corrected
                    elif connection[reads]!=reads:
                        quantification[cellname][2]+=1 #Else, add 1 to the number of reads that need to be corrected
                    if umis not in umi[reads]: #If the read is not a duplicate
                        file_matrix.loc[connection[reads],cellname]+=1 #Add one to the number of reads in the cluster for this cell
                        umi[connection[reads]].append(umis) #Add the umi to filter out the duplicates in future.
        file_matrix.to_csv('%s.connection.csv'%(self.outputfile), sep=',', encoding='utf-8') #Write Count matrix 
        with open('%s_Cell_quantification_statistics.tsv'%(self.outputfile), "w", newline="") as csvfile:
            writer=csv.writer(csvfile,delimiter='\t')
            writer.writerow(['Cellname','Total Number of reads','number of reads that do not need correction','number of reads need correction'])
            for cell in quantification.keys():
                writer.writerow([cell,quantification[cell][0],quantification[cell][1],quantification[cell][2]])
        csvfile.close()
        #This is a format that removes the zeros in the count matrices, the format looks like:
        #Cell name, Rabie cluster, Number of reads in this cell that belong to this cluster
        with open('%s.connection.ideal.format.csv'%(self.outputfile), "w", newline="") as csvfile:
            writer=csv.writer(csvfile,delimiter='\t')
            writer.writerow(['Cellname','Rabie','Counts'])
            for i in range(len(valid_cells)):
                for j in range(len(value)):
                    if file_matrix.loc[value[j],valid_cells[i]]!=0:
                        writer.writerow([valid_cells[i],value[j],file_matrix.loc[value[j],valid_cells[i]]])
        csvfile.close()
        os.remove('%s_filtered.fastq'%(self.outputfile)) #Remove intermediate fastq files

#Function for rare faction curve, must be called after the processing function, since it needs the clustering results from the processing function.
    def rarefactioncurve(self):
        starcode=open('%s.clustering.results'%(self.samplename),'r')
        starcode.line=starcode.readlines()
        truebarcode={}
        for line in starcode.line:
            connectionitem=line.split(sep='\t')[0]
            time=int(line.split(sep='\t')[1])
            cluster=line.split(sep='\t')[2]
            if time >1:
                for member in cluster.split(sep=','):
                    truebarcode[member.strip('\n')]=connectionitem.strip('\n')
        starcode.close()
        rare_faction={1000:[],5000:[],
        10000:[],50000:[],
        100000:[],500000:[],
        1000000:[]}
        rare_faction_statistics={1000:[0,0,0,0,0],5000:[0,0,0,0,0],
        10000:[0,0,0,0,0],50000:[0,0,0,0,0],
        100000:[0,0,0,0,0],500000:[0,0,0,0,0],
        1000000:[0,0,0,0,0]}#Total number of reads  number of reads number of reads with 5end handle    number of reads with 3end handle    number of reads with both handle    number of reads pass the structure filter
        for read in ParseFastq([self.fastq]):
            readid=read[0].strip('\n').split(' ')[0]
            Rabieread=read[1][0].strip('\n')
            for faction in rare_faction:
                if rare_faction_statistics[faction][0] < faction:
                    random_selector=random.randint(1,101)
                    if random_selector <=50:
                        rare_faction_statistics[faction][0]+=1
                        if self.fiveendhandle in Rabieread and self.threeendhandle in Rabieread:
                            realRabieread=Rabieread[Rabieread.find(self.fiveendhandle)+len(self.fiveendhandle):Rabieread.find(self.threeendhandle)]
                            if re.match(self.pattern,realRabieread):
                                if realRabieread in truebarcode.keys():
                                    rare_faction[faction].append(truebarcode[realRabieread])
                                rare_faction_statistics[faction][3]+=1
                                rare_faction_statistics[faction][4]+=1
                        else:
                            if self.fiveendhandle in Rabieread:
                                realRabieread=Rabieread[Rabieread.find(self.fiveendhandle)+len(self.fiveendhandle):Rabieread.find(self.fiveendhandle)+len(self.fiveendhandle)+28]
                                if re.match(self.pattern,realRabieread):
                                    if realRabieread in truebarcode.keys():
                                        rare_faction[faction].append(truebarcode[realRabieread])
                                    rare_faction_statistics[faction][1]+=1
                                    rare_faction_statistics[faction][4]+=1
                            elif self.threeendhandle in Rabieread:
                                realRabieread=Rabieread[Rabieread.find(self.threeendhandle)-28:Rabieread.find(self.threeendhandle)]
                                if re.match(self.pattern,realRabieread):
                                    if realRabieread in truebarcode.keys():
                                        rare_faction[faction].append(truebarcode[realRabieread])
                                    rare_faction_statistics[faction][2]+=1
                                rare_faction_statistics[faction][4]+=1
            if rare_faction_statistics[1000000][0] > 1000000:
                break
        with open('%s.rarefactioncurve.tsv'%(self.outputfile), "w", newline="") as csvfile:
            writer=csv.writer(csvfile,delimiter='\t')
            writer.writerow(['Number of reads','Unique sequences'])
            for sampling in rare_faction:
                writer.writerow([sampling,len(list(dict.fromkeys(rare_faction[sampling])))])
        csvfile.close()
        with open('%s.rarefactioncurve.statistics.tsv'%(self.outputfile), "w", newline="") as csvfile:
            writer=csv.writer(csvfile,delimiter='\t')
            writer.writerow(['Total number of reads','number of reads with 5end handle','number of reads with 3end handle','number of reads with both handle','number of reads pass the structure filter'])
            for sampling in rare_faction_statistics:
                writer.writerow([rare_faction_statistics[sampling][0],
                rare_faction_statistics[sampling][1],
                rare_faction_statistics[sampling][2],
                rare_faction_statistics[sampling][3],
                rare_faction_statistics[sampling][4]])
        csvfile.close()
        return 0





if __name__=="__main__":
    Process=Rabid_Barcode('S505_STAR_Reads.fastq.gz','S505')
    Process.processing()
    Process.rarefactioncurve()