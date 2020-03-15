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


def write_fastq(file,ID,seq,quality_score):
    file.write('%s\n'%ID)
    file.write('%s\n' % seq)
    file.write('+\n')
    file.write('%s\n' % quality_score)

class Rabid_Barcode:
    cell_information={} #Cellname:[number of reads,number of reads with 5end handle,number of reads with 3end handle, number of reads with both handle, number of reads pass the structure filter]
    whitelist={}
    fastq=''
    samplename=''
    outputfile=''
    fiveendhandle='GCTAGC'
    threeendhandle='GGCGCGCC'
    pattern='[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]AT[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]AT[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]'
    def __init__(self,pathtofastq,samplename):
        self.fastq=pathtofastq
        self.samplename=samplename
        self.outputfile=samplename
    def processing(self):
        outputfile=open('%s_filtered.fastq'%(self.outputfile),'wt')
        for read in ParseFastq([self.fastq]):
            readid=read[0].strip('\n').split(' ')[0]
            cellname=read[0].strip('\n').split(' ')[1].split(':')
            cellname=cellname[0].strip('\n')+cellname[1]
            umi=read[0].split(' ')[1].split(':')
            Rabieread=read[1][0].strip('\n')
            QualityScore=read[2][0].strip('\n')
            if cellname in self.cell_information.keys():
                self.cell_information[cellname][0]+=1
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
            else:
                self.cell_information[cellname]=[0,0,0,0,0]
                self.cell_information[cellname][0]+=1
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
        with open('%s_Cell_statistics.tsv'%(self.outputfile), "w", newline="") as csvfile:
            writer=csv.writer(csvfile,delimiter='\t')
            writer.writerow(['Cellname','number of reads','number of reads with 5end handle','number of reads with 3end handle','number of reads with both handle','number of reads pass the structure filter'])#number of reads, number of reads with both handle, number of reads pass the structure filter
            for cell in self.cell_information.keys():
                writer.writerow([cell,self.cell_information[cell][0],self.cell_information[cell][1],self.cell_information[cell][2],self.cell_information[cell][3],self.cell_information[cell][4]])
        csvfile.close()
        os.system('starcode --print-clusters --dist 1 -i %s_filtered.fastq -o %s.clustering.results'%(self.outputfile,self.outputfile))
        starcode=open('%s.clustering.results'%(samplename),'r')
        starcode.line=starcode.readlines()
        connection={}
        for line in starcode.line:
            connectionitem=line.split(sep='\t')[0]
            time=int(line.split(sep='\t')[1])
            cluster=line.split(sep='\t')[2]
            if time > 1:
                for member in cluster.split(sep=','):
                    connection[member.strip('\n')]=connectionitem.strip('\n')
        quantification={}
        umi={}
        valid_cells=[cell for cell in list(self.cell_information.keys()) if self.cell_information[cell][4]>0]
        starcode.close()
        value=connection.values()
        value=list(dict.fromkeys(value))
        for values in value:
            quantification[values]=0
            umi[values]=[]
        file_matrix=pd.DataFrame(0,index=list(value),columns=list(valid_cells))
        for read in ParseFastq(pathstofastqs=['%s_filtered.fastq'%(self.outputfile)]):
            cellname=read[0].split(' ')[1].split(':')
            cellname=cellname[0]+cellname[1]
            reads=read[1][0].strip('\n')
            umis=read[0].split(' ')[1].split(':')[2]
            if cellname in valid_cells:
                if reads in connection.keys():
                    if umis not in umi[connection[reads]]:
                        file_matrix.loc[connection[reads],cellname]+=1
                        umi[connection[reads]].append(umis)
        file_matrix.to_csv('%s.connection.csv'%(self.outputfile), sep=',', encoding='utf-8')
        os.remove('%s_filtered.fastq'%(self.outputfile))
        
        
        #Quantification without correction
        uncorrected.rabie.barcodes=[]
        starcode=open('%s.clustering.results'%(samplename),'r')
        starcode.line=starcode.readlines()
        for line in starcode.line:
            connectionitem=line.split(sep='\t')[0]
            time=int(line.split(sep='\t')[1])
            cluster=line.split(sep='\t')[2]
            for member in cluster.split(sep=','):
                if member not in uncorrected.rabie.barcodes:
                    uncorrected.rabie.barcodes.append(member.strip('\n'))
            if connectionitem not in uncorrected.rabie.barcodes:
                uncorrected.rabie.barcodes.append(connectionitem.strip('\n'))
        uncorrected.quantification={}
        uncorrected.umi={}
        for values in uncorrected.rabie.barcodes:
            uncorrected.quantification[values]=0
            uncorrected.umi[values]=[]  
        file_matrix=pd.DataFrame(0,index=uncorrected.rabie.barcodes,columns=list(valid_cells))
        for read in ParseFastq(pathstofastqs=['%s_filtered.fastq'%(self.outputfile)]):
            cellname=read[0].split(' ')[1].split(':')
            cellname=cellname[0]+cellname[1]
            reads=read[1][0].strip('\n')
            umis=read[0].split(' ')[1].split(':')[2]
            if cellname in valid_cells:
                if reads in uncorrected.rabie.barcodes:
                    if umis not in umi[reads]:
                        file_matrix.loc[reads,cellname]+=1
                        uncorrected.umi[reads].append(umis)
        file_matrix.to_csv('%s.connection.csv'%(self.outputfile), sep=',', encoding='utf-8')
        
        
        
        
        
        
if __name__=="__main__":
    samplelist=os.listdir('./')
    samplelist=[x for x in samplelist if '_STAR_Reads.fastq.gz' in x]
    samplelist=[x.replace('_STAR_Reads.fastq.gz','') for x in samplelist]
    for sample in samplelist:
        Process=Rabid_Barcode('%s_STAR_Reads.fastq.gz'%(sample),sample)
        Process.processing()
    
    
    
    
    
    
    
    
    
    
    
    
    