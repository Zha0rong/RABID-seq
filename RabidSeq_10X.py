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
import bz2
import argparse
import random

parser = argparse.ArgumentParser(prog = "Rabid Seq pipeline",usage="RNAseq pipeline.")



parser.add_argument("--quantify_from_inDrop_raw_fastq_files",dest='quantify_from_inDrop_fastq_files',action="store_true",help='Quantify from 3 inDrop files')
parser.add_argument("--quantify_from_inDrop_demultiplexed_fastq_files",dest='quantify_from_filtered_fastq_files',action="store_true",help='Quantify from 1 inDrop fastq files')
parser.add_argument("--quantify_from_multiple_samples",dest='quantify_from_multiple_samples',action="store_true",help='Quantify from multiple samples.')
parser.add_argument("--check_library_diversity",dest='check_library_diversity',action="store_true",help='Check Rabid barcodes diversity in Plasmid Library.')

parser.add_argument('-R1','--Cellbarcode1',dest='Cellbarcode1',action="store",required="--quantify_from_inDrop_fastq_files" in sys.argv)
parser.add_argument('-R3','--Read',dest='Read',action="store",required="--quantify_from_inDrop_fastq_files" in sys.argv
                                                                       or "--quantify_from_filtered_fastq_files" in sys.argv or "--check_library_diversity" in sys.argv)
parser.add_argument('-o','--output',dest='outputdirectory',action="store",required="--quantify_from_inDrop_fastq_files" in sys.argv
                                                                                   or "--quantify_from_filtered_fastq_files" in sys.argv or "--check_library_diversity" in sys.argv
                    or "--quantify_from_multiple_samples" in sys.argv)
parser.add_argument('-n','--name',dest='name',action="store",required="--quantify_from_inDrop_fastq_files" in sys.argv
                                                                      or "--quantify_from_filtered_fastq_files" in sys.argv
                    or "--quantify_from_multiple_samples" in sys.argv  or "--check_library_diversity" in sys.argv)
parser.add_argument('-s','--sheet',dest='sheet',action="store",required="--quantify_from_multiple_samples" in sys.argv)

parser.add_argument('-l','--levenshtein',dest='distance',action='store',default=1)


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
        sys.exit('The format of the file %s is not recognized.'%(str(pathstofastqs)))
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


class Rabid_Seq_Processor:
    '''This class is used to read in the inDrop fastq data (3 gzipped file: Cell Barcode 1, Cell Barcode 2 + UMI, Rabid (RNA) Reads).
       The output will be '''
    CB1=''
    Read=''
    samplename=''
    outputdirectory=''
    fiveendhandle = 'GCTAGC'
    threeendhandle = 'GGCGCGCC'
    pattern = '[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]AT[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]AT[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]'
    cell_information = {}  # Cellname:[number of reads,number of reads with 5end handle,number of reads with 3end handle, number of reads with both handle, number of reads pass the structure filter]
    whitelist = {}
    distance=0
    def __init__(self,CB,Read,samplename,outputdirectory,distance):
        self.CB=CB
        self.Read = Read
        self.samplename = samplename
        self.outputdirectory = outputdirectory
        self.distance=distance
        #if os.path.isdir(self.outputdirectory) is False:
            #sys.exit('RabidSeq pipeline exiting, the output directory does not exist')
        #if os.path.isfile(self.CB1) is False:
            #sys.exit('RabidSeq pipeline exiting, the CB1 fastq file does not exist')
        #if os.path.isfile(self.CB2_UMI) is False:
            #sys.exit('RabidSeq pipeline exiting, the CB2+UMI fastq file does not exist')
        #if os.path.isfile(self.Read) is False:
            #sys.exit('RabidSeq pipeline exiting, the Read fastq file does not exist')
        #if self.CB1.endswith('.gz') and self.CB2_UMI.endswith('.gz') and self.Read.endswith('.gz') is False:
            #sys.exit('RabidSeq pipeline exiting, the fastq files needs to be gzipped (Update supporting other format coming soon).')

    def _write_fastq(self,file, ID, seq, quality_score):
        file.write('%s\n' % ID)
        file.write('%s\n' % seq)
        file.write('+\n')
        file.write('%s\n' % quality_score)
    def _ParseFastq(self):
        Allfiles=[self.CB,self.Read]
        Allfiles = [gzip.open(files) for files in Allfiles]
        while True:
            try:
                names = [next(read).decode().split(' ')[0].strip('\n') for read in Allfiles]
                Sequence = [next(read).decode().strip('\n') for read in Allfiles]
                Blank = [next(read).decode().strip('\n') for read in Allfiles]
                qualityscore = [next(read).decode().strip('\n') for read in Allfiles]
                assert all(name.strip('\n') == names[0].strip('\n') for name in names)
                if names:
                    try:
                        yield [names[0].split(' ')[0], Sequence, qualityscore]
                    except:
                        return
                else:
                    break
            except StopIteration:
               break
        for read in Allfiles:
            read.close()
    def _Parse_filtered_fastq(self,filtered_fastq=None):
        if filtered_fastq is None:
             fastq = ['%s/%s_filtered.fastq' % (self.outputdirectory, self.samplename)]
        else:
            fastq=[filtered_fastq]
        fastq = [open(files) for files in fastq]
        while True:
            try:
                names = [next(read).strip('\n') for read in fastq]
                Sequence = [next(read).strip('\n') for read in fastq]
                Blank = [next(read).strip('\n') for read in fastq]
                qualityscore = [next(read).strip('\n') for read in fastq]
                assert all(name.strip('\n') == names[0].strip('\n') for name in names)
                if names:
                    try:
                        yield [names[0], Sequence, qualityscore]
                    except:
                        return
                else:
                    break
            except StopIteration:
               break
        for read in fastq:
            read.close()
    def Extract_and_Filtering(self):
        outputfile = open('%s/%s_filtered.fastq' % (self.outputdirectory,self.samplename), 'wt',encoding='utf-8')
        for read in self._ParseFastq():
            name=read[0]
            CellBarcode=read[1][0][0:16]
            UMI=read[1][0][16::]
            Rabid_sequence=read[1][1]
            Rabid_sequence_quality=read[2][1]
            cellname=CellBarcode
            name=name+' '+CellBarcode+':'+UMI
            if cellname in self.cell_information.keys():
                self.cell_information[cellname][0] += 1
                if self.fiveendhandle in Rabid_sequence and self.threeendhandle in Rabid_sequence:
                    self.cell_information[cellname][3] += 1
                    realRabieread = Rabid_sequence[
                                    Rabid_sequence.find(self.fiveendhandle) + len(self.fiveendhandle):Rabid_sequence.find(
                                        self.threeendhandle)]
                    realQualityScore = Rabid_sequence_quality[
                                       Rabid_sequence.find(self.fiveendhandle) + len(self.fiveendhandle):Rabid_sequence.find(
                                           self.threeendhandle)]
                    if re.match(self.pattern, realRabieread):
                        self.cell_information[cellname][4] += 1
                        self._write_fastq(outputfile, name, realRabieread, realQualityScore)
                else:
                    if self.fiveendhandle in Rabid_sequence:
                        self.cell_information[cellname][1] += 1
                        realRabieread = Rabid_sequence[
                                        Rabid_sequence.find(self.fiveendhandle) + len(self.fiveendhandle):Rabid_sequence.find(
                                            self.fiveendhandle) + len(self.fiveendhandle) + 28]
                        realQualityScore = Rabid_sequence_quality[
                                           Rabid_sequence.find(self.fiveendhandle) + len(self.fiveendhandle):Rabid_sequence.find(
                                               self.fiveendhandle) + len(self.fiveendhandle) + 28]
                        if re.match(self.pattern, realRabieread):
                            self.cell_information[cellname][4] += 1
                            self._write_fastq(outputfile, name, realRabieread, realQualityScore)
                    elif self.threeendhandle in Rabid_sequence:
                        self.cell_information[cellname][2] += 1
                        realRabieread = Rabid_sequence[
                                        Rabid_sequence.find(self.threeendhandle) - 28:Rabid_sequence.find(self.threeendhandle)]
                        realQualityScore = Rabid_sequence_quality[
                                           Rabid_sequence.find(self.threeendhandle) - 28:Rabid_sequence.find(self.threeendhandle)]
                        if re.match(self.pattern, realRabieread):
                            self.cell_information[cellname][4] += 1
                            self._write_fastq(outputfile, name, realRabieread, realQualityScore)
            else:
                self.cell_information[cellname]=[0,0,0,0,0]
                self.cell_information[cellname][0]+=1
                if self.fiveendhandle in Rabid_sequence and self.threeendhandle in Rabid_sequence:
                    self.cell_information[cellname][3]+=1
                    realRabieread=Rabid_sequence[Rabid_sequence.find(self.fiveendhandle)+len(self.fiveendhandle):Rabid_sequence.find(self.threeendhandle)]
                    realQualityScore=Rabid_sequence_quality[Rabid_sequence.find(self.fiveendhandle)+len(self.fiveendhandle):Rabid_sequence.find(self.threeendhandle)]
                    if re.match(self.pattern,realRabieread):
                        self.cell_information[cellname][4]+=1
                        self._write_fastq(outputfile,name,realRabieread,realQualityScore)
                else:
                    if self.fiveendhandle in Rabid_sequence:
                        self.cell_information[cellname][1]+=1
                        realRabieread=Rabid_sequence[Rabid_sequence.find(self.fiveendhandle)+len(self.fiveendhandle):Rabid_sequence.find(self.fiveendhandle)+len(self.fiveendhandle)+28]
                        realQualityScore=Rabid_sequence_quality[Rabid_sequence.find(self.fiveendhandle)+len(self.fiveendhandle):Rabid_sequence.find(self.fiveendhandle)+len(self.fiveendhandle)+28]
                        if re.match(self.pattern,realRabieread):
                            self.cell_information[cellname][4]+=1
                            self._write_fastq(outputfile,name,realRabieread,realQualityScore)
                    elif self.threeendhandle in Rabid_sequence:
                        self.cell_information[cellname][2]+=1
                        realRabieread=Rabid_sequence[Rabid_sequence.find(self.threeendhandle)-28:Rabid_sequence.find(self.threeendhandle)]
                        realQualityScore=Rabid_sequence_quality[Rabid_sequence.find(self.threeendhandle)-28:Rabid_sequence.find(self.threeendhandle)]
                        if re.match(self.pattern,realRabieread):
                            self.cell_information[cellname][4]+=1
                            self._write_fastq(outputfile,name,realRabieread,realQualityScore)
        outputfile.close()
        with open('%s/%s_Cell_statistics.tsv'%(self.outputdirectory,self.samplename), "w", newline="") as csvfile:
            writer=csv.writer(csvfile,delimiter='\t')
            writer.writerow(['Cellname','number of reads','number of reads with 5end handle','number of reads with 3end handle','number of reads with both handle','number of reads pass the structure filter'])#number of reads, number of reads with both handle, number of reads pass the structure filter
            for cell in self.cell_information.keys():
                writer.writerow([cell,self.cell_information[cell][0],self.cell_information[cell][1],self.cell_information[cell][2],self.cell_information[cell][3],self.cell_information[cell][4]])
        csvfile.close()
    def Filtering_from_filtered_fastq(self):
        outputfile = open('%s/%s_filtered.fastq' % (self.outputdirectory,self.samplename), 'wt')
        for read in self._Parse_filtered_fastq(filtered_fastq=self.Read):
            name=read[0]
            CellBarcode=read[0].split(sep=' ')[1].split(sep=':')[0]
            Rabid_sequence=read[1][0]
            Rabid_sequence_quality=read[2][0]
            cellname=CellBarcode
            if cellname in self.cell_information.keys():
                self.cell_information[cellname][0] += 1
                if self.fiveendhandle in Rabid_sequence and self.threeendhandle in Rabid_sequence:
                    self.cell_information[cellname][3] += 1
                    realRabieread = Rabid_sequence[
                                    Rabid_sequence.find(self.fiveendhandle) + len(
                                        self.fiveendhandle):Rabid_sequence.find(
                                        self.threeendhandle)]
                    realQualityScore = Rabid_sequence_quality[
                                       Rabid_sequence.find(self.fiveendhandle) + len(
                                           self.fiveendhandle):Rabid_sequence.find(
                                           self.threeendhandle)]
                    if re.match(self.pattern, realRabieread):
                        self.cell_information[cellname][4] += 1
                        self._write_fastq(outputfile, name, realRabieread, realQualityScore)
                else:
                    if self.fiveendhandle in Rabid_sequence:
                        self.cell_information[cellname][1] += 1
                        realRabieread = Rabid_sequence[
                                        Rabid_sequence.find(self.fiveendhandle) + len(
                                            self.fiveendhandle):Rabid_sequence.find(
                                            self.fiveendhandle) + len(self.fiveendhandle) + 28]
                        realQualityScore = Rabid_sequence_quality[
                                           Rabid_sequence.find(self.fiveendhandle) + len(
                                               self.fiveendhandle):Rabid_sequence.find(
                                               self.fiveendhandle) + len(self.fiveendhandle) + 28]
                        if re.match(self.pattern, realRabieread):
                            self.cell_information[cellname][4] += 1
                            self._write_fastq(outputfile, name, realRabieread, realQualityScore)
                    elif self.threeendhandle in Rabid_sequence:
                        self.cell_information[cellname][2] += 1
                        realRabieread = Rabid_sequence[
                                        Rabid_sequence.find(self.threeendhandle) - 28:Rabid_sequence.find(
                                            self.threeendhandle)]
                        realQualityScore = Rabid_sequence_quality[
                                           Rabid_sequence.find(self.threeendhandle) - 28:Rabid_sequence.find(
                                               self.threeendhandle)]
                        if re.match(self.pattern, realRabieread):
                            self.cell_information[cellname][4] += 1
                            self._write_fastq(outputfile, name, realRabieread, realQualityScore)
            else:
                self.cell_information[cellname] = [0, 0, 0, 0, 0]
                self.cell_information[cellname][0] += 1
                if self.fiveendhandle in Rabid_sequence and self.threeendhandle in Rabid_sequence:
                    self.cell_information[cellname][3] += 1
                    realRabieread = Rabid_sequence[Rabid_sequence.find(self.fiveendhandle) + len(
                        self.fiveendhandle):Rabid_sequence.find(self.threeendhandle)]
                    realQualityScore = Rabid_sequence_quality[Rabid_sequence.find(self.fiveendhandle) + len(
                        self.fiveendhandle):Rabid_sequence.find(self.threeendhandle)]
                    if re.match(self.pattern, realRabieread):
                        self.cell_information[cellname][4] += 1
                        self._write_fastq(outputfile, name, realRabieread, realQualityScore)
                else:
                    if self.fiveendhandle in Rabid_sequence:
                        self.cell_information[cellname][1] += 1
                        realRabieread = Rabid_sequence[Rabid_sequence.find(self.fiveendhandle) + len(
                            self.fiveendhandle):Rabid_sequence.find(self.fiveendhandle) + len(self.fiveendhandle) + 28]
                        realQualityScore = Rabid_sequence_quality[Rabid_sequence.find(self.fiveendhandle) + len(
                            self.fiveendhandle):Rabid_sequence.find(self.fiveendhandle) + len(self.fiveendhandle) + 28]
                        if re.match(self.pattern, realRabieread):
                            self.cell_information[cellname][4] += 1
                            self._write_fastq(outputfile, name, realRabieread, realQualityScore)
                    elif self.threeendhandle in Rabid_sequence:
                        self.cell_information[cellname][2] += 1
                        realRabieread = Rabid_sequence[
                                        Rabid_sequence.find(self.threeendhandle) - 28:Rabid_sequence.find(
                                            self.threeendhandle)]
                        realQualityScore = Rabid_sequence_quality[
                                           Rabid_sequence.find(self.threeendhandle) - 28:Rabid_sequence.find(
                                               self.threeendhandle)]
                        if re.match(self.pattern, realRabieread):
                            self.cell_information[cellname][4] += 1
                            self._write_fastq(outputfile, name, realRabieread, realQualityScore)
        outputfile.close()
        with open('%s/%s_Cell_statistics.tsv' % (self.outputdirectory, self.samplename), "w", newline="") as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow(
                ['Cellname', 'number of reads', 'number of reads with 5end handle', 'number of reads with 3end handle',
                 'number of reads with both handle',
                 'number of reads pass the structure filter'])  # number of reads, number of reads with both handle, number of reads pass the structure filter
            for cell in self.cell_information.keys():
                writer.writerow([cell, self.cell_information[cell][0], self.cell_information[cell][1],
                                 self.cell_information[cell][2], self.cell_information[cell][3],
                                 self.cell_information[cell][4]])
        csvfile.close()

    def Clustering(self):
        distance=self.distance
        os.system(
            'starcode --print-clusters --dist %s -i %s/%s_filtered.fastq -o %s/%s.clustering.results' % (
            str(distance),self.outputdirectory, self.samplename,self.outputdirectory, self.samplename))
    def Correction_and_generate_table(self):
        if os.path.isfile('%s/%s.clustering.results'%(self.outputdirectory, self.samplename)) is False:
            sys.exit('RabidSeq pipeline exiting, the Clustering results from STARCODE is missing.')
        starcode=open('%s/%s.clustering.results'%(self.outputdirectory, self.samplename),'r')
        starcode.line=starcode.readlines()
        Rabid_correction_database={}
        for line in starcode.line:
            Rabid_Cluster_Center=line.split(sep='\t')[0]
            Rabid_Frequency=int(line.split(sep='\t')[1])
            Rabid_Cluster=line.split(sep='\t')[2]
            if Rabid_Frequency > 1:
                for member in Rabid_Cluster.split(sep=','):
                    Rabid_correction_database[member.strip('\n')]=Rabid_Cluster_Center.strip('\n')
        umi_list={} #For demulitplexing
        valid_cells=[cell for cell in list(self.cell_information.keys()) if self.cell_information[cell][4]>0]
        starcode.close()
        Rabid=Rabid_correction_database.values()
        Rabid=list(dict.fromkeys(Rabid))

        for rabie in Rabid:
            umi_list[rabie]=[]
        file_matrix=pd.DataFrame(0,index=list(Rabid),columns=list(valid_cells))

        for read in self._Parse_filtered_fastq():
            CELLNAME=read[0].split(sep=' ')[1].split(sep=':')[0]+read[0].split(sep=' ')[1].split(sep=':')[1]
            UMI=read[0].split(sep=' ')[1].split(sep=':')[2]
            UMI=CELLNAME+UMI
            reads=read[1][0].strip('\n')
            if CELLNAME in valid_cells:
                if reads in Rabid_correction_database.keys():
                    if UMI not in umi_list[Rabid_correction_database[reads]]:
                        file_matrix.loc[Rabid_correction_database[reads],CELLNAME]+=1
                        umi_list[Rabid_correction_database[reads]].append(UMI)
        with open('%s/%s.table.csv'%(self.outputdirectory,self.samplename), "w", newline="") as csvfile:
            writer=csv.writer(csvfile,delimiter='\t')
            writer.writerow(['Cellname','Rabid','Counts'])
            for i in range(len(valid_cells)):
                for j in range(len(Rabid)):
                    if file_matrix.loc[Rabid[j],valid_cells[i]]!=0:
                        writer.writerow([valid_cells[i],Rabid[j],file_matrix.loc[Rabid[j],valid_cells[i]]])
        csvfile.close()

class multi_Rabid_Seq_Processer:
    Overall_sample_name=''
    Individual_sample={}
    output_dir=''
    datatype='' '''Raw or Filterd'''
    distance=1
    cell_information = {}
    fiveendhandle = 'GCTAGC'
    threeendhandle = 'GGCGCGCC'
    pattern = '[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]AT[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]AT[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]'
    def _write_fastq(self,file, ID, seq, quality_score):
        file.write('%s\n' % ID)
        file.write('%s\n' % seq)
        file.write('+\n')
        file.write('%s\n' % quality_score)
    def Parse_Fastq(self,CB1,CB2,Read):
        Allfiles=[CB1,CB2,Read]
        Allfiles = [gzip.open(files) for files in Allfiles]
        while True:
            try:
                names = [next(read).decode('UTF-8').split(' ')[0].strip('\n') for read in Allfiles]
                Sequence = [next(read).decode('UTF-8').strip('\n') for read in Allfiles]
                Blank = [next(read).decode('UTF-8').strip('\n') for read in Allfiles]
                qualityscore = [next(read).decode('UTF-8').strip('\n') for read in Allfiles]
                assert all(name.strip('\n') == names[0].strip('\n') for name in names)
                if names:
                    try:
                        yield [names[0].split(' ')[0], Sequence, qualityscore]
                    except:
                        return
                else:
                    break
            except StopIteration:
               break
        for read in Allfiles:
            read.close()
    def Parse_filtered_fastq(self,filtered_fastq):
        fastq=[filtered_fastq]
        fastq = [open(files) for files in fastq]
        while True:
            try:
                names = [next(read).strip('\n') for read in fastq]
                Sequence = [next(read).strip('\n') for read in fastq]
                Blank = [next(read).strip('\n') for read in fastq]
                qualityscore = [next(read).strip('\n') for read in fastq]
                assert all(name.strip('\n') == names[0].strip('\n') for name in names)
                if names:
                    try:
                        yield [names[0], Sequence, qualityscore]
                    except:
                        return
                else:
                    break
            except StopIteration:
               break
        for read in fastq:
            read.close()

    def __init__(self,sample_sheet,outputdir,Overall_sample_name,distance):
        self.output_dir=outputdir
        sample_info=pd.read_csv(sample_sheet,sep=',')
        if 'Read1' in list(sample_info.columns) and 'Read2' in list(sample_info.columns):
            self.datatype='Raw'
        else:
            self.datatype='Filtered'
        sample_set=list(sample_info['Individual_Sample_Name'])
        self.distance=distance
        self.Overall_sample_name=Overall_sample_name
        for sample in sample_set:
            if self.datatype=='Raw':
                R1=list(sample_info.loc[sample_info['Individual_Sample_Name']==sample,"Read1"])[0]
                R2=list(sample_info.loc[sample_info['Individual_Sample_Name']==sample,"Read2"])[0]
                R3=list(sample_info.loc[sample_info['Individual_Sample_Name']==sample,"Read3"])[0]
                self.Individual_sample[sample]=[R1,R2,R3]
            else:
                R3=list(sample_info.loc[sample_info['Individual_Sample_Name']==sample,"Read3"])[0]
                self.Individual_sample[sample]=[R3]

    def Extract_Filtering(self):
        outputfile = open('%s/%s_filtered.fastq' % (self.output_dir, self.Overall_sample_name), 'wt')
        if self.datatype=='Raw':
            for sample in self.Individual_sample.keys():
                R1=self.Individual_sample[sample][0]
                R2=self.Individual_sample[sample][1]
                R3=self.Individual_sample[sample][2]
                for read in self.Parse_Fastq(R1,R2,R3):
                    name = read[0]
                    CellBarcode1 = read[1][0]
                    CellBarcode2 = read[1][1][0:8]
                    UMI = read[1][1][8::]
                    Rabid_sequence = read[1][2]
                    Rabid_sequence_quality = read[2][2]
                    cellname = sample+'_'+CellBarcode1 + CellBarcode2
                    name = name + ' ' +sample+':'+ CellBarcode1 + ':' + CellBarcode2 + ':' + UMI
                    if cellname in self.cell_information.keys():
                        self.cell_information[cellname][0] += 1
                        if self.fiveendhandle in Rabid_sequence and self.threeendhandle in Rabid_sequence:
                            self.cell_information[cellname][3] += 1
                            realRabieread = Rabid_sequence[
                                            Rabid_sequence.find(self.fiveendhandle) + len(
                                                self.fiveendhandle):Rabid_sequence.find(
                                                self.threeendhandle)]
                            realQualityScore = Rabid_sequence_quality[
                                               Rabid_sequence.find(self.fiveendhandle) + len(
                                                   self.fiveendhandle):Rabid_sequence.find(
                                                   self.threeendhandle)]
                            if re.match(self.pattern, realRabieread):
                                self.cell_information[cellname][4] += 1
                                self._write_fastq(outputfile, name, realRabieread, realQualityScore)
                        else:
                            if self.fiveendhandle in Rabid_sequence:
                                self.cell_information[cellname][1] += 1
                                realRabieread = Rabid_sequence[
                                                Rabid_sequence.find(self.fiveendhandle) + len(
                                                    self.fiveendhandle):Rabid_sequence.find(
                                                    self.fiveendhandle) + len(self.fiveendhandle) + 28]
                                realQualityScore = Rabid_sequence_quality[
                                                   Rabid_sequence.find(self.fiveendhandle) + len(
                                                       self.fiveendhandle):Rabid_sequence.find(
                                                       self.fiveendhandle) + len(self.fiveendhandle) + 28]
                                if re.match(self.pattern, realRabieread):
                                    self.cell_information[cellname][4] += 1
                                    self._write_fastq(outputfile, name, realRabieread, realQualityScore)
                            elif self.threeendhandle in Rabid_sequence:
                                self.cell_information[cellname][2] += 1
                                realRabieread = Rabid_sequence[
                                                Rabid_sequence.find(self.threeendhandle) - 28:Rabid_sequence.find(
                                                    self.threeendhandle)]
                                realQualityScore = Rabid_sequence_quality[
                                                   Rabid_sequence.find(self.threeendhandle) - 28:Rabid_sequence.find(
                                                       self.threeendhandle)]
                                if re.match(self.pattern, realRabieread):
                                    self.cell_information[cellname][4] += 1
                                    self._write_fastq(outputfile, name, realRabieread, realQualityScore)
                    else:
                        self.cell_information[cellname] = [0, 0, 0, 0, 0]
                        self.cell_information[cellname][0] += 1
                        if self.fiveendhandle in Rabid_sequence and self.threeendhandle in Rabid_sequence:
                            self.cell_information[cellname][3] += 1
                            realRabieread = Rabid_sequence[Rabid_sequence.find(self.fiveendhandle) + len(
                                self.fiveendhandle):Rabid_sequence.find(self.threeendhandle)]
                            realQualityScore = Rabid_sequence_quality[Rabid_sequence.find(self.fiveendhandle) + len(
                                self.fiveendhandle):Rabid_sequence.find(self.threeendhandle)]
                            if re.match(self.pattern, realRabieread):
                                self.cell_information[cellname][4] += 1
                                self._write_fastq(outputfile, name, realRabieread, realQualityScore)
                        else:
                            if self.fiveendhandle in Rabid_sequence:
                                self.cell_information[cellname][1] += 1
                                realRabieread = Rabid_sequence[Rabid_sequence.find(self.fiveendhandle) + len(
                                    self.fiveendhandle):Rabid_sequence.find(self.fiveendhandle) + len(
                                    self.fiveendhandle) + 28]
                                realQualityScore = Rabid_sequence_quality[Rabid_sequence.find(self.fiveendhandle) + len(
                                    self.fiveendhandle):Rabid_sequence.find(self.fiveendhandle) + len(
                                    self.fiveendhandle) + 28]
                                if re.match(self.pattern, realRabieread):
                                    self.cell_information[cellname][4] += 1
                                    self._write_fastq(outputfile, name, realRabieread, realQualityScore)
                            elif self.threeendhandle in Rabid_sequence:
                                self.cell_information[cellname][2] += 1
                                realRabieread = Rabid_sequence[
                                                Rabid_sequence.find(self.threeendhandle) - 28:Rabid_sequence.find(
                                                    self.threeendhandle)]
                                realQualityScore = Rabid_sequence_quality[
                                                   Rabid_sequence.find(self.threeendhandle) - 28:Rabid_sequence.find(
                                                       self.threeendhandle)]
                                if re.match(self.pattern, realRabieread):
                                    self.cell_information[cellname][4] += 1
                                    self._write_fastq(outputfile, name, realRabieread, realQualityScore)
        else:
            for sample in self.Individual_sample.keys():
                R3=self.Individual_sample[sample][0]
                for read in self.Parse_filtered_fastq(R3):
                    name = read[0]
                    CellBarcode1 = read[0].split(sep=' ')[1].split(sep=':')[0]
                    CellBarcode2 = read[0].split(sep=' ')[1].split(sep=':')[1]
                    UMI = read[0].split(sep=' ')[1].split(sep=':')[2]
                    Rabid_sequence = read[1][0]
                    Rabid_sequence_quality = read[2][0]
                    cellname = sample+'_'+CellBarcode1 + CellBarcode2
                    name = read[0].split(sep=' ')[0] + ' ' +sample+':'+ CellBarcode1 + ':' + CellBarcode2 + ':' + UMI
                    if cellname in self.cell_information.keys():
                        self.cell_information[cellname][0] += 1
                        if self.fiveendhandle in Rabid_sequence and self.threeendhandle in Rabid_sequence:
                            self.cell_information[cellname][3] += 1
                            realRabieread = Rabid_sequence[
                                            Rabid_sequence.find(self.fiveendhandle) + len(
                                                self.fiveendhandle):Rabid_sequence.find(
                                                self.threeendhandle)]
                            realQualityScore = Rabid_sequence_quality[
                                               Rabid_sequence.find(self.fiveendhandle) + len(
                                                   self.fiveendhandle):Rabid_sequence.find(
                                                   self.threeendhandle)]
                            if re.match(self.pattern, realRabieread):
                                self.cell_information[cellname][4] += 1
                                self._write_fastq(outputfile, name, realRabieread, realQualityScore)
                        else:
                            if self.fiveendhandle in Rabid_sequence:
                                self.cell_information[cellname][1] += 1
                                realRabieread = Rabid_sequence[
                                                Rabid_sequence.find(self.fiveendhandle) + len(
                                                    self.fiveendhandle):Rabid_sequence.find(
                                                    self.fiveendhandle) + len(self.fiveendhandle) + 28]
                                realQualityScore = Rabid_sequence_quality[
                                                   Rabid_sequence.find(self.fiveendhandle) + len(
                                                       self.fiveendhandle):Rabid_sequence.find(
                                                       self.fiveendhandle) + len(self.fiveendhandle) + 28]
                                if re.match(self.pattern, realRabieread):
                                    self.cell_information[cellname][4] += 1
                                    self._write_fastq(outputfile, name, realRabieread, realQualityScore)
                            elif self.threeendhandle in Rabid_sequence:
                                self.cell_information[cellname][2] += 1
                                realRabieread = Rabid_sequence[
                                                Rabid_sequence.find(self.threeendhandle) - 28:Rabid_sequence.find(
                                                    self.threeendhandle)]
                                realQualityScore = Rabid_sequence_quality[
                                                   Rabid_sequence.find(self.threeendhandle) - 28:Rabid_sequence.find(
                                                       self.threeendhandle)]
                                if re.match(self.pattern, realRabieread):
                                    self.cell_information[cellname][4] += 1
                                    self._write_fastq(outputfile, name, realRabieread, realQualityScore)
                    else:
                        self.cell_information[cellname] = [0, 0, 0, 0, 0]
                        self.cell_information[cellname][0] += 1
                        if self.fiveendhandle in Rabid_sequence and self.threeendhandle in Rabid_sequence:
                            self.cell_information[cellname][3] += 1
                            realRabieread = Rabid_sequence[Rabid_sequence.find(self.fiveendhandle) + len(
                                self.fiveendhandle):Rabid_sequence.find(self.threeendhandle)]
                            realQualityScore = Rabid_sequence_quality[Rabid_sequence.find(self.fiveendhandle) + len(
                                self.fiveendhandle):Rabid_sequence.find(self.threeendhandle)]
                            if re.match(self.pattern, realRabieread):
                                self.cell_information[cellname][4] += 1
                                self._write_fastq(outputfile, name, realRabieread, realQualityScore)
                        else:
                            if self.fiveendhandle in Rabid_sequence:
                                self.cell_information[cellname][1] += 1
                                realRabieread = Rabid_sequence[Rabid_sequence.find(self.fiveendhandle) + len(
                                    self.fiveendhandle):Rabid_sequence.find(self.fiveendhandle) + len(
                                    self.fiveendhandle) + 28]
                                realQualityScore = Rabid_sequence_quality[Rabid_sequence.find(self.fiveendhandle) + len(
                                    self.fiveendhandle):Rabid_sequence.find(self.fiveendhandle) + len(
                                    self.fiveendhandle) + 28]
                                if re.match(self.pattern, realRabieread):
                                    self.cell_information[cellname][4] += 1
                                    self._write_fastq(outputfile, name, realRabieread, realQualityScore)
                            elif self.threeendhandle in Rabid_sequence:
                                self.cell_information[cellname][2] += 1
                                realRabieread = Rabid_sequence[
                                                Rabid_sequence.find(self.threeendhandle) - 28:Rabid_sequence.find(
                                                    self.threeendhandle)]
                                realQualityScore = Rabid_sequence_quality[
                                                   Rabid_sequence.find(self.threeendhandle) - 28:Rabid_sequence.find(
                                                       self.threeendhandle)]
                                if re.match(self.pattern, realRabieread):
                                    self.cell_information[cellname][4] += 1
                                    self._write_fastq(outputfile, name, realRabieread, realQualityScore)

        outputfile.close()
        with open('%s/%s_Cell_statistics.tsv' % (self.output_dir, self.Overall_sample_name), "w", newline="") as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow(
                ['Cellname', 'number of reads', 'number of reads with 5end handle', 'number of reads with 3end handle',
                 'number of reads with both handle',
                 'number of reads pass the structure filter'])  # number of reads, number of reads with both handle, number of reads pass the structure filter
            for cell in self.cell_information.keys():
                writer.writerow([self.Overall_sample_name+'_'+cell, self.cell_information[cell][0], self.cell_information[cell][1],
                                 self.cell_information[cell][2], self.cell_information[cell][3],
                                 self.cell_information[cell][4]])
        csvfile.close()
    def Clustering(self):
        distance=self.distance
        os.system(
            'starcode --print-clusters --dist %s -i %s/%s_filtered.fastq -o %s/%s.clustering.results' % (
            str(distance),self.output_dir, self.Overall_sample_name,self.output_dir, self.Overall_sample_name))
    def Correction_and_generate_table(self):
        if os.path.isfile('%s/%s.clustering.results'%(self.output_dir, self.Overall_sample_name)) is False:
            sys.exit('RabidSeq pipeline exiting, the Clustering results from STARCODE is missing.')
        starcode=open('%s/%s.clustering.results'%(self.output_dir, self.Overall_sample_name),'r')
        starcode.line=starcode.readlines()
        Rabid_correction_database={}
        for line in starcode.line:
            Rabid_Cluster_Center=line.split(sep='\t')[0]
            Rabid_Frequency=int(line.split(sep='\t')[1])
            Rabid_Cluster=line.split(sep='\t')[2]
            if Rabid_Frequency > 1:
                for member in Rabid_Cluster.split(sep=','):
                    Rabid_correction_database[member.strip('\n')]=Rabid_Cluster_Center.strip('\n')
        umi_list={}
        valid_cells=[cell for cell in list(self.cell_information.keys()) if self.cell_information[cell][4]>0]
        starcode.close()
        Rabid=Rabid_correction_database.values()
        Rabid=list(dict.fromkeys(Rabid))
        for rabie in Rabid:
            umi_list[rabie]=[]
        file_matrix=pd.DataFrame(0,index=list(Rabid),columns=list(valid_cells))

        for read in self.Parse_filtered_fastq('%s/%s_filtered.fastq'%(self.output_dir,self.Overall_sample_name)):
            CELLNAME=read[0].split(sep=' ')[1].split(sep=':')[0]+'_'+read[0].split(sep=' ')[1].split(sep=':')[1]+read[0].split(sep=' ')[1].split(sep=':')[2]
            UMI=read[0].split(sep=' ')[1].split(sep=':')[2]
            UMI=CELLNAME+UMI
            reads=read[1][0].strip('\n')
            if CELLNAME in valid_cells:
                if reads in Rabid_correction_database.keys():
                    if UMI not in umi_list[Rabid_correction_database[reads]]:
                        file_matrix.loc[Rabid_correction_database[reads],CELLNAME]+=1
                        umi_list[Rabid_correction_database[reads]].append(UMI)
        with open('%s/%s.table.csv'%(self.output_dir,self.Overall_sample_name), "w", newline="") as csvfile:
            writer=csv.writer(csvfile,delimiter='\t')
            writer.writerow(['Cellname','Rabid','Counts'])
            for i in range(len(valid_cells)):
                for j in range(len(Rabid)):
                    if file_matrix.loc[Rabid[j],valid_cells[i]]!=0:
                        writer.writerow([self.Overall_sample_name+'_'+valid_cells[i],Rabid[j],file_matrix.loc[Rabid[j],valid_cells[i]]])
        csvfile.close()





class Rabid_Barcodes_QC:
    #This class only checks diversity of Rabid Barcodes.
    Read=''
    samplename = ''
    outputdirectory = ''
    fiveendhandle = 'GCTAGC'
    threeendhandle = 'GGCGCGCC'
    pattern = '[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]AT[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]AT[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]'
    cell_information = [0, 0, 0, 0, 0]

    def __init__(self,Read,samplename,outputdirectory):
        self.Read=Read
        self.samplename=samplename
        self.outputdirectory=outputdirectory
    def processing(self):
        outputfile=open('%s/%s_filtered.fastq'%(self.outputdirectory,self.samplename),'wt')
        for read in ParseFastq([self.Read]):
            readid=read[0].strip('\n')
            Rabieread=read[1][0].strip('\n')
            QualityScore=read[2][0].strip('\n')
            self.cell_information[0]+=1
            if self.fiveendhandle in Rabieread and self.threeendhandle in Rabieread:
                self.cell_information[3]+=1
                realRabieread=Rabieread[Rabieread.find(self.fiveendhandle)+len(self.fiveendhandle):Rabieread.find(self.threeendhandle)]
                realQualityScore=QualityScore[Rabieread.find(self.fiveendhandle)+len(self.fiveendhandle):Rabieread.find(self.threeendhandle)]
                if re.match(self.pattern,realRabieread):
                    self.cell_information[4]+=1
                    write_fastq(outputfile,readid.strip('\n'),realRabieread,realQualityScore)
            else:
                if self.fiveendhandle in Rabieread:
                    self.cell_information[1]+=1
                    realRabieread=Rabieread[Rabieread.find(self.fiveendhandle)+len(self.fiveendhandle):Rabieread.find(self.fiveendhandle)+len(self.fiveendhandle)+28]
                    realQualityScore=QualityScore[Rabieread.find(self.fiveendhandle)+len(self.fiveendhandle):Rabieread.find(self.fiveendhandle)+len(self.fiveendhandle)+28]
                    if re.match(self.pattern,realRabieread):
                        self.cell_information[4]+=1
                        write_fastq(outputfile,readid.strip('\n'),realRabieread,realQualityScore)
                elif self.threeendhandle in Rabieread:
                    self.cell_information[2]+=1
                    realRabieread=Rabieread[Rabieread.find(self.threeendhandle)-28:Rabieread.find(self.threeendhandle)]
                    realQualityScore=QualityScore[Rabieread.find(self.threeendhandle)-28:Rabieread.find(self.threeendhandle)]
                    if re.match(self.pattern,realRabieread):
                        self.cell_information[4]+=1
                        write_fastq(outputfile,read[0].strip('\n'),realRabieread,realQualityScore)
        outputfile.close()

        os.system('starcode --print-clusters --dist 0 -i %s/%s_filtered.fastq -o %s/%s.clustering.results'%(self.outputdirectory,self.samplename,self.outputdirectory,self.samplename))

        os.system('gzip %s/%s_filtered.fastq'%(self.outputdirectory,self.samplename))
        with open('%s/%s.statistics.tsv'%(self.outputdirectory,self.samplename), "w", newline="") as csvfile:
            writer=csv.writer(csvfile,delimiter='\t')
            writer.writerow(['number of reads','Number of reads with only 5end handle','Number of reads with only 3end hadnle','Number of reads with both handles','Number of reads pass structure filter'])
            writer.writerow(self.cell_information)
        csvfile.close()

    def rarefactioncurve(self):
        starcode=open('%s/%s.clustering.results'%(self.outputdirectory,self.samplename),'r')
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
        1000000:[],5000000:[],
        10000000:[],50000000:[],
        100000000:[],500000000:[]}
        rare_faction_statistics={1000:[0,0,0,0,0],5000:[0,0,0,0,0],
        10000:[0,0,0,0,0],50000:[0,0,0,0,0],
        100000:[0,0,0,0,0],500000:[0,0,0,0,0],
        1000000:[0,0,0,0,0],5000000:[0,0,0,0,0],
        10000000:[0,0,0,0,0],50000000:[0,0,0,0,0],
        100000000:[0,0,0,0,0],500000000:[0,0,0,0,0]}#Total number of reads  number of reads number of reads with 5end handle    number of reads with 3end handle    number of reads with both handle    number of reads pass the structure filter
        for read in ParseFastq([self.Read]):
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


            if rare_faction_statistics[500000000][0] > 500000000:
                break

        with open('%s/%s.rarefactioncurve.tsv'%(self.outputdirectory,self.samplename), "w", newline="") as csvfile:
            writer=csv.writer(csvfile,delimiter='\t')
            writer.writerow(['Number of reads','Unique sequences'])
            for sampling in rare_faction:
                writer.writerow([sampling,len(list(dict.fromkeys(rare_faction[sampling])))])
        csvfile.close()
        with open('%s/%s.rarefactioncurve.statistics.tsv'%(self.outputdirectory,self.samplename), "w", newline="") as csvfile:
            writer=csv.writer(csvfile,delimiter='\t')
            writer.writerow(['Total number of reads','number of reads with 5end handle','number of reads with 3end handle','number of reads with both handle','number of reads pass the structure filter'])
            for sampling in rare_faction_statistics:
                writer.writerow([rare_faction_statistics[sampling][0],
                rare_faction_statistics[sampling][1],
                rare_faction_statistics[sampling][2],
                rare_faction_statistics[sampling][3],
                rare_faction_statistics[sampling][4]])
        csvfile.close()


if __name__=="__main__":
    options, arg = parser.parse_known_args()
    if options.quantify_from_inDrop_fastq_files:
        if os.path.isdir(options.outputdirectory) is False:
            sys.exit('RabidSeq pipeline exiting, the output directory does not exist')
        if os.path.isfile(options.Cellbarcode1) is False:
            sys.exit('RabidSeq pipeline exiting, the CB1 fastq file does not exist')
        if os.path.isfile(options.Cellbarcode2andUMI) is False:
            sys.exit('RabidSeq pipeline exiting, the CB2+UMI fastq file does not exist')
        if os.path.isfile(options.Read) is False:
            sys.exit('RabidSeq pipeline exiting, the Read fastq file does not exist')
        if options.Cellbarcode1.endswith('.gz') and options.Cellbarcode2andUMI.endswith('.gz') and options.Read.endswith('.gz') is False:
            sys.exit('RabidSeq pipeline exiting, the fastq files needs to be gzipped (Update supporting other format coming soon).')
        process=Rabid_Seq_Processor(CB=options.Cellbarcode1,
                                    Read=options.Read,
                                    samplename=options.name,
                                    outputdirectory=options.outputdirectory,
                                    distance=options.distance)
        process.Extract_and_Filtering()
        process.Clustering()
        process.Correction_and_generate_table()
    if options.quantify_from_filtered_fastq_files:
        if os.path.isdir(options.outputdirectory) is False:
            sys.exit('RabidSeq pipeline exiting, the output directory does not exist')
        if os.path.isfile(options.Read) is False:
            sys.exit('RabidSeq pipeline exiting, the Read fastq file does not exist')
        process=Rabid_Seq_Processor(CB='',
                                    Read=options.Read,
                                    samplename=options.name,
                                    outputdirectory=options.outputdirectory,
                                    distance=options.distance)
        process.Filtering_from_filtered_fastq()
        process.Clustering()
        process.Correction_and_generate_table()
    if options.quantify_from_multiple_samples:
        if os.path.isdir(options.outputdirectory) is False:
            sys.exit('RabidSeq pipeline exiting, the output directory does not exist')
        if os.path.isfile(options.sheet) is False:
            sys.exit('RabidSeq pipeline exiting, the sample sheet file does not exist')
        process=multi_Rabid_Seq_Processer(sample_sheet=options.sheet,
                                          outputdir=options.outputdirectory,
                                          Overall_sample_name=options.name,
                                          distance=options.distance)
        process.Extract_Filtering()
        process.Clustering()
        process.Correction_and_generate_table()
    if options.check_library_diversity:
        check_library_diversity=Rabid_Barcodes_QC(Read=options.Read,
                                                  samplename=options.name,
                                                  outputdirectory=options.outputdirectory)
        check_library_diversity.processing()
        check_library_diversity.rarefactioncurve()

