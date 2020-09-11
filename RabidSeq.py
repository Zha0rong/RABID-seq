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

parser = argparse.ArgumentParser(prog = "Rabid Seq pipeline",usage="RNAseq pipeline.")



parser.add_argument("--quantify_from_inDrop_raw_fastq_files",dest='quantify_from_inDrop_fastq_files',action="store_true",help='Quantify from 3 inDrop files')
parser.add_argument("--quantify_from_inDrop_demultiplexed_fastq_files",dest='quantify_from_filtered_fastq_files',action="store_true",help='Quantify from 1 inDrop fastq files')
parser.add_argument("--quantify_from_multiple_samples",dest='quantify_from_multiple_samples',action="store_true",help='Quantify from multiple samples.')

parser.add_argument('-R1','--Cellbarcode1',dest='Cellbarcode1',action="store",required="--quantify_from_inDrop_fastq_files" in sys.argv)
parser.add_argument('-R2','--Cellbarcode2andUMI',dest='Cellbarcode2andUMI',action="store",required="--quantify_from_inDrop_fastq_files" in sys.argv)
parser.add_argument('-R3','--Read',dest='Read',action="store",required="--quantify_from_inDrop_fastq_files" in sys.argv
                                                                       or "--quantify_from_filtered_fastq_files" in sys.argv)
parser.add_argument('-o','--output',dest='outputdirectory',action="store",required="--quantify_from_inDrop_fastq_files" in sys.argv
                                                                                   or "--quantify_from_filtered_fastq_files" in sys.argv
                    or "--quantify_from_multiple_samples" in sys.argv)
parser.add_argument('-n','--name',dest='name',action="store",required="--quantify_from_inDrop_fastq_files" in sys.argv
                                                                      or "--quantify_from_filtered_fastq_files" in sys.argv
                    or "--quantify_from_multiple_samples" in sys.argv)
parser.add_argument('-s','--sheet',dest='sheet',action="store",required="--quantify_from_multiple_samples" in sys.argv)

parser.add_argument('-l','--levenshtein',dest='distance',action='store',default=1)


class Rabid_Seq_Processor:
    '''This class is used to read in the inDrop fastq data (3 gzipped file: Cell Barcode 1, Cell Barcode 2 + UMI, Rabid (RNA) Reads).
       The output will be '''
    CB1=''
    CB2_UMI=''
    Read=''
    samplename=''
    outputdirectory=''
    fiveendhandle = 'GCTAGC'
    threeendhandle = 'GGCGCGCC'
    pattern = '[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]AT[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]AT[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]'
    cell_information = {}  # Cellname:[number of reads,number of reads with 5end handle,number of reads with 3end handle, number of reads with both handle, number of reads pass the structure filter]
    whitelist = {}
    distance=0
    def __init__(self,CB1,CB2_UMI,Read,samplename,outputdirectory,distance):
        self.CB1=CB1
        self.CB2_UMI=CB2_UMI
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
        Allfiles=[self.CB1,self.CB2_UMI,self.Read]
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
        outputfile = open('%s/%s_filtered.fastq' % (self.outputdirectory,self.samplename), 'wt')
        for read in self._ParseFastq():
            name=read[0]
            CellBarcode1=read[1][0]
            CellBarcode2=read[1][1][0:8]
            UMI=read[1][1][8::]
            Rabid_sequence=read[1][2]
            Rabid_sequence_quality=read[2][2]
            cellname=CellBarcode1+CellBarcode2
            name=name+' '+CellBarcode1+':'+CellBarcode2+':'+UMI
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
            CellBarcode1=read[0].split(sep=' ')[1].split(sep=':')[0]
            CellBarcode2=read[0].split(sep=' ')[1].split(sep=':')[1]
            Rabid_sequence=read[1][0]
            Rabid_sequence_quality=read[2][0]
            cellname=CellBarcode1+CellBarcode2
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
        valid_cells=[cell for cell in list(self.cell_information.keys()) if self.cell_information[cell][4]>0] # Only
        starcode.close()
        Rabid=Rabid_correction_database.values()
        Rabid=list(dict.fromkeys(Rabid))
        for rabie in Rabid:
            umi_list[rabie]=[]
        file_matrix=pd.DataFrame(0,index=list(Rabid),columns=list(valid_cells))

        for read in self._Parse_filtered_fastq():
            CELLNAME=read[0].split(sep=' ')[1].split(sep=':')[0]+read[0].split(sep=' ')[1].split(sep=':')[1]
            UMI=read[0].split(sep=' ')[1].split(sep=':')[2]
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
    def ParseFastq(self,CB1,CB2,Read):
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
            R1=list(sample_info.loc[sample_info['Individual_Sample_Name']==sample,"Read1"])[0]
            R2=list(sample_info.loc[sample_info['Individual_Sample_Name']==sample,"Read2"])[0]
            R3=list(sample_info.loc[sample_info['Individual_Sample_Name']==sample,"Read3"])[0]
            self.Individual_sample[sample]=[R1,R2,R3]
        print(self.Individual_sample)
        print(self.Overall_sample_name)

    def Extract_Filtering(self):
        return 0



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
        process=Rabid_Seq_Processor(CB1=options.Cellbarcode1,
                                    CB2_UMI=options.Cellbarcode2andUMI,
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
        process=Rabid_Seq_Processor(CB1='',
                                    CB2_UMI='',
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