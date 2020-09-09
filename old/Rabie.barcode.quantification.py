import os
import csv
import subprocess
import sys
from collections import defaultdict, OrderedDict
from itertools import product, combinations
import gzip
import numpy
import pandas as pd

def hammingdistance(string,reference):
    answer=0
    if (len(string)==len(reference)):
        for i in range(len(string)):
            if string[i]!=reference[i]:
                answer+=1
    else:
        answer=-1
    return answer


#Use Generator object to parse fastq file so the whole fastq file won't be load into the memory
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




if __name__=="__main__":
    namelist=os.listdir('./')
    namelist=[x for x in namelist if '_STAR_Reads.fastq.gz' in x]
    namelist=[x.replace('_STAR_Reads.fastq.gz','') for x in namelist]
    for i in range(len(namelist)):
        name=namelist[i]
        connection={}
        quantification={}
        umi={}
        cell=[]
        starcode=open('%s.clustering.results'%(name),'r')
        starcode.line=starcode.readlines()
        for line in starcode.line:
            connectionitem=line.split(sep='\t')[0]
            time=int(line.split(sep='\t')[1])
            cluster=line.split(sep='\t')[2]
            if time > 1 and len(connectionitem.strip('\n'))==28:
                for member in cluster.split(sep=','):
                    connection[member.strip('\n')]=connectionitem.strip('\n')
        starcode.close()
        value=connection.values()
        value=list(dict.fromkeys(value))
        for values in value:
            quantification[values]=0
            umi[values]=[]
        for cellname in ParseFastq(pathstofastqs=['%s.filtered.fastq'%(name)]):
        #print(cellname)
            cellname=cellname[0].split(' ')[1].split(':')
            cellname=cellname[0]+cellname[1]
            if cellname not in cell:
                cell.append(cellname)
        file_matrix=pd.DataFrame(0,index=list(value),columns=list(cell))
        for cellnames in ParseFastq(pathstofastqs=['%s.filtered.fastq'%(name)]):
            cellname=cellnames[0].split(' ')[1].split(':')
            cellname=cellname[0]+cellname[1]
            read=cellnames[1][0].strip('\n')
            umis=cellnames[0].split(' ')[1].split(':')[2]
            if read in connection.keys():
                if umis not in umi[connection[read]]:
                    file_matrix.loc[connection[read],cellname]+=1
                    umi[connection[read]].append(umis)
        file_matrix.to_csv('%s.connection.csv'%(name), sep=',', encoding='utf-8')