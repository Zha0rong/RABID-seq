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


class RabidBarcodesQc:
    
    read = ''
    sample_name = ''
    output_directory = ''
    five_end_handle = 'GCTAGC'
    three_end_handle = 'GGCGCGCC'
    pattern = '[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]AT[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]AT[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]'
    cell_information = [0, 0, 0, 0, 0]


    def __init__(self, read, sample_name, output_directory):
        self.read = read
        self.sample_name = sample_name
        self.output_directory = output_directory
        
    def processing(self):
        
        with open(f'{self.output_directory}/{self.sample_name}_filtered.fastq', 'wt') as output_file:
        
            for read in parse_fastq([self.read]):
                read_id = read[0].strip('\n')
                rabie_read = read[1][0].strip('\n')
                quality_score = read[2][0].strip('\n')
                self.cell_information[0] += 1
                
                if self.five_end_handle in rabie_read and self.three_end_handle in rabie_read:
                    self.cell_information[3]+=1
                    real_rabie_read = rabie_read[rabie_read.find(self.five_end_handle) + len(self.five_end_handle):rabie_read.find(self.three_end_handle)]
                    real_quality_score = quality_score[rabie_read.find(self.five_end_handle) + len(self.five_end_handle):rabie_read.find(self.three_end_handle)]
                    
                    if re.match(self.pattern, real_rabie_read):
                        self.cell_information[4] += 1
                        write_fastq(output_file, read_id.strip('\n'), real_rabie_read, eal_quality_score)
                else:
                    if self.five_end_handle in rabie_read:
                        self.cell_information[1] += 1
                        real_rabie_read = rabie_read[rabie_read.find(self.five_end_handle) + len(self.five_end_handle):rabie_read.find(self.five_end_handle) + len(self.five_end_handle) + 28]
                        real_quality_score = quality_score[rabie_read.find(self.five_end_handle) + len(self.five_end_handle):rabie_read.find(self.five_end_handle) + len(self.five_end_handle) + 28]
                        
                        if re.match(self.pattern,real_rabie_read):
                            self.cell_information[4]+=1
                            write_fastq(output_file,read_id.strip('\n'), real_rabie_read, real_quality_score)
                            
                    elif self.three_end_handle in rabie_read:
                        self.cell_information[2]+=1
                        real_rabie_read= rabie_read[rabie_read.find(self.three_end_handle) - 28:rabie_read.find(self.three_end_handle)]
                        real_quality_score= quality_score[rabie_read.find(self.three_end_handle) - 28:rabie_read.find(self.three_end_handle)]
                        
                        if re.match(self.pattern,real_rabie_read):
                            self.cell_information[4]+=1
                            write_fastq(output_file, read[0].strip('\n'), real_rabie_read, real_quality_score)

        os.system(f'starcode --print-clusters --dist 0 -i {self.output_directory}/{self.sample_name}_filtered.fastq -o {self.output_directory}/{self.sample_name}.clustering.results')
        os.system(f'gzip {self.output_directory}/{self.sample_name}_filtered.fastq')
        
        with open(f'{self.output_directory}/{self.sample_name}.statistics.tsv', "w", newline="") as csv_file:
            writer=csv.writer(csv_file,delimiter='\t')
            writer.writerow(['number of reads',
                             'Number of reads with only 5end handle',
                             'Number of reads with only 3end hadnle','Number of reads with both handles',
                             'Number of reads pass structure filter'])
            writer.writerow(self.cell_information)

    def rarefaction_curve(self):
        
        with open(f'{self.output_directory}/{self.sample_name}.clustering.results', 'r') as star_code:
            star_code.line = star_code.readlines()
            true_barcode = {}
            
            for line in star_code.line:
                connection_item = line.split(sep='\t')[0]
                time = int(line.split(sep='\t')[1])
                cluster = line.split(sep='\t')[2]
                
                if time > 1:
                    for member in cluster.split(sep=','):
                        true_barcode[member.strip('\n')]=connection_item.strip('\n')
                        
            
        rare_faction = {
            1000:[],
            5000:[],
            10000:[],
            50000:[],
            100000:[],
            500000:[],
            1000000:[],
            5000000:[],
            10000000:[],
            50000000:[],
            100000000:[],
            500000000:[]
        }
        
        rare_faction_statistics = {
            1000:[0, 0, 0, 0, 0],
            5000:[0, 0, 0, 0, 0],
            10000:[0, 0, 0, 0, 0],
            50000:[0, 0, 0, 0, 0],
            100000:[0, 0, 0, 0, 0],
            500000:[0, 0, 0, 0, 0],
            1000000:[0, 0, 0, 0, 0],
            5000000:[0, 0, 0, 0, 0],
            10000000:[0, 0, 0, 0, 0],
            50000000:[0, 0, 0, 0, 0],
            100000000:[0, 0, 0, 0, 0],
            500000000:[0, 0, 0, 0, 0]
        }
        
        for read in parse_fastq([self.read]):
            read_id = read[0].strip('\n').split(' ')[0]
            rabie_read = read[1][0].strip('\n')

            for faction in rare_faction:
                if rare_faction_statistics[faction][0] < faction:
                    random_selector=random.randint(1,101)

                    if random_selector <= 50:
                        rare_faction_statistics[faction][0] += 1

                        if self.five_end_handle in rabie_read and self.three_end_handle in rabie_read:
                            real_rabie_read = rabie_read[rabie_read.find(self.five_end_handle) + len(self.five_end_handle):rabie_read.find(self.three_end_handle)]

                            if re.match(self.pattern,real_rabie_read):

                                if real_rabie_read in true_barcode.keys():
                                    rare_faction[faction].append(true_barcode[real_rabie_read])

                                rare_faction_statistics[faction][3] += 1
                                rare_faction_statistics[faction][4] += 1

                        else:
                            if self.five_end_handle in rabie_read:
                                real_rabie_read = rabie_read[rabie_read.find(self.five_end_handle) + len(self.five_end_handle):rabie_read.find(self.five_end_handle) + len(self.five_end_handle) + 28]

                                if re.match(self.pattern,real_rabie_read):

                                    if real_rabie_read in true_barcode.keys():
                                        rare_faction[faction].append(true_barcode[real_rabie_read])

                                    rare_faction_statistics[faction][1] += 1
                                    rare_faction_statistics[faction][4] += 1

                            elif self.three_end_handle in rabie_read:
                                real_rabie_read = rabie_read[rabie_read.find(self.three_end_handle) - 28:rabie_read.find(self.three_end_handle)]

                                if re.match(self.pattern,real_rabie_read):

                                    if real_rabie_read in true_barcode.keys():
                                        rare_faction[faction].append(true_barcode[real_rabie_read])

                                    rare_faction_statistics[faction][2] += 1

                                rare_faction_statistics[faction][4] += 1


            if rare_faction_statistics[500000000][0] > 500000000:
                break

        with open(f'{self.output_directory}/{self.sample_name}.rarefactioncurve.tsv', "w", newline="") as csv_file:
            writer = csv.writer(csv_file,delimiter='\t')
            writer.writerow(['Number of reads','Unique sequences'])

            for sampling in rare_faction:
                writer.writerow([sampling, len(list(dict.fromkeys(rare_faction[sampling])))])

        with open(f'{self.output_directory}/{self.sample_name}.rarefactioncurve.statistics.tsv', "w", newline="") as csv_file:
            writer = csv.writer(csv_file,delimiter='\t')
            writer.writerow(['Total number of reads',
                             'number of reads with 5end handle',
                             'number of reads with 3end handle',
                             'number of reads with both handle',
                             'number of reads pass the structure filter'])

            for sampling in rare_faction_statistics:
                writer.writerow([rare_faction_statistics[sampling][0],
                rare_faction_statistics[sampling][1],
                rare_faction_statistics[sampling][2],
                rare_faction_statistics[sampling][3],
                rare_faction_statistics[sampling][4]])
