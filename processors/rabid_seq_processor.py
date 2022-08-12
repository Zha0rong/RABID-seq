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

from models.rabid_barcodes_qc import RabidBarcodesQc
from processors.multi_rabid_seq_processor import MultiRabidSeqProcesser


class RabidSeqProcessor:
    '''This class is used to read in the inDrop fastq data (3 gzipped file: Cell Barcode 1, Cell Barcode 2 + UMI, Rabid (RNA) Reads).
       The output will be '''

    cb1 = ''
    cb2_umi = ''
    read = ''
    sample_name = ''
    output_directory = ''
    five_end_handle = 'GCTAGC'
    three_end_handle = 'GGCGCGCC'
    pattern = '[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]AT[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]AT[AGC][ACT][AGT][GCT][ACG][ACT][AGT][GCT]'
    cell_information = {}
    whitelist = {}
    distance = 0

    def __init__(self, cb1, cb2_umi, read, sample_name, output_directory, distance):

        self.cb1 = cb1
        self.cb2_umi = cb2_umi
        self.read = read
        self.sample_name = sample_name
        self.output_directory = output_directory
        self.distance = distance

    def _write_fastq(self, file, id, seq, quality_score):

        file.write('%s\n' % id)
        file.write('%s\n' % seq)
        file.write('+\n')
        file.write('%s\n' % quality_score)

    def _parse_fastq(self):

        all_files = [self.cb1, self.cb2_umi, self.read]
        all_files = [gzip.open(files) for files in all_files]

        while True:
            try:
                names = [next(read).decode('UTF-8').split(' ')[0].strip('\n') for read in all_files]
                sequence = [next(read).decode('UTF-8').strip('\n') for read in all_files]
                blank = [next(read).decode('UTF-8').strip('\n') for read in all_files]
                quality_score = [next(read).decode('UTF-8').strip('\n') for read in all_files]

                assert all(name.strip('\n') == names[0].strip('\n') for name in names)

                if names:
                    try:
                        yield [names[0].split(' ')[0], sequence, quality_score]
                    except:
                        return
                else:
                    break
            except StopIteration:
                break

        for read in all_files:
            read.close()

    def _parse_filtered_fastq(self, filtered_fastq=None):

        if filtered_fastq is None:
            fastq = ['%s/%s_filtered.fastq' % (self.output_directory, self.sample_name)]
        else:
            fastq = [filtered_fastq]

        fastq = [open(files) for files in fastq]

        while True:
            try:
                names = [next(read).strip('\n') for read in fastq]
                sequence = [next(read).strip('\n') for read in fastq]
                blank = [next(read).strip('\n') for read in fastq]
                quality_score = [next(read).strip('\n') for read in fastq]

                assert all(name.strip('\n') == names[0].strip('\n') for name in names)

                if names:
                    try:
                        yield [names[0], sequence, quality_score]
                    except:
                        return
                else:
                    break
            except StopIteration:
                break

        for read in fastq:
            read.close()

    def extract_and_filtering(self):

        with open(f'{self.output_directory}/{self.sample_name}_filtered.fastq', 'wt') as output_file:

            for read in self._parse_fastq():
                name = read[0]
                cell_barcode_1 = read[1][0]
                cell_barcode_2 = read[1][1][0:8]
                umi = read[1][1][8::]
                rabid_sequence = read[1][2]
                rabid_sequence_quality = read[2][2]
                cell_name = cell_barcode_1 + cell_barcode_2
                name = f'{name} {cell_barcode_1}:{cell_barcode_2}:{umi}'

                if cell_name in self.cell_information.keys():
                    self.cell_information[cell_name][0] += 1

                    if self.five_end_handle in rabid_sequence and self.three_end_handle in rabid_sequence:
                        self.cell_information[cell_name][3] += 1
                        real_rabie_read = rabid_sequence[
                                          rabid_sequence.find(self.five_end_handle) + len(
                                              self.five_end_handle):rabid_sequence.find(
                                              self.three_end_handle)]
                        real_quality_score = rabid_sequence_quality[
                                             rabid_sequence.find(self.five_end_handle) + len(
                                                 self.five_end_handle):rabid_sequence.find(
                                                 self.three_end_handle)]

                        if re.match(self.pattern, real_rabie_read):
                            self.cell_information[cell_name][4] += 1
                            self._write_fastq(output_file, name, real_rabie_read, real_quality_score)

                    else:
                        if self.five_end_handle in rabid_sequence:
                            self.cell_information[cell_name][1] += 1
                            real_rabie_read = rabid_sequence[
                                              rabid_sequence.find(self.five_end_handle) + len(
                                                  self.five_end_handle):rabid_sequence.find(
                                                  self.five_end_handle) + len(self.five_end_handle) + 28]
                            real_quality_score = rabid_sequence_quality[
                                                 rabid_sequence.find(self.five_end_handle) + len(
                                                     self.five_end_handle):rabid_sequence.find(
                                                     self.five_end_handle) + len(self.five_end_handle) + 28]

                            if re.match(self.pattern, real_rabie_read):
                                self.cell_information[cell_name][4] += 1
                                self._write_fastq(output_file, name, real_rabie_read, real_quality_score)

                        elif self.three_end_handle in rabid_sequence:
                            self.cell_information[cell_name][2] += 1
                            real_rabie_read = rabid_sequence[
                                              rabid_sequence.find(self.three_end_handle) - 28:rabid_sequence.find(
                                                  self.three_end_handle)]
                            real_quality_score = rabid_sequence_quality[
                                                 rabid_sequence.find(self.three_end_handle) - 28:rabid_sequence.find(
                                                     self.three_end_handle)]

                            if re.match(self.pattern, real_rabie_read):
                                self.cell_information[cell_name][4] += 1
                                self._write_fastq(output_file, name, real_rabie_read, real_quality_score)
                else:
                    self.cell_information[cell_name] = [0, 0, 0, 0, 0]
                    self.cell_information[cell_name][0] += 1

                    if self.five_end_handle in rabid_sequence and self.three_end_handle in rabid_sequence:
                        self.cell_information[cell_name][3] += 1
                        real_rabie_read = rabid_sequence[rabid_sequence.find(self.five_end_handle) + len(
                            self.five_end_handle):rabid_sequence.find(self.three_end_handle)]
                        real_quality_score = rabid_sequence_quality[rabid_sequence.find(self.five_end_handle) + len(
                            self.five_end_handle):rabid_sequence.find(self.three_end_handle)]
                        if re.match(self.pattern, real_rabie_read):
                            self.cell_information[cell_name][4] += 1
                            self._write_fastq(output_file, name, real_rabie_read, real_quality_score)

                    else:
                        if self.five_end_handle in rabid_sequence:
                            self.cell_information[cell_name][1] += 1
                            real_rabie_read = rabid_sequence[rabid_sequence.find(self.five_end_handle) + len(
                                self.five_end_handle):rabid_sequence.find(self.five_end_handle) + len(
                                self.five_end_handle) + 28]
                            real_quality_score = rabid_sequence_quality[rabid_sequence.find(self.five_end_handle) + len(
                                self.five_end_handle):rabid_sequence.find(self.five_end_handle) + len(
                                self.five_end_handle) + 28]

                            if re.match(self.pattern, real_rabie_read):
                                self.cell_information[cell_name][4] += 1
                                self._write_fastq(output_file, name, real_rabie_read, real_quality_score)

                        elif self.three_end_handle in rabid_sequence:
                            self.cell_information[cell_name][2] += 1
                            real_rabie_read = rabid_sequence[
                                              rabid_sequence.find(self.three_end_handle) - 28:rabid_sequence.find(
                                                  self.three_end_handle)]
                            real_quality_score = rabid_sequence_quality[
                                                 rabid_sequence.find(self.three_end_handle) - 28:rabid_sequence.find(
                                                     self.three_end_handle)]

                            if re.match(self.pattern, real_rabie_read):
                                self.cell_information[cell_name][4] += 1
                                self._write_fastq(output_file, name, real_rabie_read, real_quality_score)

        with open(f'{self.output_directory}/{self.sample_name}_Cell_statistics.tsv', "w", newline="") as csv_file:
            writer = csv.writer(csv_file, delimiter='\t')
            writer.writerow(['Cellname',
                             'number of reads',
                             'number of reads with 5end handle',
                             'number of reads with 3end handle',
                             'number of reads with both handle',
                             'number of reads pass the structure filter'])
            for cell in self.cell_information.keys():
                writer.writerow([cell, self.cell_information[cell][0],
                                 self.cell_information[cell][1],
                                 self.cell_information[cell][2],
                                 self.cell_information[cell][3],
                                 self.cell_information[cell][4]])

    def filtering_from_filtered_fastq(self):

        with open(f'{self.output_directory}/{self.sample_name}_filtered.fastq', 'wt') as output_file:
            for read in self._parse_filtered_fastq(filtered_fastq=self.read):
                name = read[0]
                cell_barcode_1 = read[0].split(sep=' ')[1].split(sep=':')[0]
                cell_barcode_2 = read[0].split(sep=' ')[1].split(sep=':')[1]
                rabid_sequence = read[1][0]
                rabid_sequence_quality = read[2][0]
                cell_name = cell_barcode_1 + cell_barcode_2

                if cell_name in self.cell_information.keys():
                    self.cell_information[cell_name][0] += 1

                    if self.five_end_handle in rabid_sequence and self.three_end_handle in rabid_sequence:
                        self.cell_information[cell_name][3] += 1
                        real_rabie_read = rabid_sequence[
                                          rabid_sequence.find(self.five_end_handle) + len(
                                              self.five_end_handle):rabid_sequence.find(
                                              self.three_end_handle)]
                        real_quality_score = rabid_sequence_quality[
                                             rabid_sequence.find(self.five_end_handle) + len(
                                                 self.five_end_handle):rabid_sequence.find(
                                                 self.three_end_handle)]

                        if re.match(self.pattern, real_rabie_read):
                            self.cell_information[cell_name][4] += 1
                            self._write_fastq(output_file, name, real_rabie_read, real_quality_score)
                    else:
                        if self.five_end_handle in rabid_sequence:
                            self.cell_information[cell_name][1] += 1
                            real_rabie_read = rabid_sequence[
                                              rabid_sequence.find(self.five_end_handle) + len(
                                                  self.five_end_handle):rabid_sequence.find(
                                                  self.five_end_handle) + len(self.five_end_handle) + 28]
                            real_quality_score = rabid_sequence_quality[
                                                 rabid_sequence.find(self.five_end_handle) + len(
                                                     self.five_end_handle):rabid_sequence.find(
                                                     self.five_end_handle) + len(self.five_end_handle) + 28]

                            if re.match(self.pattern, real_rabie_read):
                                self.cell_information[cell_name][4] += 1
                                self._write_fastq(output_file, name, real_rabie_read, real_quality_score)

                        elif self.three_end_handle in rabid_sequence:
                            self.cell_information[cell_name][2] += 1
                            real_rabie_read = rabid_sequence[
                                              rabid_sequence.find(self.three_end_handle) - 28:rabid_sequence.find(
                                                  self.three_end_handle)]
                            real_quality_score = rabid_sequence_quality[
                                                 rabid_sequence.find(self.three_end_handle) - 28:rabid_sequence.find(
                                                     self.three_end_handle)]

                            if re.match(self.pattern, real_rabie_read):
                                self.cell_information[cell_name][4] += 1
                                self._write_fastq(output_file, name, real_rabie_read, real_quality_score)
                else:
                    self.cell_information[cell_name] = [0, 0, 0, 0, 0]
                    self.cell_information[cell_name][0] += 1

                    if self.five_end_handle in rabid_sequence and self.three_end_handle in rabid_sequence:
                        self.cell_information[cell_name][3] += 1
                        real_rabie_read = rabid_sequence[rabid_sequence.find(self.five_end_handle) + len(
                            self.five_end_handle):rabid_sequence.find(self.three_end_handle)]
                        real_quality_score = rabid_sequence_quality[rabid_sequence.find(self.five_end_handle) + len(
                            self.five_end_handle):rabid_sequence.find(self.three_end_handle)]

                        if re.match(self.pattern, real_rabie_read):
                            self.cell_information[cell_name][4] += 1
                            self._write_fastq(output_file, name, real_rabie_read, real_quality_score)
                    else:
                        if self.five_end_handle in rabid_sequence:
                            self.cell_information[cell_name][1] += 1
                            real_rabie_read = rabid_sequence[rabid_sequence.find(self.five_end_handle) + len(
                                self.five_end_handle):rabid_sequence.find(self.five_end_handle) + len(
                                self.five_end_handle) + 28]
                            real_quality_score = rabid_sequence_quality[rabid_sequence.find(self.five_end_handle) + len(
                                self.five_end_handle):rabid_sequence.find(self.five_end_handle) + len(
                                self.five_end_handle) + 28]

                            if re.match(self.pattern, real_rabie_read):
                                self.cell_information[cell_name][4] += 1
                                self._write_fastq(output_file, name, real_rabie_read, real_quality_score)

                        elif self.three_end_handle in rabid_sequence:
                            self.cell_information[cell_name][2] += 1
                            real_rabie_read = rabid_sequence[
                                              rabid_sequence.find(self.three_end_handle) - 28:rabid_sequence.find(
                                                  self.three_end_handle)]
                            real_quality_score = rabid_sequence_quality[
                                                 rabid_sequence.find(self.three_end_handle) - 28:rabid_sequence.find(
                                                     self.three_end_handle)]

                            if re.match(self.pattern, real_rabie_read):
                                self.cell_information[cell_name][4] += 1
                                self._write_fastq(output_file, name, real_rabie_read, real_quality_score)

        with open(f'{self.output_directory}/{self.sample_name}_Cell_statistics.tsv', "w", newline="") as csv_file:
            writer = csv.writer(csv_file, delimiter='\t')
            writer.writerow(
                ['Cellname', 'number of reads', 'number of reads with 5end handle', 'number of reads with 3end handle',
                 'number of reads with both handle',
                 'number of reads pass the structure filter'])
            for cell in self.cell_information.keys():
                writer.writerow([cell, self.cell_information[cell][0], self.cell_information[cell][1],
                                 self.cell_information[cell][2], self.cell_information[cell][3],
                                 self.cell_information[cell][4]])

    def clustering(self):

        distance = self.distance
        os.system(
            'starcode --print-clusters --dist %s -i %s/%s_filtered.fastq -o %s/%s.clustering.results' % (
                str(distance), self.output_directory, self.sample_name, self.output_directory, self.sample_name))

    def correction_and_generate_table(self):

        if os.path.isfile('%s/%s.clustering.results' % (self.output_directory, self.sample_name)) is False:
            sys.exit('RabidSeq pipeline exiting, the Clustering results from STARCODE is missing.')

        star_code = open(f'{self.output_directory}/{self.sample_name}.clustering.results', 'r')
        star_code.line = star_code.readlines()
        rabid_correction_database = {}

        for line in star_code.line:
            rabid_cluster_center = line.split(sep='\t')[0]
            rabid_frequency = int(line.split(sep='\t')[1])
            rabid_cluster = line.split(sep='\t')[2]

            if rabid_frequency > 1:
                for member in rabid_cluster.split(sep=','):
                    rabid_correction_database[member.strip('\n')] = rabid_cluster_center.strip('\n')

        umi_list = {}
        valid_cells = [cell for cell in list(self.cell_information.keys()) if self.cell_information[cell][4] > 0]

        star_code.close()

        rabid = rabid_correction_database.values()
        rabid = list(dict.fromkeys(rabid))

        for rabie in rabid:
            umi_list[rabie] = []

        file_matrix = pd.DataFrame(0, index=list(rabid), columns=list(valid_cells))

        for read in self._parse_filtered_fastq():

            cell_name_line = read[0].split(sep=' ')[1].split(sep=':')[0] + read[0].split(sep=' ')[1].split(sep=':')[1]
            umi_line = read[0].split(sep=' ')[1].split(sep=':')[2]
            umi_line = cell_name_line + umi_line
            reads = read[1][0].strip('\n')

            if cell_name_line in valid_cells:
                if reads in rabid_correction_database.keys():
                    if umi_line not in umi_list[rabid_correction_database[reads]]:
                        file_matrix.loc[rabid_correction_database[reads], cell_name_line] += 1
                        umi_list[rabid_correction_database[reads]].append(umi_line)

        with open(f'{self.output_directory}/{self.sample_name}.table.csv', "w", newline="") as csv_file:
            writer = csv.writer(csv_file, delimiter='\t')
            writer.writerow(['Cellname', 'Rabid', ' Counts'])

            for i in range(len(valid_cells)):
                for j in range(len(rabid)):
                    if file_matrix.loc[rabid[j], valid_cells[i]] != 0:
                        writer.writerow([valid_cells[i], rabid[j], file_matrix.loc[rabid[j], valid_cells[i]]])