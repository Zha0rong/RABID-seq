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

from processors.rabid_seq_processor import RabidSeqProcessor

parser = argparse.ArgumentParser(prog="Rabid Seq pipeline", usage="RNAseq pipeline.")

parser.add_argument("--quantify_from_inDrop_raw_fastq_files", dest='quantify_from_inDrop_fastq_files',
                    action="store_true", help='Quantify from 3 inDrop files')
parser.add_argument("--quantify_from_inDrop_demultiplexed_fastq_files", dest='quantify_from_filtered_fastq_files',
                    action="store_true", help='Quantify from 1 inDrop fastq files')
parser.add_argument("--quantify_from_multiple_samples", dest='quantify_from_multiple_samples', action="store_true",
                    help='Quantify from multiple samples.')
parser.add_argument("--check_library_diversity", dest='check_library_diversity', action="store_true",
                    help='Check Rabid barcodes diversity in Plasmid Library.')
parser.add_argument('-R1', '--Cellbarcode1', dest='Cellbarcode1', action="store",
                    required="--quantify_from_inDrop_fastq_files" in sys.argv)
parser.add_argument('-R2', '--Cellbarcode2andUMI', dest='Cellbarcode2andUMI', action="store",
                    required="--quantify_from_inDrop_fastq_files" in sys.argv)
parser.add_argument('-R3', '--Read', dest='Read', action="store",
                    required="--quantify_from_inDrop_fastq_files" in sys.argv
                             or "--quantify_from_filtered_fastq_files" in sys.argv or "--check_library_diversity" in sys.argv)
parser.add_argument('-o', '--output', dest='outputdirectory', action="store",
                    required="--quantify_from_inDrop_fastq_files" in sys.argv
                             or "--quantify_from_filtered_fastq_files" in sys.argv or "--check_library_diversity" in sys.argv
                             or "--quantify_from_multiple_samples" in sys.argv)
parser.add_argument('-n', '--name', dest='name', action="store",
                    required="--quantify_from_inDrop_fastq_files" in sys.argv
                             or "--quantify_from_filtered_fastq_files" in sys.argv
                             or "--quantify_from_multiple_samples" in sys.argv or "--check_library_diversity" in sys.argv)
parser.add_argument('-s', '--sheet', dest='sheet', action="store",
                    required="--quantify_from_multiple_samples" in sys.argv)

parser.add_argument('-l', '--levenshtein', dest='distance', action='store', default=1)


def parse_fastq(paths_to_fastqs):

    if paths_to_fastqs[0].endswith('.gz'):
        processes = [subprocess.Popen(['zcat', (fastq)], stdout=subprocess.PIPE) for fastq in paths_to_fastqs]
        total_reads = [r.stdout for r in processes]

    elif paths_to_fastqs[0].endswith('.bz2'):
        processes = [subprocess.Popen(['bzcat', (fastq)], stdout=subprocess.PIPE) for fastq in paths_to_fastqs]
        total_reads = [r.stdout for r in processes]

    elif paths_to_fastqs[0].endswith('.fastq'):
        processes = [subprocess.Popen(['cat', (fastq)], stdout=subprocess.PIPE) for fastq in paths_to_fastqs]
        total_reads = [r.stdout for r in processes]

    else:
        sys.exit(f'The format of the file {(str(paths_to_fastqs))} is not recognized.')

    while True:
        names = [next(read).decode() for read in total_reads]
        sequence = [next(read).decode() for read in total_reads]
        blank = [next(read).decode() for read in total_reads]
        quality_score = [next(read).decode() for read in total_reads]

        assert all(name == names[0] for name in names)

        if names:
            yield [names[0], sequence, quality_score]
        else:
            break

    for read in total_reads:
        read.close()


def write_fastq(file, id, seq, quality_score):

    file.write('%s\n' % id)
    file.write('%s\n' % seq)
    file.write('+\n')
    file.write('%s\n' % quality_score)


if __name__ == "__main__":

    options, arg = parser.parse_known_args()

    if options.quantify_from_inDrop_fastq_files:

        if os.path.isdir(options.output_directory) is False:
            sys.exit('RabidSeq pipeline exiting, the output directory does not exist')

        if os.path.isfile(options.Cellbarcode1) is False:
            sys.exit('RabidSeq pipeline exiting, the CB1 fastq file does not exist')

        if os.path.isfile(options.Cellbarcode2andUMI) is False:
            sys.exit('RabidSeq pipeline exiting, the CB2+UMI fastq file does not exist')

        if os.path.isfile(options.read) is False:
            sys.exit('RabidSeq pipeline exiting, the Read fastq file does not exist')

        if options.Cellbarcode1.endswith('.gz') and options.Cellbarcode2andUMI.endswith(
                '.gz') and options.read.endswith('.gz') is False:
            sys.exit(
                'RabidSeq pipeline exiting, the fastq files needs to be gzipped (Update supporting other format coming soon).')

        process = RabidSeqProcessor(cb1=options.Cellbarcode1,
                                    cb2_umi=options.Cellbarcode2andUMI,
                                    read=options.read,
                                    sample_name=options.name,
                                    output_directory=options.output_directory,
                                    distance=options.distance)
        process.extract_and_filtering()
        process.clustering()
        process.correction_and_generate_table()

    if options.quantify_from_filtered_fastq_files:

        if os.path.isdir(options.output_directory) is False:
            sys.exit('RabidSeq pipeline exiting, the output directory does not exist')

        if os.path.isfile(options.read) is False:
            sys.exit('RabidSeq pipeline exiting, the Read fastq file does not exist')

        process = Rabid_Seq_Processor(CB1='',
                                      CB2_UMI='',
                                      Read=options.read,
                                      samplename=options.name,
                                      outputdirectory=options.output_directory,
                                      distance=options.distance)
        process.filtering_from_filtered_fastq()
        process.clustering()
        process.correction_and_generate_table()

    if options.quantify_from_multiple_samples:

        if os.path.isdir(options.output_directory) is False:
            sys.exit('RabidSeq pipeline exiting, the output directory does not exist')

        if os.path.isfile(options.sheet) is False:
            sys.exit('RabidSeq pipeline exiting, the sample sheet file does not exist')

        process = multi_Rabid_Seq_Processer(sample_sheet=options.sheet,
                                            outputdir=options.output_directory,
                                            Overall_sample_name=options.name,
                                            distance=options.distance)
        process.Extract_Filtering()
        process.clustering()
        process.correction_and_generate_table()

    if options.check_library_diversity:
        check_library_diversity = Rabid_Barcodes_QC(Read=options.read,
                                                    samplename=options.name,
                                                    outputdirectory=options.output_directory)
        check_library_diversity.processing()
        check_library_diversity.rarefaction_curve()
