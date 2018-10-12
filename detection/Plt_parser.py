# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 18:58:04 2017

@author: xuedy
"""

import argparse

def gen_arg():
    parser = argparse.ArgumentParser(description='we devote to the exploration of CRISPR/Cas mechanism and optimization')
    parser.add_argument('-UID', metavar='userid', type=str,nargs='?',
                        dest='userid',
                        help='user id, necessary is -U is True, for the determination of data source')
    parser.add_argument('-D', metavar='description', type=str, nargs='?',
                        dest='description',
                        help='description of sgRNA')
    parser.add_argument('-U', metavar='upload', type=bool, nargs='?',choices=[True,False],
                        dest='upload',
                        help='if you are willing to upload part of your data and share them to all the researchers of CRIPSR/Cas')
    parser.add_argument('-O', metavar='output', type=str, nargs='?',
                        dest='output',
                        help='output folder')
    parser.add_argument('-L', metavar='label', type=str, nargs='?',
                        dest='label',
                        help='the name of CRISPR/Cas which you use in the experiment')
    parser.add_argument('-r1', metavar='read1', type=str, nargs='?',
                        dest='read1',
                        help='treated sequencing data 1, it need to be noticed that different methods require different types of input data')
    parser.add_argument('-r2', metavar='read2', type=str, nargs='?',
                        dest='read2',
                        help='treated sequencing data 2, if sequencing data is single-end, read1 is only necessary')
    subparsers =parser.add_subparsers(help = 'the experimental method of the data')
    parser_a = subparsers.add_parser('GUIDE-seq',help = 'genome-wide, unbiased identification of DSBs enabled by sequencing')
    parser_a.add_argument('-m', metavar='method', type=str, nargs='?',
                        dest='method',
                        default = 'GUIDE-seq',
                        help='DO NOT FILL ANY VALUE FOR THIS PARAMETER')
    parser_a.add_argument('-F', metavar='reference', type=str,nargs='?',
                        dest='reference',
                        help='/path/to/reference_genome.fa')
    parser_a.add_argument('-R', metavar='target', type=str, nargs='?',
                        dest='target',
                        help='target sequence including PAM')
    parser_a.add_argument('-bar1', metavar='barcode1', type=str, nargs='?',
                        dest='barcode1',
                        help='barcode 1, necessary in GUIDE-seq')
    parser_a.add_argument('-bar2', metavar='barcode2', type=str, nargs='?',
                        dest='barcode2',
                        help='barcode 2, necessary in GUIDE-seq')
    parser_a.add_argument('-cbar1', metavar='cbarcode1', type=str, nargs='?',
                        dest='cbarcode1',
                        help='control barcode 1, necessary in GUIDE-seq')
    parser_a.add_argument('-cbar2', metavar='cbarcode2', type=str, nargs='?',
                        dest='cbarcode2',
                        help='control barcode 2, necessary in GUIDE-seq')
    parser_a.add_argument('-ind1', metavar='index1', type=str, nargs='?',
                        dest='index1',
                        help='index 1, necessary in GUIDE-seq')
    parser_a.add_argument('-ind2', metavar='index2', type=str, nargs='?',
                        dest='index2',
                        help='index 2, necessary in GUIDE-seq')
    parser_a.add_argument('--d-minreads', metavar='demultiplex_min_reads', type=int, nargs='?',
                        dest='demultiplex_min_reads',
                        default=1000,
                        help='demultiplex_min_reads, parameters necessary in GUIDE-seq, default is 1000')
    parser_b = subparsers.add_parser('CIRCLE-seq',help = 'circularization for in vitro reporting of cleavage effects by sequencing')
    parser_b.add_argument('-m', metavar='method', type=str, nargs='?',
                        dest='method',
                        default = 'CIRCLE-seq',
                        help='DO NOT FILL ANY VALUE FOR THIS PARAMETER')
    parser_b.add_argument('-F', metavar='reference', type=str,nargs='?',
                        dest='reference',
                        help='/path/to/reference_genome.fa')
    parser_b.add_argument('-R', metavar='target', type=str, nargs='?',
                        dest='target',
                        help='target sequence including PAM')
    parser_b.add_argument('-cr1', metavar='cread1', type=str, nargs='?',
                        dest='cread1',
                        help='control sequencing data 1')
    parser_b.add_argument('-cr2', metavar='cread2', type=str, nargs='?',
                        dest='cread2',
                        help='control sequencing data 2')
    parser_b.add_argument('-rt', metavar='read_threshold', type=int, nargs='?',
                        dest='read_threshold',
                        default=4,
                        help='The minimum number of reads at a location for that location to be called as a site\nnecessary in CIRCLE-seq with default value 4')
    parser_b.add_argument('-ws', metavar='window_size', type=int, nargs='?',
                        dest='window_size',
                        default=3,
                        help='size of the sliding window\nnecessary in CIRCLE-seq with default value 3')
    parser_b.add_argument('-mqt', metavar='mapq_threshold', type=int, nargs='?',
                        dest='mapq_threshold',
                        default=50,
                        help='Minimum read mapping quality score\nnecessary in CIRCLE-seq with default value 50')
    parser_b.add_argument('-st', metavar='start_threshold', type=int, nargs='?',
                        dest='start_threshold',
                        default=1,
                        help='Tolerance for breakpoint location\nnecessary in CIRCLE-seq with default value 1')
    parser_b.add_argument('-gt', metavar='gap_threshold', type=int, nargs='?',
                        dest='gap_threshold',
                        default=3,
                        help='Number of tolerated gaps in the fuzzy target search setp\nnecessary in CIRCLE-seq with default value 3')
    parser_b.add_argument('-mt', metavar='mismatch_threshold', type=int, nargs='?',
                        dest='mismatch_threshold',
                        default=6,
                        help='Number of tolerated gaps in the fuzzy target search setp\nnecessary in CIRCLE-seq with default value 6')
    parser_b.add_argument('-ma', metavar='merged_analysis', type=bool, nargs='?',choices=[True,False],
                        dest='merged_analysis',
                        default=True,
                        help='Whether or not the paired read merging step should takingTrue\nnecessary in CIRCLE-seq with default value True')
    parser_c = subparsers.add_parser('SITE-seq',help = 'a biochemical method that uses the selective enrichment and identification of tagged genomic DNA ends by sequencing\n it should be noticed that the data provided by the paper has some problem and the data do not accord the mechanism of SITE-seq\nso we can not test the correctness of this part')
    parser_c.add_argument('-m', metavar='method', type=str, nargs='?',
                        dest='-method',
                        default = 'SITE-seq',
                        help='DO NOT FILL ANY VALUE FOR THIS PARAMETER')
    parser_c.add_argument('-R', metavar='target', type=str, nargs='?',
                        dest='target',
                        help='target sequence including PAM')
    parser_c.add_argument('-F', metavar='reference', type=str,nargs='?',
                        dest='reference',
                        help='/path/to/reference_genome.fa')
    """
    parser_d = subparsers.add_parser('DIGENOME-seq',help = 'in vitro Cas9-digested whole-genome sequencing, to profile genome-wide Cas9 off-target effects in human cells')
    parser_d.add_argument('-m', metavar='method', type=str, nargs='?',
                        dest='-method',
                        default = 'DIGENOME-seq',
                        help='DO NOT FILL ANY VALUE FOR THIS PARAMETER')
    parser_d.add_argument('-cr', metavar='cutoff_ratio', type=float, nargs='?',
                        dest='cutoff_ratio',
                        default = 20.0,
                        help='with default value 20.0')
    parser_d.add_argument('-cc', metavar='cutoff_count', type=int, nargs='?',
                        dest='cutoff_count',
                        default = 2,
                        help='with default value 2')
    parser_d.add_argument('-sr', metavar='sum_range', type=int, nargs='?',
                        dest='sum_range',
                        default = 1,
                        help='with default value 1')
    parser_d.add_argument('-s', metavar='step', type=int, nargs='?',
                        dest='step',
                        default = 1,
                        help='with default value 1')
    parser_d.add_argument('-ov', metavar='overhang', type=int, nargs='?',
                        dest='overhang',
                        default = 1,
                        help='e.g. SpCas9=1, ZFN=4, Cpf1=2, with default value 1')
    parser_e = subparsers.add_parser('BLESS',help = 'direct in situ breaks labeling, enrichment on streptavidin and next-generation sequencing')
    parser_e.add_argument('-m', metavar='method', type=str, nargs='?',
                        dest='-method',
                        default = 'BLESS',
                        help='DO NOT FILL ANY VALUE FOR THIS PARAMETER')
    parser_e.add_argument('-F', metavar='reference', type=str,nargs='?',
                        dest='reference',
                        help='absolute/path/to/reference_genome.fa')
    parser_e.add_argument('-C', metavar='control', type=str,nargs='?',
                        dest='reference',
                        help='absolute/path/to/reference_genome.fa')
    parser_e.add_argument('-t', metavar='threads', type=int, nargs='?',
                        dest='threads',
                        default = 1,
                        help='number of threads used during analysis, with default value 1')
    """
    
    return parser