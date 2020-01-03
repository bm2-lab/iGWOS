# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 12:44:46 2017

@author: xuedy
"""

def yaml_structure(method):
    sourcei = {'C':{'target':'',
                    'read1':'',
                    'read2':'',
                    'controlread1':'',
                    'controlread2':'',
                    'description':''},
               'G':{'forward':'',
                    'reverse':'',
                    'index1':'',
                    'index2':''}
               }
    sourcep = {'C':{'read_threshold':'',
                    'window_size':'',
                    'mapq_threshold':'',
                    'start_threshold':'',
                    'gap_threshold':'',
                    'mismatch_threshold':'',
                    'merged_analysis':''},
               'G':{'demultiplex_min_reads':'',
                    'target':'',
                    'description':'',
                    'c_barcode1':'',
                    'c_barcode2':'',
                    'barcode1':'',
                    'barcode2':''}
               }
    input_inf = sourcei[method]
    para = sourcep[method]
    paras = para.keys()
    input_label = input_inf.keys()
    return input_inf,para,paras,input_label

def genyaml_C(input_inf,input_label,para,paras,output_add,ref,label):
    f = open('run/run.yaml','w')
    f.write('reference_genome: '+ref+'\n')
    f.write('analysis_folder: '+output_add+'\n\n')
    f.write('bwa: bwa\nsamtools: samtools\n\n')
    #print(para)
    #print(paras)
    for i in paras:
        f.write(i+': '+str(para[i])+'\n')
    f.write('\nsamples:\n')
    f.write(' '*4+label+':\n')
    for i in input_label:
        f.write(' '*8+i+': '+input_inf[i]+'\n')
    f.close()


def genyaml_G(input_inf,input_label,para,output_add,ref,label):
    f = open('run/run.yaml','w')
    f.write('reference_genome: '+ref+'\n')
    f.write('output_folder: '+output_add+'\n\n')
    f.write('bwa: bwa\nbedtools: bedtools\n\n')
    f.write('demultiplex_min_reads: '+str(para['demultiplex_min_reads'])+'\n\n')
    f.write('undemultiplexed:\n')
    for i in input_label:
        f.write(' '*4+i+': '+input_inf[i]+'\n')
    f.write('\nsamples:\n'+' '*4+'control:\n')
    f.write(' '*8+'target:\n')
    f.write(' '*8+'barcode1: '+para['c_barcode1']+'\n')
    f.write(' '*8+'barcode2: '+para['c_barcode2']+'\n')
    f.write(' '*8+'description: Control\n\n')
    f.write(' '*4+label+':\n')
    f.write(' '*8+'target: '+para['target']+'\n')
    f.write(' '*8+'barcode1: '+para['barcode1']+'\n')
    f.write(' '*8+'barcode2: '+para['barcode2']+'\n')
    f.write(' '*8+'description: '+para['description']+'\n')
    f.close()