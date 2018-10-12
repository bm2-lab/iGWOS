# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 16:02:49 2017

@author: xuedy
"""

import subprocess
import os
try:
    from pyfaidx import Fasta
except:
    print('*'*20+'\ncan not import pyfaidx\n'+'*'*20)
from P3_file_upload import get_time
from S1 import find_initial_read_pileups, call_site_seq_features

def at():
    thread = os.popen("grep 'processor' /proc/cpuinfo | sort -u | wc -l").read().rstrip('\n')
    thread = int(thread)
    return round(thread/2)

def site_data_analysis(read1, read2, output, target, ref, p=at()):
    """
    if file.endswith('.sra'):
        print ('[{0}][INFO][SITEseq] unpack the sequencing data'.format(get_time()))
        try:
            cmd = 'fastq-dump {0}'.format(file)
        except:
            print('there is something wrong with fastq-dump')
        subprocess.call(cmd, executable='/bin/bash', shell=True)
        file = file.rstrip('.sra')+'.fastq'
    """
    name = read1.rstrip('.fastq').split('_')[0]
    cmd = []
    cmd.append("bowtie2 -p {0} -x hg38 --gbar 30 -1 {1} -2 {2} -S {3}".format(p,read1,read2,name+'.sam'))
    cmd.append("samtools view -bS {0} > {1}".format(name+'.sam',name+'.bam'))
    cmd.append("samtools sort {0} -o {1}".format(name+'.bam',name+'s'))
    cmd.append("samtools index {0}".format(name+'s.bam'))
    for i in cmd:
        try:
            subprocess.call(i, executable='/bin/bash', shell=True)
        except:
            print('there is something wrong when executing the order "{0}" , maybe the path of {1} not in the $PATH'.format(i,i.split(' ')[0]))
        print ('[{0}][INFO][SITEseq] {1}'.format(get_time(),i))
    res = find_initial_read_pileups(name+'s.bam')
    try:
        ind = Fasta(ref)
    except:
        print ('something wrong with the reference')
    print ('[{0}][INFO][SITEseq] analysing'.format(get_time()))
    call_site_seq_features(res,name+'s.bam',ind,output+'/intermediate_outcome.txt')
    u_d = trans_SITEseq_output_to_standard_output(output+'/intermediate_outcome.txt',output+'/outcome.txt',target)
    offtargets = [[i[2],i[4]] for i in u_d]

    return offtargets


def trans_SITEseq_output_to_standard_output(input_file,output_file,target):

    with open(input_file,'r') as f:
        upload_data = []
        for l in f:
            if l.startswith('>'):
                ll = l.rstrip('\n').split('|')
                num = ll[-1]
                num = int(num)
                chromo = ll[1].split(':')[0]
                loc = int(ll[1].split(':')[1])
                chain = '*'
            else:
                off_seq = l.rstrip('\n')
                try:
                    sgRNA,mismatch_n,chain = findsgRNA(off_seq,target)
                    if chain == '+':
                        l_b_sg = len(off_seq.split(sgRNA)[0])
                        loc = loc-20+l_b_sg
                    else:
                        sgRNAs = rev(sgRNA)
                        l_b_sg = len(off_seq.split(sgRNAs)[1])
                        loc = loc+20-l_b_sg
                except:
                    print('there is a sequence with no sgRNA in it')
            upload_data.append([chromo,loc,off_seq,chain,num,mismatch_n])
    with open(output_file,'w') as fo:
        for i in upload_data:
            chromo,loc,off_seq,chain,num,mismatch_n = i
            if chain == '+':
                start_loc = loc
                end_loc = loc+23
            else:
                start_loc = loc-22
                end_loc = loc+1
            fo.write('{0}:{1}-{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(chromo,start_loc,end_loc,num,chain,sgRNA,mismatch_n,target))
            
    return upload_data

def count_different(seq1,seq2):
    n = 0
    for i in range(20):
        if seq1[i] != seq2[i]:
            n += 1
    return n


def fastq_merge(file,out_name):
    data = []
    for i in file:
        with open(i,'r') as f:
            for l in f:
                data.append(l)
    with open(out_name,'w') as o_f:
        for i in data:
            o_f.write(i)

def rev(seq):

    seq = seq.replace('A', 'X')
    seq = seq.replace('T', 'A')
    seq = seq.replace('X', 'T')
    seq = seq.replace('C', 'X')
    seq = seq.replace('G', 'C')
    seq = seq.replace('X', 'G')
    seq = seq[::-1]

    return seq


def findsgRNA(seq,target):
    
    if seq.islower():
        seq = seq.upper()
    """
    seq_u = seq[:20]
    seq_d = seq[20:]
    if 'GG' in seq_d[:8]:
        PAM = seq_d[:8].split('GG')[0]
        l = len(PAM)+2
        ls = l-23
        sgRNA = seq_u[ls:]+PAM+'GG'
    elif 'CC' in seq_u[-8:]:
        PAM = seq_u[-8:].split('CC')[-1]
        l = len(PAM)+2
        ls = 23-l
        sgRNA = 'CC'+PAM+seq_d[:ls]
        sgRNA = rev(sgRNA)
    """

    score_l = []
    seqs = rev(seq)
    for i in range(18):
        gg_loc = 21+i
        if seq[gg_loc:gg_loc+2] == 'GG':
            score = count_different(seq[gg_loc-21:gg_loc-1],target)
            score_l.append([seq[gg_loc-21:gg_loc+2],score,'+'])
        if seqs[gg_loc:gg_loc+2] == 'GG':
            score = count_different(seqs[gg_loc-21:gg_loc-1],target)
            score_l.append([seqs[gg_loc-21:gg_loc+2],score,'-'])
            
    min_diff = 20
    min_diff_i = 0
    for i in range(len(score_l)):
        num = score_l[i][1]
        if num < min_diff:
            min_diff = num
            min_diff_i = i
    sgRNA,mismatch_n,chain = score_l[min_diff_i]
    
    try:
        return sgRNA,mismatch_n,chain
    except:
        print('there is a sequence with no PAM in it')
        
def SITE_ot_2_visl(offtargets,target):
    visual = []
    for i in offtargets:
        chromo,loc,off_seq,chain,num,mismatch_n = i
        d = {}
        d['reads'] = num
        d['realigned_ref_seq'] = target
        d['seq'] = off_seq
        visual.append(d)
    return visual
