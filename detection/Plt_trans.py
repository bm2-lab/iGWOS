# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 21:28:09 2017

@author: xuedy
"""

import time
from S2 import findsgRNA

def write_upf(upload_data,date,uid,species,description,target_seq,methods):
    ulf_name = 'upload_file_{0}_{1}.txt'.format(uid,str(round(time.time()))[-4:])
    f = open('uploadfile/upload_file_name.txt','w')
    f.write(ulf_name)
    f.close()
    u = open('uploadfile/'+ulf_name,'w')
    for i in upload_data:
        u.write(','.join([date,uid,methods,species,description,target_seq,i[0],str(i[1])])+'\n')
    u.close()

def trans_guideseq(uid,yaml_file='run/run.yaml'):
    date = time.strftime('%Y-%m-%d',time.localtime(time.time()))
    y = open(yaml_file,'r')
    input_inf = y.readlines()
    y.close()
    species = input_inf[0].rstrip('.fa\n').split('/')[-1]
    output_path = input_inf[1].rstrip('\n').lstrip('output_folder: ')
    target_seq = input_inf[-4].lstrip('target: ').rstrip('\n')
    description = input_inf[-1].lstrip('description: ').rstrip('\n')
    name = input_inf[-5].lstrip(' ').rstrip(':\n')
    upload_data = []
    f = open(output_path+'/identified/{0}_identifiedOfftargets.txt'.format(name),'r')
    base_num = 1
    for l in f:
        if l.startswith('chr'):
            ll = l.split('\t')
            off_seq = ll[21]
            num = ll[11]
            if off_seq != '':
                if off_seq[:-3] == target_seq[:-3]:
                    base_num = float(num)
                else:
                    upload_data.append([off_seq,int(num)])
    f.close()
    for i in range(len(upload_data)):
        upload_data[i][1] = round(upload_data[i][1]/base_num,4)
    write_upf(upload_data,date,uid,species,description,target_seq,methods='GUIDE-seq')

def trans_circleseq(uid,yaml_file='run/run.yaml'):
    date = time.strftime('%Y-%m-%d',time.localtime(time.time()))
    y = open(yaml_file,'r')
    input_inf = y.readlines()
    y.close()
    species = input_inf[0].rstrip('.fa\n').split('/')[-1]
    output_path = input_inf[1].rstrip('\n').lstrip('analysis_folder: ')
    target_seq = input_inf[-6].lstrip('target: ').rstrip('\n')
    description = input_inf[-1].lstrip('description: ').rstrip('\n')
    name = input_inf[-7].lstrip(' ').rstrip(':\n')
    upload_data = []
    f = open(output_path+'/identified/{0}_identified_matched.txt'.format(name),'r')
    base_num = 1
    for l in f:
        ll = l.split('\t')
        off_seq = ll[11]
        num = ll[4]
        if off_seq != '':
            if off_seq[:-3] == target_seq[:-3]:
                base_num = float(num)
            else:
                upload_data.append([off_seq,int(num)])
    f.close()
    for i in range(len(upload_data)):
        upload_data[i][1] = round(upload_data[i][1]/base_num,4)
    write_upf(upload_data,date,uid,species,description,target_seq,methods='CIRCLE-seq')

def trans_siteseq(uid,args):
    species = args.reference.split('/')[-1].rstrip('.fa')
    description = args.description
    date = time.strftime('%Y-%m-%d', time.localtime(time.time()))
    f = open(args.output+'/outcome.txt','r')
    upload_data = []
    base_num = 1
    for l in f:
        if l.startswith('>'):
            num = l.rstrip('\n').split('|')[-1]
            num = float(num)
        else:
            off_seq = l.rstrip('\n')
            off_seq = findsgRNA(off_seq,args.target)
            upload_data.append([off_seq,int(num)])
    if target_seq == None:
        max_n = 0
        n = 0
        for i in range(len(upload_data)):
            if len(upload_data[i][1]) > n:
                n = len(upload_data[i][1])
                max_n = i
        target_seq = upload_data[max_n][0]
        base_num = upload_data[max_n][1]
    else:
        target_seq = args.target
        for i in range(len(upload_data)):
            if upload_data[i][0] == target_seq:
                base_num = upload_data[i][1]
                max_n = i
    upload_data.pop(max_n)
    for i in range(len(upload_data)):
        upload_data[i][1] = round(upload_data[i][1] / base_num, 4)

    write_upf(upload_data,date,uid,species,description,target_seq,methods='SITE-seq')

    
def circle_standard(input_file,output_file,target):
    data = []
    n = 0
    idx = 1
    with open(input_file,'r') as f:
        for l in f:
            if n == 0:
                n += 1
            else:
                ll = l.split('\t')
                r_idx = str(idx)
                loc = 'chr'+ll[3]
                num = ll[4]
                chain = ll[12]
                mismatch_n = ll[11]
                sgRNA = ll[10]
                data.append([loc,num,chain,sgRNA,mismatch_n,target])
                idx += 1
    with open(output_file,'w') as f_o:
        f_o.write('Location\tReads\tStrand\tCleavage_seq\tMismatch\tTagret_seq\n')
        for i in data:
            f_o.write('\t'.join(i)+'\n')

def guide_standard(input_file,output_file,target):
    data = []
    idx = 1
    with open(input_file,'r') as f:
        n = 0
        for l in f:
            if n == 0:
                n = 1
            else:
                ll = l.split('\t')
                sgRNA = ll[21]
                if sgRNA != '':
                    r_idx = str(idx)
                    loc = ll[0]
                    num = ll[11]
                    chain = ll[29]
                    mismatch_n = ll[22]
                    sgRNA = ll[21]
                    data.append([loc,num,chain,sgRNA,mismatch_n,target])
                    idx += 1
    with open(output_file,'w') as f_o:
        f_o.write('Location\tReads\tStrand\tCleavage_seq\tMismatch\tTagret_seq\n')
        for i in data:
            f_o.write('\t'.join(i)+'\n')


def trans_pack(method,uid):
    trpk = {'GUIDE-seq':trans_guideseq(uid),
            'CIRCLE-seq':trans_circleseq(uid),
            'SITE-seq':trans_siteseq(uid)}
    return trpk[method]
