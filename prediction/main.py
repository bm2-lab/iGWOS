from parse import args
import os
import pandas as pd
import numpy as np
from pyfaidx import Fasta
import tensorflow as tf
from DeepCRISPR.deepcrispr import DCModelOfftar
from CFD.otscore import calcCfdScore
from MIT.otscore import calcMitScore
from CROPIT.otscore import calcCropitScore


# first, get POT list of gRNA based on Cas-OFFinder.
def cas_input(genome,gRNAs,mismatch,gpu):
    f=open('data/cas_input.txt','w')
    f.write(genome+'\n')
    f.write('NNNNNNNNNNNNNNNNNNNNNGG\n')
    for i in gRNAs:
        f.write(i[:20] + 'NNN '+str(mismatch)+'\n')
    f.close()

    os.system("CUDA_VISIBLE_DEVICES="+str(gpu)+" ./cas-offinder data/cas_input.txt G data/cas_output.txt")


def pot(gid,gRNAs):
    f_cas=pd.read_csv('data/cas_output.txt',sep='\t',names=['pattern','Chr','Start','OTS','Strand','Mismatch'])
    f_cas['pattern'] = f_cas.pattern.apply(lambda x: x[:20])
    f_cas.Start= f_cas.Start.apply(lambda x: x+1)
    gRNA_dic={'sgID':gid, 'gRNA':gRNAs}
    f_gRNA = pd.DataFrame(gRNA_dic,columns=['sgID','gRNA'])
    f_gRNA['pattern'] = f_gRNA.gRNA.apply(lambda x: x[:20])
    f_pot = pd.merge(f_cas, f_gRNA, how='left', on='pattern')
    f_pot.drop(['pattern'], axis=1, inplace=True)
    f_pot = f_pot.reindex(columns=['sgID', 'gRNA', 'OTS', 'Chr', 'Strand', 'Start', 'Mismatch'])
    f_pot.OTS = f_pot.OTS.str.upper()
    f_pot.to_csv('data/pot.tab', sep='\t', index=False)
    f_gRNA = f_pot[f_pot.gRNA==f_pot.OTS]
    f_gRNA.drop('OTS', axis=1, inplace=True)
    f_gRNA.to_csv('data/gRNA.tab', sep='\t', index=False)
    return f_gRNA, f_pot


# encode the ots
def epi(gline, ctcf, dnase, h3k4me3, rrbs, span):
    rna = gline[1]
    Chr = gline[2]
    Strand = gline[3]
    Start = gline[4]
    ls = []
    for epi in ['A', 'C', 'G', 'T']:
        ls.append([1 if rna[i] == epi else 0 for i in range(span)])
    for epi in ['ctcf', 'dnase', 'h3k4me3', 'rrbs']:
        f = eval(epi)
        if Strand == '-':
            seq = f[Chr][Start:Start + span].seq[::-1]
        else:
            seq = f[Chr][Start:Start + span].seq
        ls.append([1 if seq[i] == 'A' else 0 for i in range(span)])
    return ls


def encode(f1, en_path, cid, func=epi):
    span = 23
    ctcf = Fasta('{0}/{1}_ctcf.fa'.format(en_path, cid))
    dnase = Fasta('{0}/{1}_dnase.fa'.format(en_path, cid))
    h3k4me3 = Fasta('{0}/{1}_h3k4me3.fa'.format(en_path, cid))
    rrbs = Fasta('{0}/{1}_rrbs.fa'.format(en_path, cid))

    f2 = pd.read_csv('data/pot.tab', usecols=[0, 2, 3, 4, 5], sep='\t', low_memory=False)
    input = f2.apply(lambda row: epi(row, ctcf, dnase, h3k4me3, rrbs, span), axis=1).tolist()
    x_ot_off_target = np.array(input).reshape(f2.shape[0], 8, 1, span)

    f3 = pd.read_csv('data/pot.tab', sep='\t', usecols=[0, 1], low_memory=False)
    f3 = pd.merge(f3, f1, how='inner', on=['sgID', 'gRNA'])
    input = f3.apply(lambda row: epi(row, ctcf, dnase, h3k4me3, rrbs, span), axis=1).tolist()
    x_sg_off_target = np.array(input).reshape(f3.shape[0], 8, 1, span)

    # fpkl = open('data/encode_off.pkl', 'wb')
    # pickle.dump([x_sg_off_target, x_ot_off_target], fpkl)
    # fpkl.close()
    return x_sg_off_target, x_ot_off_target


def deepots(gpu, f1, x_sg_off_target, x_ot_off_target, step = 1000):
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu)
    # f = open('data/encode_off.pkl', 'rb')
    # x_sg_off_target, x_ot_off_target = pickle.load(f)
    # f.close()
    sess = tf.InteractiveSession()
    # using regression model, otherwise classification model
    off_target_model_dir = 'DeepCRISPR/trained_models/offtar_pt_cnnofftar_pt_cnn'
    is_reg = False
    dcmodel = DCModelOfftar(sess, off_target_model_dir, is_reg)
    offnum = len(x_ot_off_target)
    predicted_off_target = []
    for i in range(0, offnum, step):
        x_sg = x_sg_off_target[i:i + step]
        x_off = x_ot_off_target[i:i + step]
        predicted_off = dcmodel.offtar_predict(x_sg, x_off)
        predicted_off_target.extend(list(predicted_off))

    print("Num of predicted OTS", len(predicted_off_target))
    f1['DeepCRISPR'] = pd.Series(predicted_off_target)
    f1.to_csv('data/deepcrispr.tab', sep='\t', index=False)
    return f1


def score_join(gRNA_path,f):
    #f = f[f.Mismatch > 0]

    #CFD, MIT, CROP-IT score
    f['CFD'] = f.apply(lambda row: calcCfdScore(row['gRNA'], row['OTS']), axis=1)
    f['MIT'] = f.apply(lambda row: calcMitScore(row['gRNA'], row['OTS']), axis=1)
    f['CROP-IT'] = f.apply(lambda row: calcCropitScore(row['gRNA'], row['OTS']), axis=1)

    #crisproff score
    os.system("./crisproff.sh "+gRNA_path)
    fcroff = pd.read_csv('data/crisproff.tab', sep='\t', low_memory=False)
    f = pd.merge(f, fcroff, how="left", on=['gRNA','OTS','Chr','Strand','Start'])

    #uCRISPR score
    os.system("./ucrispr.sh ")
    fucr = pd.read_csv('data/ucrispr.out', sep=' ', low_memory=False)
    f = pd.concat([f, fucr['uCRISPR']], axis=1)

    # save join result
    f.to_csv('data/join.tab', sep="\t", index=False)
    print(f.describe())
    return(f)


def integ(f,output):

    # normalize
    #f['CRISPRseek'] = f['CRISPRseek'] / 100
    #f['CCTop'] = f['CCTop'] / 224
    f['CROP-IT'] = f['CROP-IT'] / 650
    f['MIT'] = f['MIT'] / 100
    f['uCRISPR'] = f['uCRISPR'] / 6.5
    f['CRISPRoff'] = (f['CRISPRoff'] + 40) / 75

    # integrate prediction to iGWOS
    tlst = ["CRISPRoff", "uCRISPR", "MIT", "CFD", "DeepCRISPR", "CROP-IT"]
    weightper = [102.732, 99.673, 98.702, 92.501, 87.856, 82.08]
    weight = [weightper[i] / 100 for i in range(0, len(weightper))]
    f['iGWOS_weight'] = np.sum(f[tlst[i]] * weight[i] for i in range(6))

    weight_rank=range(1,7)[::-1]
    f['iGWOS_rank'] = np.sum(f[tlst[i]] * weight_rank[i] for i in range(6))
    f.to_csv('{0}/igwos.tab'.format(output), sep="\t", index=False)


gRNA_path=args.gRNA
cell=args.cell
gen_path=args.genome
mismatch=args.mismatch
gpu=args.gpu
cid_path=args.cid
en_path=args.encode
out_path=args.output

if gen_path[-1] == '/':
    gen_path = gen_path[:-1]
if en_path[-1] == '/':
    en_path = en_path[:-1]
if out_path[-1] == '/':
    out_path = out_path[:-1]

if not os.path.exists(out_path):
    os.mkdir(out_path)

# read fa file of gRNAs.
fa_gRNA = Fasta(gRNA_path, sequence_always_upper=True)
gid = list(fa_gRNA.keys())
gRNAs = [fa_gRNA[i][:].seq for i in gid]
fa_gRNA.close()

# first, get POT list of gRNA based on Cas-OFFinder.
cas_input(gen_path, gRNAs, mismatch, gpu)
f_gRNA, f_pot = pot(gid, gRNAs)

# second, encode ots and predict with deepcrispr
f_cid=pd.read_csv(cid_path,sep='\t',names=['cid','cell'])
cid=f_cid.cid[f_cid.cell==cell].tolist()[0]

sg_ots,ot_ots = encode(f_gRNA, en_path, cid)
f_deep = deepots(gpu,f_pot,sg_ots,ot_ots)

# third, join tools' prediction scores
f_join = score_join(gRNA_path, f_deep)

#fourth, integrate to iGWOS score
integ(f_join, out_path)