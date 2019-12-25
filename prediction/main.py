from parse import args
import os
import pandas as pd
import numpy as np
import pickle
from pyfaidx import Fasta
import tensorflow as tf
from DeepCRISPR.deepcrispr import DCModelOfftar
from CFD.otscore import calcCfdScore
from MIT.otscore import calcMitScore
from CROPIT.otscore import calcCropitScore
from CCTop.otscore import calcCcTopScore
from sklearn.ensemble import AdaBoostClassifier

# first, get POT list of gRNA based on Cas-OFFinder.
def cas_input(genome,gRNAs,mismatch):
    f=open('data/cas_input.txt','w')
    f.write(genome+'\n')
    f.write('NNNNNNNNNNNNNNNNNNNNNGG\n')
    for i in gRNAs:
        f.write(i[:20] + 'NNN '+str(mismatch)+'\n')
    f.close()
    print("Make gRNA sequence file for Cas-OFFinder prediction in data/cas_input.txt")
    os.system(" ./cas-offinder data/cas_input.txt G data/cas_output.txt")
    print("Obtain candidate off-target sites by Cas-OFFinder in data/cas_output.txt")
    #os.system("CUDA_VISIBLE_DEVICES="+str(gpu)+" ./cas-offinder data/cas_input.txt G data/cas_output.txt")


def pot(gid,gRNAs):
    f_cas=pd.read_csv('data/cas_output.txt',sep='\t',names=['pattern','Chr','Start','OTS','Strand','Mismatch'])
    f_cas['pattern'] = f_cas.pattern.apply(lambda x: x[:20])
    f_cas.Start= f_cas.Start.apply(lambda x: x+1)
    gRNA_dic={'sgID':gid, 'gRNA':gRNAs}
    f_gRNA = pd.DataFrame(gRNA_dic,columns=['sgID','gRNA'])
    f_gRNA['pattern'] = f_gRNA.gRNA.apply(lambda x: x[:20])
    f_pot = pd.merge(f_cas, f_gRNA, how='left', on='pattern')
    f_pot = f_pot.drop(['pattern'], axis=1)
    f_pot = f_pot.reindex(columns=['sgID', 'gRNA', 'OTS', 'Chr', 'Strand', 'Start', 'Mismatch'])
    f_pot.OTS = f_pot.OTS.str.upper()
    f_pot.to_csv('data/pot.tab', sep='\t', index=False)
    print("Make formatted candidate off-target file in data/pot.tab")
    f_gRNA = f_pot[f_pot.gRNA==f_pot.OTS]
    f_gRNA = f_gRNA.drop('OTS', axis=1)
    f_gRNA.to_csv('data/grna.tab', sep='\t', index=False)
    print("Make formatted gRNA file in data/grna.tab")
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
    print("Encode gRNAs and candidate off-target sites")

    f2 = pd.read_csv('data/pot.tab', usecols=[0, 2, 3, 4, 5], sep='\t', low_memory=False)
    input = f2.apply(lambda row: epi(row, ctcf, dnase, h3k4me3, rrbs, span), axis=1).tolist()
    x_ot_off_target = np.array(input).reshape(f2.shape[0], 8, 1, span)

    f3 = pd.read_csv('data/pot.tab', sep='\t', usecols=[0, 1], low_memory=False)
    f3 = pd.merge(f3, f1, how='inner', on=['sgID', 'gRNA'])
    input = f3.apply(lambda row: epi(row, ctcf, dnase, h3k4me3, rrbs, span), axis=1).tolist()
    x_sg_off_target = np.array(input).reshape(f3.shape[0], 8, 1, span)
    # return x_sg_off_target, x_ot_off_target
    fpkl = open('data/encode_off.pkl', 'wb')
    pickle.dump([x_sg_off_target, x_ot_off_target], fpkl)
    fpkl.close()


def deepots(f1, step = 1000):
    #os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu)
    f = open('data/encode_off.pkl', 'rb')
    x_sg_off_target, x_ot_off_target = pickle.load(f)
    f.close()
    print("Predict with DeepCRISPR")
    sess = tf.InteractiveSession()
    # using regression model, otherwise classification model
    off_target_model_dir = 'DeepCRISPR/trained_models/offtar_pt_cnn'
    is_reg = False
    dcmodel = DCModelOfftar(sess, off_target_model_dir, is_reg)
    offnum = len(x_ot_off_target)
    predicted_off_target = []
    for i in range(0, offnum, step):
        x_sg = x_sg_off_target[i:i + step]
        x_off = x_ot_off_target[i:i + step]
        predicted_off = dcmodel.offtar_predict(x_sg, x_off)
        predicted_off_target.extend(list(predicted_off))

    print("Number of predicted OTS", len(predicted_off_target))
    f1['DeepCRISPR'] = pd.Series(predicted_off_target)
    f1.to_csv('data/deepcrispr.tab', sep='\t', index=False)
    print("Obtain prediction result with DeepCRISPR in data/deepcrispr.tab")
    return f1

def igwosv(gRNA_path, f,output):
    #f = f[f.Mismatch > 0]
    # CFD, MIT, Cropit, and CCTop score
    f['CFD'] = f.apply(lambda row: calcCfdScore(row['gRNA'], row['OTS']), axis=1)
    f['MIT'] = f.apply(lambda row: calcMitScore(row['gRNA'], row['OTS']), axis=1)
    f['CROP-IT'] = f.apply(lambda row: calcCropitScore(row['gRNA'], row['OTS']), axis=1)
    f['CCTop'] = f.apply(lambda row: calcCcTopScore(row['gRNA'], row['OTS']), axis=1)
    # CRISPRoff score
    print('Calculate CRISPRoff score')
    os.system("./crisproff.sh "+gRNA_path)
    fcroff = pd.read_csv('data/crisproff.tab', sep='\t', low_memory=False)
    f = pd.merge(f, fcroff, how="left", on=['gRNA','OTS','Chr','Strand','Start'])
    # #uCRISPR score
    # os.system("./ucrispr.sh ")
    # fucr = pd.read_csv('data/ucrispr.out', sep=' ', low_memory=False)
    # f = pd.concat([f, fucr['uCRISPR']], axis=1)

    # get training data of vitro
    fv = open('data/data_vitro.pkl', 'rb')
    XV,YV = pickle.load(fv)
    fv.close()
    #form test data
    CRISPRoff = f['CRISPRoff']
    CFD = f['CFD']
    MIT = f['MIT']
    Cropit = f['CROP-IT']
    CCTop = f['CCTop']
    X = np.stack((CRISPRoff, CFD, MIT, Cropit, CCTop), axis=1)
    #Adaboost classifier to predict ots score
    print("Integrate CRISPRoff, CFD, MIT, Cropit, and CCTop prediciton scores")
    clf = AdaBoostClassifier(random_state=1, n_estimators=50, algorithm='SAMME.R')
    f['iGWOS'] = clf.fit(XV, YV).predict_proba(X)[:, 1]
    #f['iGWOS'] = f.iGWOS.apply(lambda x: x ** 3)
    order = ['sgID', 'gRNA', 'OTS', 'Chr', 'Strand', 'Start', 'Mismatch', 'iGWOS']
    f = f[order]
    f.to_csv('{0}/igwosv.tab'.format(output), sep="\t", index=False)
    print("Output iGWOS predicition result in {0}/igwosv.tab".format(output))
    print(f.describe())


def igwosc(gRNA_path,f,output):
    # f = f[f.Mismatch > 0]
    f['CFD'] = f.apply(lambda row: calcCfdScore(row['gRNA'], row['OTS']), axis=1)
    f['CROP-IT'] = f.apply(lambda row: calcCropitScore(row['gRNA'], row['OTS']), axis=1)
    # CRISPRoff score
    print('Calculate CRISPRoff score')
    os.system("./crisproff.sh " + gRNA_path)
    fcroff = pd.read_csv('data/crisproff.tab', sep='\t', low_memory=False)
    f = pd.merge(f, fcroff, how="left", on=['gRNA', 'OTS', 'Chr', 'Strand', 'Start'])
    # get training data of cell
    fc = open('data/data_cell.pkl', 'rb')
    XC, YC = pickle.load(fc)
    fc.close()
    # form test data
    CRISPRoff = f['CRISPRoff']
    DeepCRISPR = f['DeepCRISPR']
    CFD = f['CFD']
    Cropit = f['CROP-IT']
    X = np.stack((CRISPRoff, DeepCRISPR, CFD, Cropit), axis=1)
    # Adaboost classifier to predict ots score
    print("Integrate CRISPRoff, DeepCRISPR, CFD, and Cropit prediciton scores")
    clf = AdaBoostClassifier(random_state=1, n_estimators=50, algorithm='SAMME.R')
    f['iGWOS'] = clf.fit(XC, YC).predict_proba(X)[:, 1]
    #f['iGWOS'] = f.iGWOS.apply(lambda x: x**3)
    order=['sgID', 'gRNA', 'OTS', 'Chr', 'Strand', 'Start', 'Mismatch','iGWOS']
    f=f[order]
    f.to_csv('{0}/igwosc.tab'.format(output), sep="\t", index=False)
    print("Output iGWOS predicition result in {0}/igwosc.tab".format(output))
    print(f.describe())


# parse arguments
gRNA_path=args.gRNA
gen_path=args.genome
genome = os.path.split(gen_path)[1]
mismatch=args.mismatch
#gpu=args.gpu
out_path=args.output
cp=args.circos

if gen_path[-1] == '/':
    gen_path = gen_path[:-1]
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
cas_input(gen_path, gRNAs, mismatch)
f_gRNA, f_pot = pot(gid, gRNAs)

# parse technique type
ttype=args.type

if ttype=='VITRO':
    # in-vitro CIRCLE-seq with CRISPRoff, CFD, MIT, Cropit, and CCTop
    f_igwos=igwosv(gRNA_path,f_pot,out_path)
    if cp == 1:
        print("visualize the genome-wide off-target profile with the circos plot")
        os.system("./circos.sh {0}/igwosv.tab {1} {2}".format(out_path, genome,out_path))
elif ttype=='CELL':
    cell=args.cell
    cid_path=args.cid
    en_path=args.encode
    if en_path[-1] == '/':
        en_path = en_path[:-1]
    # encode ots and predict with deepcrispr
    f_cid=pd.read_csv(cid_path,sep='\t',names=['cid','cell'])
    cid=f_cid.cid[f_cid.cell==cell].tolist()[0]
    encode(f_gRNA, en_path, cid)
    f_deep = deepots(f_pot)
    # cell-based technique with CRISPRoff, DeepCRISPR,  CFD, and Cropit
    f_igwos =igwosc(gRNA_path,f_deep,out_path)
    if cp==1:
        print("visualize the genome-wide off-target profile with the circos plot")
        os.system("./circos.sh {0}/igwosc.tab {1} {2}".format(out_path,genome,out_path))