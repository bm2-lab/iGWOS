import re
import os
from pyfaidx import Fasta
from otscore import calcCcTopScore
for file in os.listdir('../../OT/'):
    tech = file.split('.')[0]
    namelst = ['sgID', 'seq', 'mismatch', 'chr', 'strand', 'start', 'end', 'span', 'score']
    f = open('OT/' + file, 'w')
    f.write('{0}\n'.format('\t'.join(namelst)))
    f1 = open('../../OT/'+file,'r')
    f1.readline()
    f2 = Fasta('../../fa/'+tech,sequence_always_upper=True)
    sgs = {}
    for sg in f2.keys():
        sgs[sg] = []
    for line in f1.readlines():
        l = line.strip().split('\t')
        l1 = l[0].split('-')
        sg = '-'.join(l1[:-1])
        sgseq = f2[sg][:].seq.upper()
        ot = l[1].upper()
        score = str(calcCcTopScore(sgseq,ot))
        l.append(score)
        # print(l)
        f.write('\t'.join(l) + '\n')
        # print('\t'.join(l))
    f1.close()
    f.close()