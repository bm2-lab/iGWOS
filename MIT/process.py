# /usr/bin/python2.7
from otscore import calcMitScore
for tech in ["5gRNA","12gRNA"]:
    namelst = ['sgID', 'gRNA', 'OTS', 'Mismatch','Prediction']
    f = open('OT/{0}.tab'.format(tech) , 'w')
    f.write('{0}\n'.format('\t'.join(namelst)))
    f1 = open('../Data/{0}/POTseq.tab'.format(tech),'r')
    f1.readline()
    for line in f1.readlines():
        l = line.strip().split('\t')
        sgseq = l[1]
        ot = l[2]
        score=str(calcMitScore(sgseq,ot))
        l.append(score)
        f.write('\t'.join(l) + '\n')
    f1.close()
    f.close()
