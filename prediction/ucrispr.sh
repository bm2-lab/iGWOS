#!/bin/bash

pot='data/pot.tab'
tname='ucrispr'

cp $pot /usr/local/tool/uCRISPR/data/
cd /usr/local/tool/uCRISPR

awk 'BEGIN{FS="\t";OFS=" "}NR>1{print $2,$3}' data/pot.tab > data/$tname.dat

./uCRISPR -off ./data/$tname.dat > data/$tname.out

#awk 'NR==FNR{a[i]=$0;i++}NR>FNR{print a[j]"\t"$3;j++}' data/pot.tab data/$tname.out > data/$tname.tab

sed -i 's/uCRISPR_score/uCRISPR/g' data/$tname.out

cp data/$tname.out /root/project/data/
