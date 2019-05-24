#!/bin/bash

for tech in /home/yjf/off-target/off-target/evaluate_tech/tech_OTS/bench_POT/*
do

tname=${tech##*/}
tname=${tname%.tab}
echo $tname
awk 'BEGIN{FS="\t";OFS=" "}NR>1{print $2,$3}' $tech > mydat/$tname.dat
./uCRISPR -off ./mydat/$tname.dat > mydat/$tname.out

awk 'NR==FNR{a[i]=$0;i++}NR>FNR{print a[j]"\t"$3;j++}' $tech mydat/$tname.out > POT/$tname.tab
sed -i 's/uCRISPR_score/uCRISPR/g' POT/$tname.tab
wc -l $tech
wc -l POT/$tname.tab
head POT/$tname.tab
done