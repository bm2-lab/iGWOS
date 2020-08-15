#!/bin/bash

## default data/grna.fa
gfile=$1
mis=$2
gname=${gfile##*/}

## Creating the index structure for genome sequences.
# risearch2.x -c /data/genome/fasta/hg19.fa -o CRISPRoff/hg19.suf

## off-target prediction
if [ ! -d "CRISPRoff/data/" ];then
	mkdir CRISPRoff/data
fi

cp $gfile CRISPRoff/data/
cd CRISPRoff/

if [ ! -d "outgz/" ];then
	mkdir outgz
fi

if [ ! -d "result/" ];then
	mkdir result
fi

rm -f outgz/*.out.gz
rm -f result/*.tsv

risearch2.x -q data/$gname -i hg19.suf -s 1:20 -m $mis:0 -e 10000 -l 0 --noGUseed -p3
mv *.out.gz outgz/


## CRISPRoff score
python2 CRISPRspec_CRISPRoff_pipeline.py --guides data/$gname --risearch_results_folder outgz/ --no_azimuth  --CRISPRoff_scores_folder result/


## formatted OTS
echo -e 'gRNA\tOTS\tChr\tStrand\tStart\tCRISPRoff'> data/crisproff.tab
for file in result/*.tsv
do
	gseq=${file##*/}
	grna=${gseq:0:23}
	#echo $grna
	awk 'BEGIN{FS="\t";OFS="\t"}NR>3{if($1~/chr[0-9XY]*$/ && $4~/GG$/) print "'$grna'",toupper($4),$1,$6,$2+1,$5}' $file >> data/crisproff.tab
done

cp data/crisproff.tab ../data/

