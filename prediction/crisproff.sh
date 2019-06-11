#!/bin/bash

## default data/grna.fa
gfile=$1
gname=${gfile##*/}

## Creating the index structure for genome sequences.
# risearch2.x -c /data/genome/fasta/hg19.fa -o CRISPRoff/hg19.suf

## off-target prediction
cp $gfile CRISPRoff/grna/
cd CRISPRoff/
if [ ! -d "outgz/" ];then
	mkdir outgz
fi

if [ ! -d "result/" ];then
	mkdir result
fi

rm -f outgz/*.out.gz
rm -f result/*.tsv
risearch2.x -q grna/$gname -i hg19.suf -s 1:20 -m 5:0 -e 10000 -l 0 --noGUseed -p3
mv *.out.gz outgz/


## CRISPRoff score
python2 CRISPRspec_CRISPRoff_pipeline.py --guides grna/$gname --risearch_results_folder outgz/ --no_azimuth  --CRISPRoff_scores_folder result/


## formatted OTS
echo -e 'gRNA\tOTS\tChr\tStrand\tStart\tCRISPRoff'> pot.tab
for file in result/*.tsv
do
	gseq=${file##*/}
	grna=${gseq:0:23}
	#echo $grna
	awk 'BEGIN{FS="\t";OFS="\t"}NR>3{if($1~/chr[0-9XY]*$/ && $4~/GG$/) print "'$grna'",toupper($4),$1,$6,$2+1,$5}' $file >> pot.tab
done

