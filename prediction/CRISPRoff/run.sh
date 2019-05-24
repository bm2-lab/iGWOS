
## Creating the index structure for genome sequences.
#risearch2.x -c hg19.fa -o /home/yjf/off-target/off-target/tool/CRISPRoff/crisproff-1.1.1/hg19.suf


## off-target prediction
#for tech in `ls ../bench_fa`
#do
#
#mkdir outgz/$tname
#mkdir result/$tname
#tname=${tech%%.fa}
#echo $tname
#risearch2.x -q ../bench_fa/$tech -i hg19.suf -s 1:20 -m 5:0 -e 10000 -l 0 --noGUseed -p3
#mv *.out.gz outgz/$tname/
#
#done


## CRISPRoff score
#for tech in `ls outgz`
#do
#
#echo $tech
#python2 CRISPRspec_CRISPRoff_pipeline.py --guides ../bench_fa/$tech.fa --risearch_results_folder outgz/$tech/ --no_azimuth  --CRISPRoff_scores_folder result/$tech/
#
#done


# formatted OTS
for tech in `ls result`
do

echo -e 'gRNA\tOTS\tChr\tStrand\tStart\tCRISPRoff'> OTS/${tech}.tab
for file in result/$tech/*
do

gseq=${file##*/}
grna=${gseq:0:23}
echo $grna
awk 'BEGIN{FS="\t";OFS="\t"}NR>3{if($1~/chr[0-9XY]*$/ && $4~/GG$/) print "'$grna'",toupper($4),$1,$6,$2+1,$5}' $file >> OTS/${tech}.tab

done

awk 'END{print NR-1}' OTS/${tech}.tab
done