cat /dev/null > data/circos_input.txt
awk 'BEGIN{FS="\t";OFS=" "}NR>1{print "hs"substr($4,4),$6,$6+22,$8}' $1 >> data/circos_input.txt
cd data
sed -i "s/hg[0-9]\+/$2/g" circos.conf
circos -conf circos.conf 