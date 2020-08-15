#!/bin/bash
awk 'BEGIN{FS="\t";OFS=" "}NR>1{print "hs"substr($4,4),$6,$7,$9}' $1 |sort -r -k 4 |awk 'NR<201{print $0}' > data/circos_input.txt
sed -i "s/hg[0-9]\+/$2/g" circos.conf
circos -conf circos.conf --outputdir $3 -debug_group summary,timer > data/circos.log