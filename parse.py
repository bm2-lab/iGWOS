import argparse
# iGWOS description and version
parser = argparse.ArgumentParser(prog='iGWOS', description="Predict genome-wide CRISPR-Cas9 off-target sites with iGWOS.")
parser.add_argument('-v','--version', action='version', version='%(prog)s 1.0')
# gRNAs file input for off-target prediction
parser.add_argument('-gRNA',help='gRNAs file in Fasta format', default='data/grna.fa')
# run cas-offinder to get the candidate off-target sites.
parser.add_argument('-g', dest='genome', help='genome folder for candidate off-target searching, default=genome/hg19', default='genome/hg19')
parser.add_argument('-m', dest='mismatch', help='maximum mismatch allowed in off-target prediction, default=5', type=int, default=5, choices=range(7))
# deepcrispr prediction with gRNA encoded in certain cell-types.
parser.add_argument('-cell', help='cell-type of gRNAs',default='K562')
parser.add_argument('-cid', help='cell-id file, formed like data/encode_hg19.tab', default='data/encode_hg19.tab')
parser.add_argument('-e', dest='encode', help='epigenomic encode folder, default=/data/genome/encode/fa/',default='/data/genome/encode/fa')

# circos plot
parser.add_argument('-circos',dest='circos',help='whether to draw a circos plot to visualize the top 200 risky predicted off-target profile, default=1',type=int, default=1,choices=[0,1])
parser.add_argument('-gpu', help='select a gpu device to perform cas-offinder and/or deepcrispr, default=0', type=int, default=0)
# output file
parser.add_argument('-o',dest='output',help='output folder, default=output/',default='output')


# encode all the options
args = parser.parse_args()