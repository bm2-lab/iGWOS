import argparse

parser = argparse.ArgumentParser(prog='iGWOS', description="Prediction of CRISPR-Cas9 off-target sites with iGWOS.")
parser.add_argument('-v','--version', action='version', version='%(prog)s 1.0')
# gRNAs file input for off-target prediction
parser.add_argument('-gRNA',help='gRNAs file in Fasta format', default='data/grna.fa')
# run deepcrispr with gRNA in certain cell-type
parser.add_argument('-cell', help='the cell-type of performed gRNAs',default='K562')
# run cas-offinder to get the candidate off-target sites.
parser.add_argument('-g', dest='genome', help='the genome folder for candidate off-target searching, default=genome/hg19', default='genome/hg19')
parser.add_argument('-m', dest='mismatch', help='the maximum mismatch allowed in off-target prediction, default=5', type=int, default=5, choices=range(6))
parser.add_argument('-gpu', help='select a gpu device to perform cas-offinder and deepcrispr, default=0', type=int, default=0)
# deepcrispr prediction with gRNA encoded.
parser.add_argument('-cid', help='the cell-id file, formed like data/encode_hg19.tab', default='data/encode_hg19.tab')
parser.add_argument('-e', dest='encode', help='the encode folder, default=/data/genome/encode',default='/data/genome/encode')
# output file
parser.add_argument('-o',dest='output',help='the output folder, default=data/',default='data')
# encode all the options
args = parser.parse_args()