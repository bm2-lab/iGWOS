import argparse
# iGWOS description and version
parser = argparse.ArgumentParser(prog='iGWOS', description="Prediction of CRISPR-Cas9 off-target sites with iGWOS.")
parser.add_argument('-v','--version', action='version', version='%(prog)s 1.0')
# gRNAs file input for off-target prediction
parser.add_argument('-gRNA',help='gRNAs file in Fasta format', default='data/grna.fa')
# run cas-offinder to get the candidate off-target sites.
parser.add_argument('-g', dest='genome', help='the genome folder for candidate off-target searching, default=genome/hg19', default='genome/hg19')
parser.add_argument('-m', dest='mismatch', help='the maximum mismatch allowed in off-target prediction, default=5', type=int, default=5, choices=range(6))
parser.add_argument('-gpu', help='select a gpu device to perform cas-offinder and/or deepcrispr, default=0', type=int, default=0)
# output file
parser.add_argument('-o',dest='output',help='the output folder, default=data/',default='data')

# OTS prediction on in-vitro or cell-based technique
subparsers = parser.add_subparsers(title='subcommands', description='select the type of OTS detection technique', dest='type')
# sub-command in-vitro
parser_v = subparsers.add_parser('VITRO', help='in-vitro CIRCLE-seq')
#parser_v.add_argument('-x', type=int, help='x value')
# sub-command cell-based
parser_c = subparsers.add_parser('CELL', help='cell-based techniques')
# deepcrispr prediction with gRNA encoded in certain cell-types.
parser_c.add_argument('-cell', help='the cell-type of gRNAs',default='K562')
parser_c.add_argument('-cid', help='the cell-id file, formed like data/encode_hg19.tab', default='data/encode_hg19.tab')
parser_c.add_argument('-e', dest='encode', help='the epigenomic encode folder, default=/data/genome/encode/fa/',default='/data/genome/encode/fa')

# encode all the options
args = parser.parse_args()