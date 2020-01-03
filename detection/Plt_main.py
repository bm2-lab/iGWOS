# -*- coding: utf-8 -*-
"""
Created on Sun Sep  3 14:35:08 2017

@author: xuedy
"""

import os
import sys
from Plt_parser import gen_arg, check_parser
from Plt_process import process_pack
from Plt_trans import trans_pack
from P3_file_upload import get_time
#import subprocess

if __name__ == '__main__':
    parser = gen_arg()
    args = parser.parse_args()

    ky = check_parser(args)
    if ky == 1:

        output_path = args.output
        if output_path.endswith('/'):
            output_path = output_path.rstrip('/')
            args.output = output_path

        # trans relative path to absolute path
        if not output_path.startswith('/'):
            output_path = os.getcwd()+'/'+output_path
            args.output = output_path
        """
        if args.reference == 'hg38':
            args.reference = 'data/hg38.fa'
        elif args.reference == 'hg19':
            args.reference = 'data/hg19.fa'
        """

        if os.path.exists(output_path):
            if len(os.listdir(output_path)) != 0:
                print('{0} is not empty, if insist set {1} as output path (yes/y no/n)'.format(output_path, output_path))
                j = sys.stdin.readline()
                j = j.strip()
                if j in ('yes', 'y'):
                    methods = args.type
                    ky = process_pack(methods)(args)
                    if ky == 1:
                        print('[{0}][INFO][{1}] finished'.format(get_time(),methods))
                    else:
                        print('[{0}][INFO][{1}] failed'.format(get_time(),methods))
            else:
                methods = args.type
                ky = process_pack(methods)(args)
                if ky == 1:
                    print('[{0}][INFO][{1}] finished'.format(get_time(),methods))
                else:
                    print('[{0}][INFO][{1}] failed'.format(get_time(),methods))
        else:
            os.makedirs(output_path)
            methods = args.type
            ky = process_pack(methods)(args)
            if ky == 1:
                print('[{0}][INFO][{1}] finished'.format(get_time(),methods))
            else:
                print('[{0}][INFO][{1}] failed'.format(get_time(),methods))
    else:
        print('file {0} does not exist'.format(' , '.join(ky)))
