# -*- coding: utf-8 -*-
"""
Created on Sun Sep  3 14:35:08 2017

@author: xuedy
"""

import os
from Plt_parser import gen_arg
from Plt_process import process_pack
from Plt_trans import trans_pack
from P3_file_upload import ssh_scp_put, get_time
#import subprocess

parser = gen_arg()
args = parser.parse_args()


path = args.output
"""
if not path.startswith('/'):
    path = os.getcwd()+'/'+path
"""
if not os.path.exists(path):
    os.makedirs(path)
    
methods = args.method
process_pack(methods)(args)

if args.upload == True:
    uid = args.userid
    trans_pack(methods,uid)
    with open('uploadfile/upload_file_name.txt','r') as f_n:
        upf = f_n.readlines()[0]
    ssh_scp_put('xuedy','xuedy',upf)

print('[{0}][INFO][{1}] finished'.format(get_time(),methods))