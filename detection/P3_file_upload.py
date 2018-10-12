# -*- coding: utf-8 -*-
"""
Created on Sun Aug 27 18:55:18 2017

@author: xuedy
"""

import time
#from ftplib import FTP
import paramiko
import os
import datetime

def get_time():
    now = datetime.datetime.now()
    rt = now.strftime('%Y-%m-%d %H:%M:%S').split(' ')
    d = rt[0].split('-')
    t = rt[1].split(':')
    if int(t[0]) > 21:
        t[0] = str(int(t[0])-12)
        ts = ':'.join(t)+'PM'
    elif int(t[0]) < 13:
        ts = ':'.join(t)+'AM'
    else:
        t[0] = str(int(t[0])-12)
        ts = '0'+':'.join(t)+'PM' 
    ds = '/'.join(d[1:])
    return '{0} {1}'.format(ds,ts)

date = time.strftime('%Y-%m-%d',time.localtime(time.time()))
"""
import subprocess
def up_load(file,uid,pw):
    cmd = "scp {0} {1}@192.168.1.101:/home/xuedy/database".format(file,uid)
    subprocess.call(cmd, executable='/bin/bash', shell=True)
    cmd = pw
    subprocess.call(cmd, executable='/bin/bash', shell=True)
    
"""    

def ssh_scp_put(user,password,file):
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect('10.11.41.109', 10122, user, password)
    a = ssh.exec_command('date')
    stdin, stdout, stderr = a
    print (stdout.read())
    sftp = paramiko.SFTPClient.from_transport(ssh.get_transport())
    sftp = ssh.open_sftp()
    sftp.put(os.getcwd()+'/'+file, '/home/xuedy/database/'+file)