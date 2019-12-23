import subprocess

def gen_conf(file='circos_input_file.txt',ref='hg19'):
    
    if ref == 'hg19':
        karyotype = 'karyotype.human.hg19.txt'
    else:
        print('no ref for circos existed')
    
    conf = [
            'karyotype = {0}'.format(karyotype),
            'chromosomes_units = 1000000',
            '<plots>',
            '<plot>',
            'type = text',
            'file = {0}'.format(file),
            'r1 = 1r',
            'r0 = 0.90r',
            'color = chr4',
            'label_font = bold',
            'label_size = 18p',
            'padding = 4p',
            '<rules>',
            '</rules>',
            '</plot>',
            '</plots>',
            '<<include etc/ideogram.conf>>',
            '<<include etc/ticks.conf>>',
            '<image>',
            '<<include etc/image.conf>>',
            '</image>',
            '<<include etc/colors_fonts_patterns.conf>>',
            '<<include etc/housekeeping.conf>>'
            ] 
    with open('circos.conf','w') as f:
        for i in conf:
            f.write(i)
            f.write('\n')
    
    
def gen_circos_data(stand_file):
    
    circos_data = []
    
    with open(stand_file,'r') as r_f:
        for l in r_f:
            ls = l.split('\t')
            loc = ls[0]
            sgRNA = ls[3]
            num = ls[1]
            chromosome,s_e = loc.split(':')
            start,end = s_e.split('-')
            description = sgRNA+'_'+num
            chromosome = chromosome.replace('chr','hs')
            circos_data.append('\t'.join([chromosome,start,end,description]))
    with open('circos_input_file.txt','w') as w_f:
        for i in circos_data:
            w_f.write(i)
            w_f.write('\n')
            
            
def circos(stand_file,output_path):
    gen_circos_data(stand_file)
    gen_conf()

    cmd = 'circos -conf circos.conf -debug_group summary,timer > run.out'
    subprocess.call(cmd, executable='/bin/bash', shell=True)
    if output_path.endswith('/'):
        output_path = output_path
    else:
        output_path += '/'
    cmd = 'mv circos.png circos.svg -t {0}'.format(output_path)
    subprocess.call(cmd, executable='/bin/bash', shell=True)
            
    