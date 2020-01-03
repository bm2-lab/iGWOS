import subprocess

def gen_conf(file='circos_input.txt',ref='hg19'):
    
    if ref == 'hg19':
        karyotype = 'data/karyotype.human.hg19.txt'
    else:
        print('No reference genome existing for Circos plot')
    
    conf = [
            'karyotype = {0}'.format(karyotype),
            'chromosomes_units = 1000000',
            '<plots>',
            'file = {0}'.format(file),
            'type = scatter',
            '<plot>',
            'r1 = 0.99r',
            'r0 = 0.60r',
            'color = ddgrey',
            'stroke_color     = black',
            'stroke_thickness = 1',
            'glyph            = circle',
            'glyph_size       = 15p',
            'max   = 100000',
            'min   = 0',
			'<backgrounds>',
			'<background>',
			'color     = vlred',
			'y0        = 0.8r',
			'</background>',
			'<background>',
			'color     = vlyellow',
			'y1        = 0.8r',
			'y0        = 0.4r',
			'</background>',
			'<background>',
			'color     = vvlgreen',
			'y1        = 0.4r',
			'</background>',
			'</backgrounds>',
			'<axes>',
			'<axis>',
			'color     = lred',
			'thickness = 1',
			'spacing   = 0.05r',
			'y0        = 0.8r',
			'</axis>',
			'<axis>',
			'color     = dred',
			'thickness = 2',
			'spacing   = 0.1r',
			'y0        = 0.8r',
			'</axis>',
			'<axis>',
			'color     = lyellow',
			'thickness = 1',
			'spacing   = 0.05r',
			'y1        = 0.8r',
			'y0        = 0.4r',
			'</axis>',
			'<axis>',
			'color     = dyellow',
			'thickness = 2',
			'spacing   = 0.1r',
			'y1        = 0.8r',
			'y0        = 0.4r',
			'</axis>',
			'<axis>',
			'color     = lgreen',
			'thickness = 1',
			'spacing   = 0.05r',
			'y1        = 0.4r',
			'</axis>',
			'<axis>',
			'color     = dgreen',
			'thickness = 2',
			'spacing   = 0.1r',
			'y1        = 0.4r',
			'</axis>',
			'</axes>',
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
            if l.startswith('Location') == False:
                ls = l.split('\t')
                loc = ls[0]
                sgRNA = ls[3]
                num = ls[1]
                chromosome,s_e = loc.split(':')
                start,end = s_e.split('-')
                description = num
                chromosome = chromosome.replace('chr','hs')
                circos_data.append(' '.join([chromosome,start,end,description]))
    with open('circos_input.txt','w') as w_f:
        for i in circos_data:
            w_f.write(i)
            w_f.write('\n')
            
            
def circos(stand_file,output_path):
    gen_circos_data(stand_file)
    gen_conf()

    cmd = 'circos -conf circos.conf -debug_group summary,timer > run.out'
    subprocess.call(cmd, executable='/bin/bash', shell=True)
    if output_path.endswith('/') == False:
        output_path += '/'
    cmd = 'mv circos.png circos.svg -t {0}'.format(output_path+'visualization/')
    subprocess.call(cmd, executable='/bin/bash', shell=True)
            
    
