import subprocess
import os


def filterBackgroundSites(bedtools_path, sample_path, control_path, outfile):
    # modifed original code
    """
    output_folder = os.path.dirname(outfile)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    sample_r_path = sample_path.rstrip('.txt') + '_ex.txt'
    control_r_path = control_path.rstrip('.txt') + '_ex.txt'
    #print sample_r_path
    m0_command = 'cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33 {0}'.format(sample_path)
    m1_command = 'cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33 {0}'.format(control_path)
    with open(sample_r_path, 'w') as m_f:
        subprocess.call(m0_command.split(), stdout=m_f)
    with open(control_r_path, 'w') as m_f:
        subprocess.call(m1_command.split(), stdout=m_f)


    bedtools_filter_command = '{0} intersect -a {1} -b {2}'.format(bedtools_path, sample_r_path, control_r_path)

    with open(outfile, 'w') as outfile:
        subprocess.call(bedtools_filter_command.split(), stdout=outfile)
    """

    # written code
    def fbs_python(sample_path, control_path, outfile, strict=False):
        with open(sample_path,'r') as f:
            s = []
            for l in f:
                s.append(l.strip())

        with open(control_path,'r') as f:
            c = []
            for l in f:
                c.append(l.strip())
        filtered = []
        header = s[0]
        if strict:
            for l in s[1:]:
                if l not in c:
                    filtered.append(l)
        else:
            c_idx = ['*'.join(i.split('\t')[:3]) for i in c[1:]]
            s_idx = ['*'.join(i.split('\t')[:3]) for i in s[1:]]
            for i in range(len(s_idx)):
                if s_idx[i] not in c_idx:
                    filtered.append(s[i+1])

        filtered = [header] + filtered
        outfile_folder = '/'.join(outfile.split('/')[:-1]) + '/'
        if os.path.exists(outfile_folder) == False:
            os.mkdir(outfile_folder)
        with open(outfile,'w') as f:
            for l in filtered:
                f.write(l + '\n')

    fbs_python(sample_path, control_path, outfile)
