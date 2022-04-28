## new script for second part of canopus (first is filtering mgf)
## Takes in mgf and various directories and alignment file
## very similar to old version

import sys, os

def launch_sirius(mgf, odir, mstype):
    '''
    :param mgf: is full path to folder with input ms files
    :param odir: is where sirius output will be
    '''
    command = f'sirius -o {odir} -i {mgf} formula -c 50 -p {mstype} '\
        '--ppm-max-ms2 5.0 --elements-enforced HCNOPS '\
            ' zodiac structure canopus '\

    os.system(command)

def read_alignment(alif):
    '''
    :param alif: path to alignment file
    :returns: dict of metabolites with id as key and ms1 as value
    '''
    odict = {}
    with open(alif) as alifil:
        lin = alifil.readline().strip()
        ct = 0
        while lin:
            linl = lin.split('\t')
            if linl[0][0].isdigit() and linl[31] != '' and '2+' not in linl[4] and '2-' not in linl[4]:
                odict[linl[0]] = [linl[30], linl[31], linl[4], linl[2]]
            ct += 1
            lin = alifil.readline().strip()

    return odict

def make_mgf_for_compound(mgfinpdir, mgfstr, alid, ct):
    '''
    :param mgfinpdir: is a string, path to where mgfs will be output
    :param mgfstr: is a string of lines to output
    :param alid: is a dict, output of read_alignment
    :param ct: int, will be part of directory name 
    output is: a new ms file for the compound in mgf directory
    returns: name of mgf outputted
    '''
    ##making new directory
    mgfl = mgfstr.split('\n')
    scans = mgfl[1].strip().split('=')[1]

    ##outputting new ms files
    mgfnam = f'{mgfinpdir}/{ct}_scan{scans}.ms'
    with open(mgfnam, 'w') as mgfout: 
        mgfout.write(f'>compound scans{scans}\n')
        for i, lin in enumerate(mgfl):
            if 'PEPMASS' in lin:
                mass = lin.split('=')[1]
                mgfout.write(f'>parentmass {mass}\n')
            elif 'ION=' in lin:
                ion = lin.split('=')[1]
                mgfout.write(f'>ionization {ion}\n')
                mgfout.write('>ms1\n')
                newms1 = alid[scans][0].replace(' ', '\n').replace(':', ' ')
                mgfout.write(f'{newms1}\n>ms2\n')
                for index in range(i+1, (len(mgfl)-1)):
                    newline = mgfl[index] + '\n'
                    if 'END' not in newline:
                        mgfout.write(newline)

    return mgfnam

def main(mgfinp, sirdir, mgfpardir, alifil, mstype):

    alid = read_alignment(alifil)

    ct = 1
    with open(mgfinp, 'r') as inp:
        inlin = inp.readline()
        while inlin:
            if inlin.startswith('BEGIN'):
                ostring = inlin
            elif inlin.startswith('END'):
                ostring += inlin
                inlin = inp.readline()
                ostring += inlin
                scans = ostring.split('\n')[1].strip().split('=')[1]
                if scans in alid.keys():
                    nam = make_mgf_for_compound(mgfpardir, ostring, alid, ct)
                    ct += 1
                
            else:
                ostring += inlin
            inlin = inp.readline()
        
    launch_sirius(mgfpardir, sirdir, mstype)
        
if __name__ == '__main__':

    if len(sys.argv) != 6:
        sys.exit('ARGS: 1) full path to input mgf\n '\
            '2) sirius dir, in which to output sirius result directories\n'\
                'NEEDS TO BE EMPTY DIRECTORY\n'\
                '3) directory for holding individual mgfs (cannot be sirius dir)\n'\
                    '4) Path to alignment file corresponding to input mgf\n'\
                        '5)  What MS type? qtof OR orbitrap')

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])   
    print('Done!') 