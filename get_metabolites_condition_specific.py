import sys, re
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
from statistics import median

def get_groups(coln, omics):
    '''
    :param coln: a list of strings
    :param omics: string "metab" or "trans"
    This function is hard coded to look for groups at
    [32:] of this list
    RETURNS: dict of columns in each condition
    '''
    outd = {}
    condl = []
    ind = 32 if omics == 'metab' else 1
    newl = [i.split('_')[0] for i in coln[ind:]]
    for i in newl:
        if i not in condl:
            condl.append(i)

    for con in condl:
        matching = []
        recon = re.compile(f'^{con}')
        for i in range(ind, len(coln)):
            if re.match(recon, coln[i]):
                matching.append(i)
        
        outd[con] = matching 

    return outd

def find_dists(infil, omics):
    '''
    :param infil: a string, path to input alignment file
    :param omics: string
    returns: float of 95th percentile of enrichment values
    '''
    fcl =[]
    with open(infil, 'r') as inf:
        inl = inf.readline()
        conds = get_groups(inl.strip().split('\t'), omics)
        while inl:
            if (omics == 'metab' and inl[0].isdigit()) or (omics == 'trans' and inl[0] == 'B'):
                tmpd = {}
                inlist = inl.strip().split('\t')
                for k, v in conds.items():
                    pkareas = [float(inlist[col]) for col in v]
                    med = median(pkareas)
                    tmpd[k] = med

                for cond, med in tmpd.items():
                    ## Recording fold change of max value / second max value
                    ## (How much higher is the value of the condition the metab
                    ## is highest in from the value of the condition the metab is 
                    ## second highest in)
                    allotherareas = [tmpd[k] for k in tmpd.keys() if k != cond]
                    fc = med / max(allotherareas)
                    if fc > 1:
                        fcl.append(fc)
            inl = inf.readline()
    
    return np.percentile(fcl, 95)

def get_specificity(inlist, conds, cutoff):
    '''
    :param inlist: list of peak areas for a metabolite
    :param conds: dict of conditions (from get_groups)
    :param cutoff: float for fold change cutoff
    Func checks for specificy of a metabolite
    If it is highly expressed in one condition only
    RETURNS: False OR a (string, float) tuple
    '''
    hits = 0
    goodcond = ''
    medarea = {}
    for k, v in conds.items():
        pkareas = [float(inlist[col]) for col in v]
        med = median(pkareas)
        medarea[k] = med

    for cond, area in medarea.items():
        allotherareas = [medarea[k] for k in medarea.keys() if k != cond]
        ## 14 is ROUGHLY equal to peak aboundance of 1e5. This checks to see that
        ## if a metab is in other conditions, it's lowly abundand.
        if area / max(allotherareas) >= cutoff and max(allotherareas) <= 14:
            return cond, area, area/max(allotherareas)

    return False

def get_cond_specific_metabs(infil, omics):
    '''
    :param infil: a string, path to input alignment file
    :param omics: a string
    This file is hard-coded to look for condition columns at [32:-1]
    RETURNS: dictionary of condition specific metabs
    conditions as key and RT, mz, adduct, average peak area, ms1 and ms2 as values
    '''
    outd = {}
    cutoff = find_dists(infil, omics)
    print('The Fold Change cutoff is: ', str(cutoff))

    with open(infil, 'r') as inf:
        inl = inf.readline()
        conds = get_groups(inl.strip().split('\t'), omics)
        print(conds)
        while inl:
            if (omics == 'metab' and inl[0].isdigit()) or (omics == 'trans' and inl[0] == 'B'):
                result = get_specificity(inl.strip().split('\t'), conds, cutoff)
                if result:
                    avgarea = result[1]
                    cond = result[0]
                    foldchange = result[2]
                    inlst = inl.strip().split('\t')                
                    if omics == 'trans':
                        outl = [inlst[0], str(avgarea), str(foldchange)]
                    if omics == 'metab':
                        alignid = inlst[0]; rt = inlst[1]
                        mz = inlst[2]; ms1 = inlst[30]
                        adduct = inlst[4]; ms2 = inlst[31]
                        outl = [alignid, rt, mz, adduct, ms1, ms2, str(avgarea), str(foldchange)]
                    if (omics == 'metab' and ms2 != '') or omics == 'trans':
                        if cond not in outd.keys():
                            outd[cond] = [outl]
                        else:
                            outd[cond].append(outl)

            inl = inf.readline()
    
    return outd

def write_out(uniqued, odir, omics):
    '''
    :param uniqued: a dict with conditions as keys and metab info as values
    :param odir: string, path to output file
    :param omics: trans or metab
    OUTPUTS: metabolite info for each metabolite
    '''
    if uniqued and omics == 'metab':
        with open(odir + '_enriched_metabolites.tab', 'w') as ofil:
            ofil.write('Alignment.ID\tRT\tmz\tAdduct\tMS1\tMS2\tAverageArea\tFoldChange\tCondition\n')
            for k, v in uniqued.items():
                for metab in v:
                    ostr = '\t'.join(metab)
                    ofil.write(f'{ostr}\t{k}\n')
    elif uniqued and omics == 'trans':
        with open(odir + '_enriched_transcripts.tab', 'w') as ofil:
            ofil.write('Gene\tAverageArea\tFoldChange\tCondition\n')
            for k, v in uniqued.items():
                for trans in v:
                    ostr = '\t'.join(trans)
                    ofil.write(f'{ostr}\t{k}\n')

    else:
        print('No condition-specific metabolites :(')

def make_hist(uniqued, odir):
    '''
    :param uniqued: dict with conditions as keys and metab info as values
    :param odir: is a string, path to output dir
    OUTPUTS: histogram with unique metabolite counts
    per condition
    '''
    ll = []; kl = []
    for k, v in uniqued.items():
        ll.append(len(v))
        kl.append(k.replace('.', ' '))

    print(kl); print(ll)
    y_pos = np.arange(len(kl))
    plt.figure()
    plt.bar(y_pos, ll, align='center', alpha=0.5)
    plt.xticks(y_pos, kl, rotation = 45)
    plt.ylabel('Number of Condition-Specific Metabolites')
    plt.savefig(odir + '_cond_specific_metabs_barchart.pdf', bbox_inches = 'tight')


def main(infil, odir, omics):
    '''
    all inputs are strings
    '''
    metabs = get_cond_specific_metabs(infil, omics)
    write_out(metabs, odir, omics)
    make_hist(metabs, odir)
    print('Done!')

if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit('ARGS: 1) input file (alignment file) to find condition\n'\
            'specific metabolites in 2) Desired output filename (and directory)\n'\
                '3) What omics type? trans OR metab')

    main(sys.argv[1], sys.argv[2], sys.argv[3])
