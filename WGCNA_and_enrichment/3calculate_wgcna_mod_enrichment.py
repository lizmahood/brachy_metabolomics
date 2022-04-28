import sys, os, pandas as pd
from types import MethodDescriptorType
from fisher import FisherExactTest
from calculate_canopus_enrichment import fisher
from sre_compile import isstring
from statsmodels.stats.multitest import multipletests

## Find out what inputs are
## For each WGCNA module, make tuple of: (+mod+class, +mod-class, -mod+class, -mod-class)
## Put it into fisher
## Output: df with ModName, ClassName, 'NumClassMod', 'NumNotClassMod', 'NumClassNotMod', 'NumNotClassNotMod','En', 'pValue'

def make_fisher_input(modfile, canopus, pkarea):
    '''
    :param modfile: string, path to WGCNA module file
    :param canopus: string, path to filtered canopus summary
    :param pkarea: string, path to peak area file
    '''

    moddf = pd.read_csv(modfile, sep = '\t')
    canopusdf = pd.read_csv(canopus, sep = '\t')
    pkareadf = pd.read_csv(pkarea, sep = '\t')

    pkareadf = pkareadf.dropna(axis = 0, subset = ['MS.MS.spectrum'])
    ## The MS.MS.spectrum col isn't important here. Just need to keep 
    ## pkareadf a df
    pkareadf = pkareadf[['Alignment.ID', 'MS.MS.spectrum']]

    canopusdf = canopusdf[['name', 'class']]

    allcanopus = pkareadf.merge(canopusdf, how = 'left', right_on = 'name', left_on = 'Alignment.ID')

    both = moddf.merge(allcanopus, how = 'right', left_on = 'X', right_on = 'name')
    both['class'] = both['class'].fillna('None')

    ## Dict of dicts to contain fisher inputs per class + module
    odict = {}
    uniqmod = both['selectColors'].unique().tolist()

    for mod in uniqmod:
        moddict = {}
        canopus_mod = both[(both['selectColors'] == mod)]
        canopus_notmod = both[(both['selectColors'] != mod)]

        ## getting list and valuecounts of the DAM's classes
        uniqclass = canopus_mod['class'].unique().tolist()
        mod_counts = canopus_mod['class'].value_counts()
        notmod_counts = canopus_notmod['class'].value_counts()

        for clas in uniqclass:
            cls_mod = mod_counts[clas]
            nocls_mod = canopus_mod.shape[0] - cls_mod
            try:
                cls_nomod = notmod_counts[clas]
            except:
                cls_nomod = 0
            nocls_nomod = canopus_notmod.shape[0] - cls_nomod
            moddict[clas] = [cls_mod, nocls_mod, cls_nomod, nocls_nomod]

        odict[mod] = moddict
        clean_dict = {k: odict[k] for k in odict if isstring(k)}


    return clean_dict

def main(modfile, canopus, pkarea, odir):
    '''
    all inputs strings
    '''
    fisher_dict = make_fisher_input(modfile, canopus, pkarea)
    fullout = pd.DataFrame(columns = ['NumClassDAM', 'NumNotClassDAM', 
    'NumClassNotDAM', 'NumNotClassNotDAM','En', 'pValue', 'Mod', 'Class'])

    for modn, mdict in fisher_dict.items():
        for cnam, ctups in mdict.items():
            pv, en = fisher(ctups)
            ctups.append(pv)
            ctups.append(en)

        out = pd.DataFrame.from_dict(mdict, orient = 'index')
        out.columns = ['NumClassDAM', 'NumNotClassDAM', 'NumClassNotDAM', 'NumNotClassNotDAM','En', 'pValue']
        out['Mod'] = modn
        out['Class'] = list(mdict.keys())
        out = out.sort_values('pValue')
        fullout = fullout.append(out)
    
    correctedp = multipletests(fullout['pValue'].to_list(), method = 'fdr_bh')
    fullout['correctedp'] = correctedp[1]
    fullout.to_csv(f'{odir}_WGCNA_module_class_enrichment.tsv', sep = '\t', index = False)

if __name__ == '__main__':

    if len(sys.argv) != 5:
        sys.exit('ARGS: 1) WGCNA module output file 2) canopus file\n'\
            ' 3) peak area file 4) output path and prefix')

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

    print('Done!')