import sys, os, glob
import pandas as pd
from joypy import joyplot 
import matplotlib.pyplot as plt

## Output with unknown classes: 4-col file of probability
## of the associated class for each instance. Ridge plot
## of probabilities per class

## Output with known class: 2-col file of TP/FN and probabilities
## for each instance, ridge plot of TP/FN probabilities

## Output with additional false positive class (can't be the same
## as the actual class): Additional file (one col) of FP probabilites
## and a density plot of these

def parse_canopus_outputs(c_summary, f_summary):
    '''
    :param c_summary: canopus_summary for all compounds
    :param f_summary: formula summary file for all compounds
    '''
    canopus = pd.read_table(c_summary)
    forms = pd.read_table(f_summary)
    both = pd.merge(canopus, forms, left_on = 'name', right_on = 'id')
    if 'FEATURE' in both.iloc[0]['name']:
        both['newname'] = both['name'].str.slice(2, -10)
    else:
        both['ind'] = both['name'].str.split('_').str[1]
        both['scans'] = both['name'].str.split('_').str[2]
        both['newname'] = both['ind'] + '_' + both['scans']
    out = both[['newname', 'precursorFormula', 'class']].copy()
    odict = out.set_index('newname').T.to_dict('list')

    return odict

def make_joyplot(toplot, odir, title):
    '''
    :param toplot: df with columns 'class' and 'probs'
    :param odir: string, complete name of output plot
    :param title: string, to pass into the title
    '''
    plt.figure()
    length = len(pd.unique(toplot['class']))
    joyplot(toplot, by='class',figsize=(8, length), ylim = 'own', overlap = 0, fade = True)
    plt.title(title)
    plt.savefig(odir)

def main(sirdir, actual_c, fp_c, mode, c_summary, f_summary, odir):
    '''
    :param sirdir: string, path to sirius output folders
    :param actual_c: string, actual class or 'Unknown'
    :param fp_c: string, false positive class or 'none'
    :param mode: string, 'pos' or 'neg'
    :param c_summary: string, path to canopus summary file
    :param f_summary: string, path to formula summary file
    :param odir: string, desired directory and prefix for all output file names
    '''
    sumd = parse_canopus_outputs(c_summary, f_summary)
    alldirs = list(os.walk(f'{sirdir}'))[0][1]
    if actual_c != 'Unknown':
        tp = {}; fn = {}
        if fp_c != 'none': 
            fp = {}

    ## reaading in canopus (or canopus_neg) file
    firstdir = list(sumd.keys())[0]
    index = [idx for idx, s in enumerate(alldirs) if firstdir in s][0]

    if mode == 'pos':
        all_classes = pd.read_table(f'{sirdir}/{alldirs[index]}/canopus.tsv')
    else:
        all_classes = pd.read_table(f'{sirdir}/{alldirs[index]}/canopus_neg.tsv')

    for dr in alldirs:
        ## going into each results directory and getting the probability of the actual class
        if dr in sumd.keys():
            form = sumd[dr][0]
            idd = dr.split('scan')[1]
            subdir = list(os.walk(f'{sirdir}/{dr}'))[0][1][0]
            ftps = list(os.walk(f'{sirdir}/{dr}/{subdir}/canopus/'))[0][2]
            try:
                ftp = [x for x in ftps if form in x][0]
                probs = pd.read_table(f'{sirdir}/{dr}/{subdir}/canopus/{ftp}', header = None)
                all_classes = all_classes.assign(probabilities = probs.values)
                if actual_c != 'Unknown':
                    thisprob = all_classes.loc[all_classes['name'] == actual_c]['probabilities'].values[0]
                    if actual_c == sumd[dr][1]:
                        tp[idd] = thisprob
                    else:
                        fn[idd] = thisprob
                        if fp_c == sumd[dr][1]:
                            fprob = all_classes.loc[all_classes['name'] == fp_c]['probabilities'].values[0]
                            fp[idd] = fprob
                else: 
                    try:
                        thisprob = all_classes.loc[all_classes['name'] == sumd[dr][1]]['probabilities'].values[0]
                        sumd[dr].append(thisprob)
                    except:
                        continue
            except:
                print('Canopus formula doesn\'t agree with sirius formula')
                print(ftp)
                print(form)
    
    if actual_c == 'Unknown':
        df = pd.DataFrame.from_dict(sumd, orient = 'index', columns = ['formula', 'class', 'prob'])
        newnames = [x.split('scan')[1] for x in list(df.index.values)]
        df['newnames'] = newnames
        df = df.set_index('newnames')
        df.to_csv(f'{odir}_probs_all_classes.tsv', sep = '\t')
        toplot = df[['class', 'prob']]
        title = 'Posterior Probability per Class'
        make_joyplot(toplot, f'{odir}_all_classes.pdf', title)
    
    else:
        tpdf = pd.DataFrame.from_dict(tp, orient = 'index', columns = ['Prob'])
        tpdf['class'] = 'TP'
        fndf = pd.DataFrame.from_dict(fn, orient = 'index', columns = ['Prob'])
        fndf['class'] = 'FN'
        dists = tpdf.append(fndf)
        dists.to_csv(f'{odir}_TP_FN_{actual_c}.tsv', sep = '\t')
        title = f'TP and FN Probabilities of Class: {actual_c}'
        make_joyplot(dists, f'{odir}_TP_FN_{actual_c}.pdf', title)
        if fp_c != 'none':
            fpout = pd.DataFrame(fp, columns = ['FalseProb'])
            fpout.to_csv(f'{odir}_FP_{fp_c}.tsv')

    print('Done!')

if __name__ == '__main__':

    if len(sys.argv) != 8:
        sys.exit('ARGS: 1) sirius output directory 2) Is there an actual class of these compounds?\n'\
            'Put name of class (with _ for spaces) OR Unknown 3) Is there a class you\'re looking \n'\
                'for False Positives of? Name of class OR none 4) pos OR neg 5) path to canopus \n'\
                    'summary file 6) Path to formula summary file 7) Output directory and prefix of all files')

    actual_c = sys.argv[2].replace('_', ' ')
    fp_c = sys.argv[3].replace('_', ' ')

    main(sys.argv[1], actual_c, fp_c, sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])

    print('Done!')    
