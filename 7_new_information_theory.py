from __future__ import division
import pandas as pd
pd.options.mode.chained_assignment = None # to get rid of FP SettingWithCopyWarning
import sys, math, numpy as np
from statistics import mean
from scipy import stats

def make_df(inp, typ):
    '''
    :param inp: string, path to input
    :param typ: string, metab or trans 
    :returns: input, with nonessential cols removed
    '''
    df = pd.read_table(inp)
    if typ == 'metab':
        df.drop(df.columns[1:32],axis=1,inplace=True)
        df = df.replace(0.01, 0)
    return df

def make_pij(inp_df):
    '''
    :param inp_df: df, output of make_df
    :returns: df with normalized metabolite values
    '''
    num_pks = []
    ncol = inp_df.shape[1]
    for col in range(1, ncol):
        total = inp_df.iloc[:, col].sum()
        inp_df.iloc[:, col] = inp_df.iloc[:, col]/total
        num_pks.append(len(inp_df[inp_df.iloc[:, col] > 0]))

    inp_df['Pi'] = inp_df.iloc[:, 1:ncol].mean(axis = 1)

    return inp_df, num_pks

def find_shannon(inp_df):
    '''
    :param inp_df: output of make_pij
    :returns: list of shannon entropy per sample
    '''
    shannon = []
    ncol = inp_df.shape[1]-1
    for col in range(1, ncol):
        tmp = inp_df.iloc[:,col]
        tmp = tmp[tmp > 0]
        tmp = tmp.apply(lambda x: x * np.log2(x))
        shannon.append(-sum(tmp))

    return shannon

def find_si(inp_df):
    '''
    :param inp_df: output of make_pij
    :returns: inp_df with another column, Si
    '''
    ncol = inp_df.shape[1]-1
    values = inp_df.iloc[:, 1:ncol].div(inp_df.Pi, axis = 0)
    sidf = values * np.log2(values)
    inp_df['Si'] = sidf.mean(axis = 1)

    return inp_df

def find_di(inp_df):
    '''
    :param inp_df: output of find_si
    :returns: list of di per sample
    '''
    di = []
    ncol = inp_df.shape[1]-2
    for col in range(1, ncol):
        di_temp = inp_df.iloc[:, col] * inp_df['Si']
        di.append(sum(di_temp))

    return di

def find_rdpi(control, test, modus = 'mean'):
    '''
    :param control: df of the control columns
    :param test: df of the test columns
    :param modus: string, abssum OR mean OR dist
    :returns: df of RDPI per stress
    '''
    out = pd.DataFrame(columns = ['Names', 'RDPI'])
    for ctrl in range(control.shape[1]):
        ctrlc = control.columns[ctrl]
        for tst in range(test.shape[1]):
            tstc = test.columns[tst]
            newdf = pd.concat([control[ctrlc], test[tstc]], axis = 1)
            nonzero = newdf[(newdf.sum(axis = 1) > 0)]
            if 'Run' in ctrlc:
                newctrl = '_'.join(ctrlc.split('_')[:2])
                newtest = '_'.join(tstc.split('_')[:2])
                name = f'{newctrl}_{newtest}'
            else: name = f'{ctrlc}_{tstc}'
            if nonzero.shape[0] > 0:
                if modus == 'dist':
                    nonzero['diff'] = nonzero.iloc[:, 1] - nonzero.iloc[:, 0]
                else:
                    topp = abs(nonzero.iloc[:, 0] - nonzero.iloc[:, 1])
                    bottom = nonzero.iloc[:, 0] + nonzero.iloc[:, 1]
                    nonzero['diff'] = topp / bottom
                if modus == 'mean':
                    to_append = pd.Series([name, mean(nonzero['diff'])], index = out.columns)
                elif modus == 'abssum':
                    to_append = pd.Series([name, nonzero['diff'].sum(axis = 0)], index = out.columns)
                elif modus == 'dist':
                    to_append = pd.Series([name, [nonzero['diff'].sort_values()]], index = out.columns)
            else: 
                to_append = pd.Series([name, 0], index = out.columns)
            out = out.append(to_append, ignore_index=True)
    
    return out

def find_rdpi_superclass(control, test, aid, superclass, all_rdpi, rdpi_dist):
    '''
    :param control: df of control values
    :param test: df of test (stress) values
    :param aid: df, 1 col -- Alignment IDs
    :param superclass: df of the predicted class and superclass of metabolites
    :param all_rdpi: df
    :param rdpi_dist: df
    '''
    ccol, tcol = control.shape[1], test.shape[1]
    abunddf = pd.concat([aid, control, test], axis = 1)

    merged = pd.DataFrame.merge(abunddf, superclass, 
    how = 'left', left_on= 'Alignment.ID', right_on = 'name')

    merged['class'] = merged['class'].fillna('None')
    scs = merged['class'].unique()
    rdpi, peakarea_change, pkarea_percent = [], [], []

    metabdf = pd.DataFrame(columns = ['Alignment.ID', 'AvgPkAreaChange', 'Class'])
    clssmetdf = pd.DataFrame(columns = ['Class', 'Names', 'RDPI'])
    bigsups = []

    for sup in scs:
        changed, tmpl, perc_tmpl = {}, [], []
        supdf = merged[(merged['class'] == sup)]
        ## New check 6.20.21
        if supdf.shape[0] >= 5:
            bigsups.append(sup)
            ## Getting regular RDPI (think violin plots) per each class
            tmppdf = find_rdpi(supdf.iloc[:, 1:ccol+1], 
            supdf.iloc[:, ccol+1:tcol+ccol+1], modus = 'mean')
            tmppdf.insert(loc = 0, column = 'Class', value = sup)
            clssmetdf = clssmetdf.append(tmppdf, ignore_index = True)

            ## Getting indices of metabs present in ctrl and/or stress for
            ## logging IDs of these metabs
            tmpid_df = supdf.iloc[:, 1:tcol+ccol+1]
            tmpid_df['Mean'] = tmpid_df.mean(axis = 1)
            nonzero_inds = tmpid_df.index[tmpid_df['Mean'] > 0]
            alignids = supdf.loc[nonzero_inds, 'Alignment.ID']
            changed['Alignment.ID'] = alignids.tolist()
            print('Number of non-0 alignment IDs for this condition and superclass, via mean')
            print(len(alignids.to_list()))

            ## Now getting other metrics (similar to RDPI) per class
            for ctrl in range(1, ccol+1):
                ctrlc = supdf.columns[ctrl]
                for tst in range(ccol+1, tcol + ccol + 1):
                    tstc = supdf.columns[tst]
                    newdf = pd.concat([supdf[ctrlc], supdf[tstc]], axis = 1)
                    nonzero = newdf[(newdf.sum(axis = 1) > 0)]

                    ## debugging 
                    print(f'Number of non-0 AID for {ctrl} and {tst}')
                    print(nonzero.shape[0])
                    hmm = supdf.columns[0]
                    bug_newdf = pd.concat([supdf[hmm], supdf[ctrlc], supdf[tstc]], axis = 1)

                    bug_nonzero = bug_newdf[(bug_newdf.iloc[:, [1,2]].sum(axis = 1) > 0)] 
                    print(np.setdiff1d(alignids.tolist(), bug_nonzero['Alignment.ID'].tolist()))





                    ## For logging how much of the total RDPI is due to each superclass
                    topp = abs(nonzero.iloc[:, 0] - nonzero.iloc[:, 1])
                    bottom = nonzero.iloc[:, 0] + nonzero.iloc[:, 1]
                    nonzero['absdiff'] = topp / bottom
                    tmpl.append(nonzero['absdiff'].sum(axis = 0))
                    ## For getting the peak area changes of the metabolites in each class
                    ## And the percentile of the average peak area change
                    try:
                        avgchange = mean(nonzero.iloc[:, 1] - nonzero.iloc[:, 0])
                        changed[f'XC{ctrl}T{tst}'] = list(nonzero.iloc[:, 1] - nonzero.iloc[:, 0])
                    except: avgchange = 0
                    perc_tmpl.append(avgchange) 
            
            ## Averaging across replicates
            divs = (tmpl / all_rdpi['RDPI']) 
            if len(divs) < 5:
                print(divs)
            rdpi.append(mean(divs))
            peakarea_change.append(mean(perc_tmpl))

            all_changes = pd.DataFrame.from_dict(changed)
            all_changes_rep_means = all_changes.drop('Alignment.ID', axis = 1).mean(axis = 1)
            this_sup = pd.Series([sup] * all_changes_rep_means.shape[0], dtype='object')
            tmpout = pd.concat([all_changes['Alignment.ID'], all_changes_rep_means, this_sup], axis = 1)
            tmpout.columns = ['Alignment.ID', 'AvgPkAreaChange', 'Class']
            metabdf = metabdf.append(tmpout, ignore_index = True)

            ## Percentiles of average change per superclass
            tmpl2 = []
            for _ in range(len(perc_tmpl)):
                dist = rdpi_dist.loc[_, 'RDPI'][0].tolist()
                avg = perc_tmpl[_]
                avgperc = stats.percentileofscore(dist, avg, kind = 'mean')
                tmpl2.append(avgperc)
            pkarea_percent.append(mean(tmpl2))

    out = pd.concat([pd.Series(bigsups), pd.Series(rdpi)], axis = 1)
    out_pa = pd.concat([pd.Series(bigsups), pd.Series(pkarea_percent), 
    pd.Series(peakarea_change)], axis = 1)
    return (out, out_pa, metabdf, clssmetdf)

def put_test_with_control(inp_df, superclass = 'none'):
    '''
    :param inp_df: non-normalized input, the output of make_df
    :param ofil: string, path to output location and prefix of output name
    :returns: df of [test_sample_no, RDPI]
    '''
    if superclass != 'none':
        scdf = pd.read_table(superclass, sep = '\t')
        aid = inp_df['Alignment.ID']
        outdf = pd.DataFrame(columns = ['Stress', 'MetabClass', 'RDPI', 'Avg_PeakHeight_change', 'ChangePercentile'])
        metaboutdf = pd.DataFrame(columns = ['Alignment.ID', 'Stress', 'Class', 'AvgPkAreaChange'])
        clasoutdf = pd.DataFrame(columns = ['Class', 'Names', 'RDPI'])

    else: outdf = pd.DataFrame(columns = ['Names', 'RDPI'])
    df = inp_df.drop(inp_df.columns[0:1], axis = 1)

    ## No more hard-coded conditions
    groups = set([x.split('_')[0] for x in df.columns[1:]])
    tis = set([x.split('.')[-1] for x in groups])
    setts = set([x.split('.')[0] for x in groups])
    cold = {}
    for x in setts:
        for y in tis:
            ## Making sure that there is a control and stress condition in each set-tissue pair
            condnames = [string for string in df.columns[1:].tolist() if x in string and y in string]
            print(condnames)
            condds = set([t.split('.')[1] for t in condnames])
            if len(condds) > 1 and 'Ctrl' in condds:
                cold[f'{x}.{y}'] = []
                for col in range(len(df.columns.tolist())):
                    if x in df.columns[col] and y in df.columns[col]:
                        cold[f'{x}.{y}'].append(col)

    for k, v in cold.items():
        print(k)
        tmpdf = df.iloc[:, v]
        if len(tmpdf.columns) > 0:
            ctrl_conds = tmpdf.columns.str.contains('Ctrl')
            ctrl = tmpdf.loc[:,ctrl_conds].copy()
            tmpdf.drop(tmpdf.columns[ctrl_conds], axis = 1, inplace = True)
            conds = tmpdf.columns.str.split('_').str[0]
            uconds = np.unique(conds)
            for cond in uconds:
                tmpstress = tmpdf.loc[:, tmpdf.columns.str.contains(cond)]
                if superclass != 'none':
                    all_rdpi = find_rdpi(ctrl, tmpstress, modus = 'abssum')
                    rdpi_dist = find_rdpi(ctrl, tmpstress, modus = 'dist')
                    tmpout_sc, peakarea_sc, metabdf, tmpclassout = find_rdpi_superclass(ctrl, 
                    tmpstress, aid, scdf, all_rdpi, rdpi_dist)
                    names = pd.Series([cond] * tmpout_sc.shape[0])
                    metabnames = pd.Series([cond] * metabdf.shape[0])
                    mtempout = pd.concat([metabdf.iloc[:,0], metabnames, metabdf.iloc[:,1], metabdf.iloc[:,2]], axis = 1)
                    mtempout.columns = ['Alignment.ID', 'Stress', 'AvgPkAreaChange', 'Class']
                    tmpout = pd.concat([names, tmpout_sc, peakarea_sc.iloc[:,2], peakarea_sc.iloc[:,1]], axis = 1)
                    tmpout.columns = ['Stress', 'MetabClass', 'RDPI', 'Avg_PeakHeight_change', 'ChangePercentile']
                    metaboutdf = metaboutdf.append(mtempout, ignore_index = True)
                    clasoutdf = clasoutdf.append(tmpclassout, ignore_index=True)
                else: tmpout = find_rdpi(ctrl, tmpstress, 'mean')
                outdf = outdf.append(tmpout, ignore_index = True)
                
    if superclass != 'none':
        return (outdf, metaboutdf, clasoutdf)
    else: return outdf
            
def main(inp, ofil, typ, superclass):
    '''
    :param inp: string
    :param ofil: string
    :param typ: string
    '''
    df = make_df(inp,typ)
    print(df.head)
    if superclass != 'none':   
        rdpidf, metabdf, clasrdpidf = put_test_with_control(df, superclass)
        pd.DataFrame.to_csv(metabdf, ofil + '_pk_area_change_per_metabolite.tsv', sep = '\t', index = False)
    else: rdpidf = put_test_with_control(df, superclass)
    df, num_pks = make_pij(df)
    shannon = find_shannon(df)
    df = find_si(df)
    di = find_di(df)

    pd.DataFrame.to_csv(df, ofil + 'info_theory_per_' + typ + '.tsv', sep='\t')
    pd.DataFrame.to_csv(rdpidf, ofil + '_RDPI_per_stress_replicate.tsv', sep = '\t')
    if superclass != 'none':
        pd.DataFrame.to_csv(clasrdpidf, 
        ofil + '_RDPI_per_class_and_stress_replicate.tsv', sep = '\t')

    cols_to_write = list(df.columns.values.tolist())[1:-2]
    with open(ofil + '_shannon_di_per_sample.tsv', 'w') as out:
        out.write('SampleName\tShannon\tTotalPeaks\tdi\n')
        for val in range(len(cols_to_write)):
            out.write(f'{cols_to_write[val]}\t{shannon[val]}\t{num_pks[val]}\t{di[val]}\n')

if __name__ == '__main__':
    if len(sys.argv) != 5:
        sys.exit('ARGS: 1) Input normalized peak area alignment file, \n'\
            'with negative values turned to 0 2) output name 3) metab OR trans\n'\
            '4) Do you want RDPI broken down by superclasses? none OR path to canopus_summary file')

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

    print('Done!')
