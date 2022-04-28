import sys
import numpy as np
import pandas as pd

def get_average_pk_area_per_class(canopus, pkarea):
    '''
    :param canopus: pd df
    :param pkarea: pd df
    :returns: pd df
    '''

    canopus = canopus.loc[:, ['name', 'class']]
    pkarea = pkarea.iloc[:, np.r_[0, 32:pkarea.shape[1]]]

    both = canopus.merge(pkarea, how = 'left', left_on = 'name', right_on = 'Alignment.ID')

    ## excluding classes with <5 metabolites belonging to them
    classes = both[['class']].value_counts()
    classes = classes.rename_axis('unique_values').reset_index(name='counts')
    good = classes.loc[(classes['counts'] >= 5), 'unique_values'].to_list()
    both = both[(both['class'].isin(good))]

    ## aggregating to get mean of remaining classes per sample (column)
    out = both.groupby(['class']).mean()
    out = out.drop(['name', 'Alignment.ID'], axis = 1)

    ## cleaning up column names
    outcols = out.columns
    newcols = [x.split('_Run')[0] for x in outcols]
    out.columns = newcols
    return(out)

def main(canopus_s, pkarea_s, outn):
    '''
    All inputs strings
    '''
    canopus = pd.read_csv(canopus_s, sep = '\t')
    pkarea = pd.read_csv(pkarea_s, sep = '\t')

    out = get_average_pk_area_per_class(canopus, pkarea)

    out.to_csv(f'{outn}_mean_pk_area_per_class_for_BAR.tsv', sep = '\t')

if __name__ == '__main__':

    if len(sys.argv) != 4:
        sys.exit('ARGS: 1) Canopus summary, filtered 2) Filtered and'\
            'normalized peak area file 3) output prefix and name')

    main(sys.argv[1], sys.argv[2], sys.argv[3])

    print('Done!')