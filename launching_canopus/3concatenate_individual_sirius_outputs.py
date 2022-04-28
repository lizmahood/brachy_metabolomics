import sys, glob, os

def concat_sirius(sirdir):
    '''
    @param: sirdir is a string, path to parent sirius directory
    Returns: 2 lists of lines in sirius output files (first is
    formula file, second is canopus file)
    '''
    compounds_dirs = os.listdir(sirdir)
    canopus = []; formulas = []
    for dir in compounds_dirs:
        os.chdir(f'{sirdir}/{dir}')
        try:
            with open('./formula_identifications.tsv', 'r') as forms:
                tmpl = [lin.strip() for lin in forms]
                formulas += tmpl
        
            with open('./canopus_summary.tsv', 'r') as canop:
                tmpc = [lin.strip() for lin in canop]
                canopus += tmpc
        except:
            print(dir)
            print('This probably didn\'t complete!')

        os.chdir(f'{sirdir}')
        
    return formulas, canopus

def remove_duplicate_lines(fillst):
    '''
    @param: fillst is a list of strings, first element is "header"
    returns: list with all but the first instance of "header" removed
    '''
    header = fillst[0]
    filtered = list(filter(lambda lin: lin != header, fillst))
    return [header] + filtered

def main(sirdir, onam):
    '''
    @param: sirdir is a string, path to parent sirius directory
    @param: onam is a string, path to output dir and name
    '''
    formulas, canopus = concat_sirius(sirdir)
    formgood = remove_duplicate_lines(formulas)
    canopgood = remove_duplicate_lines(canopus)

    with open(onam + '_formulas.tsv', 'w') as oform:
        for lin in formgood:
            oform.write(lin + '\n') 

    with open(onam + '_canopus_summary.tsv', 'w') as ocanop:
        for lin in canopgood:
            ocanop.write(lin + '\n') 

if __name__ == '__main__':

    if len(sys.argv) != 3:
        sys.exit('ARGS: 1) FULL directory FROM ROOT to sirius output files, must end in /\n'\
            ' 2) output name (not directory! no path!) for concatenated sirius output files\n'\
                'output gets put into directory of sirius output files')

    main(sys.argv[1], sys.argv[2])
    print('Done!')