import sys
import numpy as np
import pandas as pd
import seaborn as sns

def read_in_mgf(mgfin, typ):
    '''
    :param mgfin: string, path to input
    :param typ: string, jgi OR ours
    '''
    odict = {}
    with open(mgfin, 'r') as mgff:
        lin = mgff.readline()
        while lin:
            if (lin.startswith('SCAN') and typ == 'ours') or (lin.startswith('FEAT') and typ == 'jgi'):
                name = lin.strip().split('=')[1]
                odict[name] = []
                peaks = []
                abunds = []
            elif lin.startswith('PEP'):
                pm = float(lin.strip().split('=')[1])
                odict[name].append(pm)
            elif lin.startswith('RT'):
                rt = float(lin.strip().split('=')[1])
                if 'SECOND' in lin:
                    rt = round(rt / 60, 5)
                odict[name].append(rt)
            elif lin[0].isdigit():
                abund = float(lin.strip().split()[1])
                if abund >= 3000:
                    peaks.append(float(lin.split()[0]))
                    abunds.append(abund)
            elif lin.startswith('END'):
                odict[name].append(peaks)
                odict[name].append(abunds)
            lin = mgff.readline()
    
    return odict

def centroid_peaks(pks):
    '''
    :param pks: 2D array of fragment peak mz and intensity
    :returns: 2D array but with peaks < 0.01 Da centroided
    '''
    inds_to_del = []
    for pk in range(pks.shape[1] -1):
        thispk = pks[0][pk]
        nextpk = pks[0][pk + 1]
        tmpinds = []
        while round(nextpk - thispk, 5) <= 0.01:
            tmpinds.extend([pk, pk + 1])
            pk += 1
            if pk < pks.shape[1] - 2:
                thispk = pks[0][pk]
                nextpk = pks[0][pk + 1]
            else:
                break
        
        if len(tmpinds) > 1:
            ## find position of peak with max intensity
            maxx = np.argmax(pks[1][tmpinds])
            ## make it so it isn't in vector of peaks to delete
            todel = np.delete(tmpinds, [maxx])
            inds_to_del.extend(todel)

    out = np.delete(pks, [inds_to_del], axis = 1)
    return out

def cosine_score(pks1, pks2):
    '''
    :param pks1: 2d array. First element is m/z and second is abund. Centroided
    :param pks2: Same
    :returns: float, cosine score
    '''
    comb_pks = []
    for index, x in enumerate(pks1[0]):
        temp = {"mass": x, "a": np.sqrt(pks1[1][index]) * np.square(pks1[0][index]), "b": 0}
        comb_pks.append(temp)
    
    for index, y in enumerate(pks2[0]):
        massExists = False

        for z in comb_pks:
            if abs(y - z['mass']) <= 0.01:
                z["b"] = np.sqrt(pks2[1][index]) * np.square(pks2[0][index])
                massExists = True

        if not massExists:
            temp = {"mass": y, "a": 0, "b": np.sqrt(pks2[1][index]) * np.square(pks2[0][index])}
            comb_pks.append(temp)
            massExists = False

    comb_pks.sort(key=lambda x: x["mass"])

    pkvectora = np.array([x['a'] for x in comb_pks], dtype = 'float')
    pkvectorb = np.array([x['b'] for x in comb_pks], dtype = 'float')

    top = np.square(np.sum(pkvectora * pkvectorb))
    bottom = np.sum(np.square(pkvectora)) * np.sum(np.square(pkvectorb))
    cosine = top/bottom
    return(cosine)

def get_highest_cosine_per_jgi_peak(jgi_peaks, our_peaks, mzthresh, rtthresh):
    '''
    :param jgi_peaks: dict of jgi's mgf
    :param our_peaks: dict of our mgf
    :param thresh: float, threshold for MS1 delta mz between jgi and our pk
    :returns: list of lists, each sublist has: jgi ID, our ID, cosine score
    '''
    outdf = pd.DataFrame(columns = ['JGI_ID', 'our_id', 'cosine'])
    plotdf = pd.DataFrame(columns = ['cosine', 'mz'])
    our_peak_df = pd.DataFrame.from_dict(our_peaks, orient = 'index')
    our_peak_df.columns = ['mz', 'rt', 'fragment_mz', 'fragment_intensity']

    for k, v in jgi_peaks.items():
        jgi_mass = v[0]
        jgi_rt = v[1]
        tmpdf = our_peak_df[(round(abs(our_peak_df['mz'] - jgi_mass), 5) <= mzthresh) &
        (round(abs(our_peak_df['rt'] - jgi_rt), 5) <= rtthresh)]

        if tmpdf.shape[0] > 0:
            jgi_centroid = centroid_peaks(np.row_stack((v[2], v[3])))
            tmp_results = pd.DataFrame(columns = ['id', 'cosine', 'mz'])
            for i in range(tmpdf.shape[0]):
                ours_centroid = centroid_peaks(np.row_stack((tmpdf.iloc[i, 2], tmpdf.iloc[i, 3])))
                cosine = cosine_score(jgi_centroid, ours_centroid)
                tmpseries = pd.Series([tmpdf.index[i], cosine, tmpdf.iloc[i, 0]], index = tmp_results.columns)
                tmp_results = tmp_results.append(tmpseries, ignore_index = True)
            try:
                max_ind = tmp_results['cosine'].idxmax()
                to_append = pd.Series([k, tmp_results.iloc[max_ind, 0], tmp_results.iloc[max_ind, 1]], 
                index = outdf.columns)
                plot_to_append = pd.Series([tmp_results.iloc[max_ind, 1], tmp_results.iloc[max_ind, 2]], 
                index = plotdf.columns)
                plotdf = plotdf.append(plot_to_append, ignore_index = True)
                outdf = outdf.append(to_append, ignore_index = True)
            except ValueError:
                print(tmp_results)
        
    return (outdf, plotdf)

def main(jgi_mgf, our_mgf, output, mzthresh, rtthresh):
    '''
    All inputs are strings
    '''
    jgid = read_in_mgf(jgi_mgf, 'jgi')
    ourd = read_in_mgf(our_mgf, 'ours')

    outdf, plotdf = get_highest_cosine_per_jgi_peak(jgid, ourd, float(mzthresh), float(rtthresh))
    outdf.to_csv(output + f'jgi_ourdata_matches_RT{rtthresh}_ms1{mzthresh}.tsv', sep = '\t')
    scatfig = sns.scatterplot(data = plotdf, x = "mz", y = "cosine").set_title(f'M/Z by '\
        f'cosine score with MS1 threshold {mzthresh}\n and RT threshold '\
            f'{rtthresh} and MS2 threshold 0.01')
    fig = scatfig.get_figure()
    fig.savefig(output + f'_cosine_by_mz_plot_RT{rtthresh}_ms1{mzthresh}.pdf')

if __name__ == '__main__': 

    if len(sys.argv) != 6:
        sys.exit('ARGS: 1) JGI mgf 2) Our corresponding mgf\n'
        '3) full path and prefix for output files 4) Threshold for MS1 mass error\n'
        '5) Threshold for max RT error')

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])

    print('Done!')
       