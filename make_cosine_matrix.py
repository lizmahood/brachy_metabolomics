import numpy as np
import pandas as pd
import ray, sys

def read_in_mgf(mgfin):
    '''
    :param mgfin: string, path to input
    '''
    odict = {}
    with open(mgfin, 'r') as mgff:
        lin = mgff.readline()
        while lin:
            if lin.startswith('SCAN'):
                name = int(lin.strip().split('=')[1])
                odict[name] = []
                peaks = []
                abunds = []
            elif lin.startswith('PEP'):
                pm = float(lin.strip().split('=')[1])
                odict[name].append(pm)
            elif lin.startswith('MSL'):
                additional = lin + mgff.readline()
                odict[name].append(additional)
            elif lin.startswith('RT'):
                rt = float(lin.strip().split('=')[1])
                odict[name].append(rt)
            elif lin.startswith('ION'):
                ion = lin.strip().split('=')[1]
                odict[name].append(ion)
            elif lin[0].isdigit():
                abund = float(lin.strip().split()[1])
                #if abund >= 1500:
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

def cosine_score_hybrid(pks1, pks2):
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

def find_pk_intersection(pks1, pks2):
    '''
    :param pks1: 2d array. First element is m/z and second is abund. Centroided
    :param pks2: Same
    :returns: float 
    '''
    wts1 = []
    wts2 = []
    ## making lists of pks1 pks/intensities, followed by pks2
    masses1l = list(pks1[0])
    masses1l.extend(pks2[0])

    intens1l = list(pks1[1])
    intens1l.extend(pks2[1])

    ancestors = (['a'] * len(pks1[0])) + (['b'] * len(pks2[0]))
    ## zipping and sorting these
    alltogether = zip(masses1l, intens1l, ancestors)
    res = sorted(alltogether, key = lambda x: x[0])

    for indd in range(len(res) -1):
        if res[indd + 1][0] - res[indd][0] <= 0.01:
            if res[indd + 1][2] != res[indd][2]:
                ## Got a match! Getting the weights
                thiswt = np.square(res[indd][0]) * np.sqrt(res[indd][1])
                nextwt = np.square(res[indd + 1][0]) * np.sqrt(res[indd + 1][1])

                if res[indd][2] == 'a':
                    wts1.append(thiswt)
                    wts2.append(nextwt)
                else:
                    wts2.append(thiswt)
                    wts1.append(nextwt) 
    
    if len(wts1) == len(wts2) and len(wts1) > 0:
        wts1 = np.asarray(wts1)
        wts2 = np.asarray(wts2)
        top = np.square(np.sum(wts1 * wts2))
    elif len(wts1) == 0:
        top = 0
    elif len(wts1) != len(wts2): 
        print('Uh-oh!')
    return top

def cosine_score_li(pks1, pks2):
    '''
    Calculates the cosine score according to Li et al Science Advances 2020
    '''
    top = find_pk_intersection(pks1, pks2)
    weights1 = np.square(pks1[0]) * np.sqrt(pks1[1])
    weights2 = np.square(pks2[0]) * np.sqrt(pks2[1])
    bottom = np.sum(np.square(weights1)) * np.sum(np.square(weights2))

    return top/bottom


#@ray.remote
def get_scores(k1, mgfd):
    '''
    :param k1: str
    :param mgfd: dict
    :returns: list of cosine scores of k1 vs all entries above it
    '''
    keyl = list(mgfd.keys())
    outl = []
    query = np.row_stack((np.array(mgfd[k1][-2],dtype = 'float'), 
                            np.array(mgfd[k1][-1], dtype = 'float')))
    querycent = centroid_peaks(query)

    ind = keyl.index(k1)

    for i in range((ind + 1), len(keyl)):
        k2 = keyl[i]
        matchcent = centroid_peaks(np.row_stack((np.array(mgfd[k2][-2],dtype = 'float'),
                                                    np.array(mgfd[k2][-1], dtype = 'float'))))
        cosine = cosine_score_li(querycent, matchcent)
        outl.append(cosine)

    ## Need to pad the beginning of the array with 0s
    ## To make outl the same length as mgfd.keys()
    nums_0 = len(list(mgfd.keys())) - len(outl)
    list_0 = [0] * nums_0
    list_0.extend(outl)
    return list_0

def main(mgfdin, outp):
    '''
    '''
    mgfd = read_in_mgf(mgfdin)

    ## I don't think this is going to work
    #ray.init()

    ## Many more things to loop through than there are cores
    futures = []
    for index, k1 in enumerate(mgfd.keys()):
        futures.append(get_scores(k1, mgfd))
        if index % 100 == 0:
            print(index)
    
    #hmm = ray.get(futures)

    ## Turning the list of lists into ndarray with upper and lower tris
    utri = np.array(futures)
    both_tri = utri + utri.T - np.diag(np.diag(utri))
    np.fill_diagonal(both_tri, 1)
    df = pd.DataFrame(both_tri, columns = list(mgfd.keys()), index = list(mgfd.keys()))

    df.to_csv(outp, sep = '\t', index = mgfd.keys())

if __name__ == '__main__':

    if len(sys.argv) != 3:
        sys.exit('ARGS: 1) Path to input mgf 2) Full path and full name of output matrix')

    main(sys.argv[1], sys.argv[2])

    print('Done!')

