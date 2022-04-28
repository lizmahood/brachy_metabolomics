###Hopefully this can be as general-purpose as possible

##Have file of alignment IDs (have user input column with align IDs)
##Have mgf (total mgf)
##Output filtered mgf of just compounds with an ID in the file
##Also possible that not all IDs have an entry in the mgf

##Need: from 4make_flavonoid_mgf.py import get_alignmentID_from_mgf, filter_ms2_peaks, write_selected_mgf
##Need to make different get_ids_from_file function.

##Use get_ids_from_file (return list of ids), then get_alignmentID_from_mgf (return dict of overall mgf),
##then filter_ms2_peaks (on whole mgfdict, returns filtered dict), then write_selected_mgf (feed it whole
##filtered mgf, and list of ids)

import sys
import numpy as np

def get_ids_of_interest(alignin, index):
    '''
    :param alignin: string, path to intput alignment file
    :param index: int, what is the index of the column with ids?
    returns: list of alignment ids (ints)
    '''
    aligndict = {}
    with open(alignin, 'r') as inf:
        inl = inf.readline().strip()
        while inl:
            inlist = inl.split('\t')
            if inlist[index][0].isdigit() and inlist[31] != '':
                k = int(inlist[index])
                add = inlist[4]
                if 'adduct' in inlist[5]:
                    addl = inlist[5].split(';')
                    connect = [x.split('adduct linked to ')[1].split('_')[0] for x in addl if 'adduct' in x]
                    connected = [int(x) for x in connect]
                    aligndict[k] = [add, addl, connected]
                else:
                    aligndict[k] = [add, 'None', 'None']
            inl = inf.readline().strip()

    return aligndict

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

def filter_out_adducts(aligndict, mgfd, threshold):
    '''
    :param aligndict: dict, output of get_align_ids_of_interest
    :param mgfd: dict, output of read_in_mgf
    :param threshold: float, above which 2 metabolites are adducts
    :returns: mgfd, but without adduct entries
    '''
    keystorem = []
    for k, v in aligndict.items():
        if aligndict[k][2] != 'None':
            if k in mgfd.keys():
                query = np.row_stack((np.array(mgfd[k][-2],dtype = 'float'), 
                                      np.array(mgfd[k][-1], dtype = 'float')))
                for i in aligndict[k][2]:
                    if i in mgfd.keys():
                        querycent = centroid_peaks(query)
                        adductcent = centroid_peaks(np.row_stack((np.array(mgfd[i][-2],dtype = 'float'),
                                                                  np.array(mgfd[i][-1], dtype = 'float'))))
                        cosine = cosine_score(querycent, adductcent)
                        if cosine >= threshold:
                            keystorem.append(i)

    new_dict = {key:val for key, val in mgfd.items() if key not in keystorem}
    return new_dict

def write_selected_mgf(newdict, aligndict, ofil):
    '''
    :param mgfdict: is the dict with values as mgf entries
    :param aligndict: dict of ids to output
    :param ofil: is a string of path to output dir
    '''
    with open(ofil, 'w') as out:
        for k, v in newdict.items():
            if k in aligndict.keys():
                towrite = f'BEGIN IONS\nSCANS={k}\nPEPMASS={v[0]}\n'\
                f'{v[1]}RTINMINUTES={v[2]}\nION={v[3]}\n'
                out.write(towrite)
                for pk, abund in zip(v[-2], v[-1]):
                    out.write(f'{pk} {abund}\n')
            
                out.write('END IONS\n\n')

def main(mgfin, alignin, index, threshold, ofil):
    '''
    All inputs strings except index
    '''
    aligndict = get_ids_of_interest(alignin, index)
    mgfd = read_in_mgf(mgfin)
    mgfdfilt = filter_out_adducts(aligndict, mgfd, threshold)
    write_selected_mgf(mgfdfilt, aligndict, ofil)

if __name__ == '__main__':

    if len(sys.argv) != 6:
        sys.exit('ARGS: 1) full mgf to select IDs from '\
            '2) file with IDs of interest 3) index of col of file(2) with IDs '\
                '4) threshold for adduct cosine score 5) full path to output directory')

    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), sys.argv[5])
    print('Done!')