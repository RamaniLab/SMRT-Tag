'''
05_postprocesing.py
Colin McNally
2021/07/13

For each sample,load accessibility footprints and call nucleosomes. 
'''
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
import pickle
import seaborn as sns
from glob import glob
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Load accessibility footprints and call nucleosomes")
    parser.add_argument('samples', nargs='+',help='Either an integer index into the reference file, or a string identifier in the format [cell].[samplename]')
    parser.add_argument('-r','--referenceFile', nargs=1, required=True, help='The sample reference file from which to load sample information')
    parser.add_argument('-o', '--outputlocation',dest='output_directory',help='Location to save the outputs')

    args = parser.parse_args()
    return args

def call_inaccessible_regions(dataPBase,samp,sampleRef):
    '''Call inacessible regions'''

    zmwinfoFile = '{0}/processed/full/{1}_{2}_full_zmwinfo.pickle'.format(dataPBase,
                                                                              sampleRef['cell'][samp],
                                                                              sampleRef['sampleName'][samp])
    zmwinfo = pd.read_pickle(zmwinfoFile)
    zmwinfo.set_index('zmw', inplace=True)
    zmwinfo.drop_duplicates(inplace=True) # remove duplicate entries


    hmmFile = '{0}/processed/HMMout/{1}_{2}_NNsingle.pickle'.format(dataPBase,
                                                                  sampleRef['cell'][samp],
                                                                  sampleRef['sampleName'][samp])
    with open(hmmFile, 'rb') as fin:
        hmmdat = pickle.load(fin)


    regions = {'zmw':[], 'length':[], 'start':[], 'end':[]}

    for zmw in tqdm(zmwinfo['zmw'], desc=str(samp), position=0, mininterval=1):
        
        ## try and load the ZMW
        try:
            hmm = hmmdat[zmw]
            goodzmw = True
        except KeyError:
            goodzmw = False
            pass
        
        
        if goodzmw:
            ## hard binarize the inaccessible regions
            inacregion = hmm
            inacregion[np.isfinite(inacregion)] = inacregion[np.isfinite(inacregion)] > 0.5

            ## locate switch regions
            inacswitch = np.diff(inacregion)
            switchp = np.where(np.logical_or(inacswitch == 1, inacswitch == -1))[0]

            if len(switchp) < 1:     
                if hmm[40] == 0: #checking if the entire molecule is inaccessible, 50 being an arbitrary position to check
                    regions['zmw'].append(zmw)
                    regions['length'].append(len(hmm) - 0)
                    regions['start'].append(np.nan)
                    regions['end'].append(np.nan)
                continue

            if inacswitch[switchp[0]] == -1:
                inInacReg = False
                regStart = -1
                regEnd = -1
            if inacswitch[switchp[0]] == 1:
                inInacReg = True
                regStart = np.nan
                regEnd = -1
                
                
            for point in switchp:
                if inacswitch[point] == -1 and not inInacReg:
                    inInacReg = True
                    regStart = point + 1
                if inacswitch[point] == 1 and inInacReg:
                    inInacReg = False
                    regEnd = point
                    regions['zmw'].append(zmw)
                    if np.isnan(regStart):
                        regions['length'].append(regEnd - 0)
                    else:
                        regions['length'].append(regEnd - regStart)
                    regions['start'].append(regStart)
                    regions['end'].append(regEnd)
            if inInacReg:
                regions['zmw'].append(zmw)
                regions['length'].append(len(hmm) - regStart)
                regions['start'].append(regStart)
                regions['end'].append(np.nan)
                

    regionD = pd.DataFrame(regions)

    if not os.path.exists(dataPBase + '/processed/inaccessibleRegions/'):
        os.makedirs(dataPBase + '/processed/inaccessibleRegions/')

    regionD.to_csv(dataPBase + '/processed/inaccessibleRegions/{0}_{1}_inacRegions.csv'.format(sampleRef['cell'][samp],sampleRef['sampleName'][samp]))


def calculate_density(dataPBase,samp,sampleRef):
    divs = np.array([130, 250, 370, 500, 630, 800, 970, 1090, 1200]) 


    zmwinfoFile = '{0}/processed/full/{1}_{2}_full_zmwinfo.pickle'.format(dataPBase,
                                                                              sampleRef['cell'][samp],
                                                                              sampleRef['sampleName'][samp])
    zmwinfo = pd.read_pickle(zmwinfoFile)
    zmwinfo.set_index('zmw', inplace=True)
    zmwinfo.drop_duplicates(inplace=True) # remove duplicate entries


    hmmFile = '{0}/processed/HMMout/{1}_{2}_NNsingle.pickle'.format(dataPBase,
                                                                  sampleRef['cell'][samp],
                                                                  sampleRef['sampleName'][samp])
    with open(hmmFile, 'rb') as fin:
        hmmdat = pickle.load(fin)


    regionFile = '{0}/processed/inaccessibleRegions/{1}_{2}_inacRegions.csv'.format(dataPBase,
                                                                                        sampleRef['cell'][samp],
                                                                                        sampleRef['sampleName'][samp])
    regiondf = pd.read_csv(regionFile, index_col=0)

    densityD = {'zmw':[], 'numRegions':[], 'nucs':[], 'subnucs':[], 'moleculeLength':[], 'calledLength':[], 'fracInaccessible':[], 'overnucs':[]}

    doneZMW = []
    curZMW = None
    for ir in np.arange(regiondf.shape[0]):
        if regiondf['zmw'][ir] != curZMW:
            if curZMW is not None:
                curNs = np.array(curNs)
                densityD['zmw'].append(curZMW)
                doneZMW.append(curZMW)
                densityD['subnucs'].append(np.sum(curNs == 0))
                densityD['overnucs'].append(np.sum(curNs == 8))   #consider changing this qualificiation for overnuc
                densityD['numRegions'].append(len(curNs))
                densityD['nucs'].append(np.sum(curNs[curNs < 8]))   
                densityD['moleculeLength'].append(zmwinfo['cclen'][curZMW])

                hmm = hmmdat[curZMW]
                densityD['calledLength'].append(np.sum(np.isfinite(hmm)))
                inacregion = hmm
                inacregion[np.isfinite(inacregion)] = inacregion[np.isfinite(inacregion)] > 0.5
                densityD['fracInaccessible'].append(np.sum(inacregion == 0) / np.sum(np.isfinite(inacregion)))

            curZMW = regiondf['zmw'][ir]
            curNs = []
        thislength = regiondf['length'][ir]
        thisnucs = np.sum(thislength > divs)
        curNs.append(thisnucs)
    # now do last molecule
    curNs = np.array(curNs)
    densityD['zmw'].append(curZMW)
    doneZMW.append(curZMW)
    densityD['subnucs'].append(np.sum(curNs == 0))
    densityD['overnucs'].append(np.sum(curNs == 9))
    densityD['numRegions'].append(len(curNs))
    densityD['nucs'].append(np.sum(curNs[curNs < 9]))
    densityD['moleculeLength'].append(zmwinfo['cclen'][curZMW])

    hmm = hmmdat[curZMW]
    densityD['calledLength'].append(np.sum(np.isfinite(hmm)))
    inacregion = hmm
    inacregion[np.isfinite(inacregion)] = inacregion[np.isfinite(inacregion)] > 0.5
    densityD['fracInaccessible'].append(np.sum(inacregion == 0) / np.sum(np.isfinite(inacregion)))
                
    # add any missing molecules that had zero inaccessible regions
    doneZMW = np.array(doneZMW)
    allzmw = np.array(list(hmmdat.keys())) #zmwinfo.index.to_numpy()
    missingZMW = allzmw[~np.isin(allzmw, doneZMW)]
    for uzmw in missingZMW:
        densityD['zmw'].append(uzmw)
        densityD['subnucs'].append(0)
        densityD['overnucs'].append(0)
        densityD['numRegions'].append(0)
        densityD['nucs'].append(0)
        densityD['moleculeLength'].append(zmwinfo['cclen'][uzmw])

        hmm = hmmdat[uzmw]
        densityD['calledLength'].append(np.sum(np.isfinite(hmm)))
        inacregion = hmm
        inacregion[np.isfinite(inacregion)] = inacregion[np.isfinite(inacregion)] > 0.5
        densityD['fracInaccessible'].append(np.sum(inacregion == 0) / np.sum(np.isfinite(inacregion)))
        
    densityDF = pd.DataFrame(densityD)
    densityDF.set_index('zmw', inplace=True)
    densityDF.sort_index(axis='index', inplace=True)
    
    if not os.path.exists(dataPBase + '/processed/density/'):
            os.makedirs(dataPBase + '/processed/density/')

    densityFile = dataPBase + '/processed/density/{0}_{1}_density.csv'.format(
                                                                            sampleRef['cell'][samp],
                                                                            sampleRef['sampleName'][samp])
    densityDF.to_csv(densityFile) #columns=['zmw', 'numRegions', 'nucs', 'subnucs', 'overnucs', 'moleculeLength', 'calledLength', 'fracInaccessible'])


def main():
    args=parse_args()

    ## read sampleRef
    sampleRef = pd.read_csv(args.referenceFile,sep=',',index_col='index')

    ## set output directory
    dataPBase=args.output_directory

    ## for each sample
    for samp in args.samples:
        call_inaccessible_regions(dataPBase,samp,sampleRef)
        calculate_density(dataPBase,samp,sampleRef)

if __name__=='__main__':
    main()