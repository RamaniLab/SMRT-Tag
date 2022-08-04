'''
04_HMM.py
Colin McNally
2021/07/13

For each sample, load output of NN inference and binarize into accessibility footprints using HMM.
'''
import numpy as np
import pandas as pd
import socket
import sys
import os
import pickle
from pomegranate import HiddenMarkovModel, State, BernoulliDistribution
import glob
import re
from tqdm import tqdm
import multiprocessing as mp
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Binarize NN infernece outputs into accessibility footprints.")
    parser.add_argument('samples', nargs='+',help='Either an integer index into the reference file, or a string identifier in the format [cell].[samplename]')
    parser.add_argument('-t', '--threshold', type=float, help='Threshold for defining IPD residuals as methylated. [default] 0.42',default=0.42)
    parser.add_argument('-j', '--threads', type=int, help='Number of threads to use. Defaults to 1')
    parser.add_argument('-r','--referenceFile', nargs=1, required=True, help='The sample reference file from which to load sample information')
    parser.add_argument('-o', '--outputlocation',dest='output_directory',help='Location to save the outputs')

    args = parser.parse_args()
    return args

def runHMM(infile):

    ## Load the forHMM file
    with open(infile,'rb') as fin:
        hmmInput = pickle.load(fin)
    
    goodZMW = list(hmmInput.keys())
    hmmRes = {}
    for zmw in tqdm(goodZMW):
        goodDat = hmmInput[zmw]['inDat']
        cclen = hmmInput[zmw]['cclen']

        if goodDat.shape[0] < 2:
            continue
            # Don't try to create HMM if the input data is for some reason super short

        AcAd = []
        InacAd = []

        for row in goodDat.itertuples(index=True):
            
            AcAd.append(State(BernoulliDistribution(row.posProb), name="Ac_{0}".format(row.Index)))
            InacAd.append(State(BernoulliDistribution(row.negProb), name="Inac_{0}".format(row.Index)))

        model = HiddenMarkovModel()
        model.add_states(AcAd)
        model.add_states(InacAd)
        model.add_transition(model.start, AcAd[0], 0.5)
        model.add_transition(model.start, InacAd[0], 0.5)

        leaveInacProb = 1/1000
        leaveAcProb = 1/1000
        for b in np.arange(goodDat.shape[0]-1):
            dist = goodDat['pos'][b+1] - goodDat['pos'][b]
            stayInacP = (1 - leaveInacProb)**dist
            model.add_transition(InacAd[b], InacAd[b+1], stayInacP)
            model.add_transition(InacAd[b], AcAd[b+1], 1 - stayInacP)
            stayAcP = (1 - leaveAcProb)**dist
            model.add_transition(AcAd[b], AcAd[b+1], stayAcP)
            model.add_transition(AcAd[b], InacAd[b+1], 1 - stayAcP)

        model.add_transition(AcAd[-1], model.end, 1)
        model.add_transition(InacAd[-1], model.end, 1)
        model.bake()

        path = model.viterbi(goodDat['methPred'])
        refbases = np.arange(cclen)
        pathRes = np.full(goodDat.shape[0], np.nan)
        for p in path[1]:
            psplit = p[1].name.split('_')
            if len(psplit) > 1:
                if psplit[0] == 'Ac':
                    pathRes[int(psplit[1])] = 1
                if psplit[0] == 'Inac':
                    pathRes[int(psplit[1])] = 0

        hmmRes[zmw] = np.interp(refbases, goodDat['pos'], pathRes, left=np.nan, right=np.nan)
    
    return hmmRes
 


def main():
    args=parse_args()

    ## read sampleRef
    sampleRef = pd.read_csv(args.referenceFile,sep=',',index_col='index')

    ## set output directory
    dataPBase=args.output_directory

    ## for each sample
    for samp in args.samples:
        ## for each piece produced by NN inference
        pieces = sorted(glob(os.path.join(dataPBase,'processed','forHMM',
                '{0}_{1}_forHMM*resid-{2}*piece*.pickle'.format(
                                                        sampleRef['cell'][samp],
                                                        sampleRef['sampleName'][samp],
                                                        args.threshold
                                                        )
                                                        )))
       
        with mp.Pool(processes=args.threads) as pool:
            _hmmRes=pool.map(runHMM,pieces)

        ## merge dictionaries
        hmmRes={}
        for d in _hmmRes:
            hmmRes.update(d)

        if not os.path.exists(dataPBase + '/processed/binarized/HMMout/'):
            os.makedirs(dataPBase + '/processed/binarized/HMMout/')


        with open(os.path.join(dataPBase, 'processed','binarized','HMMout',"{}_{}_NNsingle.pickle".format(sampleRef['cell'][samp],sampleRef['sampleName'][samp])),'wb') as fout:
            pickle.dump(hmmRes, fout, protocol=4)

if __name__=='__main__':
    main()