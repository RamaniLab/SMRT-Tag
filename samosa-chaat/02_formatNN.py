'''
02_formatNN.py
Colin McNally
2021/07/13

For each sample,load IPD pickle files and process them into a format that can be used as input for NN inference
'''
import numpy as np
import pandas as pd
from tqdm import tqdm
import pickle
import os
import sys
import socket
import argparse
import multiprocessing as mp
import itertools
import functools
import glob
import re
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Produce data formatted for NN inference from IPD pickle files")
    parser.add_argument('samples', nargs='+',help='An integer index into the reference file')
    parser.add_argument('-r','--referenceFile', nargs=1, required=True, help='The sample reference file from which to load sample information')
    parser.add_argument('-o', '--outputlocation',dest='output_directory',help='Location to save the outputs')

    args = parser.parse_args()
    return args

def reverse_complement(sequence):
    revbase = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    newSequence = ''
    for b in sequence[::-1]:
        newSequence += revbase[b]
    return newSequence

def complement(sequence):
    revbase = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    newSequence = ''
    for b in sequence:
        newSequence += revbase[b]
    return newSequence


def createNNinputs(infile,samp, dataPBase, zmwlimit=None):
    bases = ['A', 'C', 'G', 'T']
    checkbases = np.concatenate([np.arange(-3,0), np.arange(1,11)])
    contextColDic = {checkbases[x]:{bases[i]:(x * 4 + i) for i in range(4)} for x in np.arange(len(checkbases))}
    usepercs = np.arange(10,91,10)

    ipdArr = [] # the mean IPD at this Adenine
    zmwArr = [] # The ZMW hole number
    zmwPos = [] # The position of this adenine within the CCS
    sampleArr = [] # The sample number of this zmw
    numSubArr = [] # The number of subreads that contributed to this adenine IPD mean

    with open(infile, 'rb') as fin:
        ipdfull = pickle.load(fin, encoding="latin1")
    
    zmwlist = list(ipdfull.keys())
    if zmwlimit != None:
        assert zmwlimit <= len(zmwlist), "There aren't that many ZMW in the sample!"
        zmwlist = zmwlist[0:zmwlimit]
    
   
    for zmw in zmwlist:    
        ipdd = ipdfull[zmw]

        for pos in np.arange(3, len(ipdd['read'])-10):                        
            if ~np.isnan(ipdd['reverseM'][pos]) and ipdd['read'][pos] == 'A':
                ipdArr.append(ipdd['reverseM'][pos])
                zmwArr.append(zmw)
                zmwPos.append(pos)
                sampleArr.append(samp)
                numSubArr.append(ipdd['reverseNsub'][pos])
        
        for pos in np.arange(10, len(ipdd['read'])-3):
            if ~np.isnan(ipdd['forwardM'][pos]) and ipdd['read'][pos] == 'T':
                ipdArr.append(ipdd['forwardM'][pos])
                zmwArr.append(zmw)
                zmwPos.append(pos)
                sampleArr.append(samp)
                numSubArr.append(ipdd['forwardNSub'][pos])
    
    
    ipdArr = np.array(ipdArr, dtype=np.float32).reshape(-1,1)
    zmwArr = np.array(zmwArr, dtype=np.int32).reshape(-1,1)
    zmwPos = np.array(zmwPos, dtype=np.int32).reshape(-1,1)
    sampleArr = np.array(sampleArr, dtype=np.int16).reshape(-1,1)
    numSubArr = np.array(numSubArr, dtype=np.int16).reshape(-1,1)
    
    contextmat = np.full((ipdArr.shape[0], len(checkbases)*4), False, dtype=np.bool_)
    percsmat = np.full((ipdArr.shape[0], 4 * len(usepercs)), np.nan, dtype=np.float32)
    
    ic = 0
    for zmw in zmwlist:
        ipdd = ipdfull[zmw]
        
        for pos in np.arange(3, len(ipdd['read'])-10):
            if ~np.isnan(ipdd['reverseM'][pos]) and ipdd['read'][pos] == 'A':
                revcontext = ipdd['read'][pos-3:pos+11]
                for offs in checkbases:
                    contextmat[ic, contextColDic[offs][revcontext[3+offs]]] = True
                for i in range(4):
                    percsmat[ic,np.arange(len(usepercs)) + (i * len(usepercs))] = ipdd['percentiles'][bases[i]][usepercs]
                ic += 1

        for pos in np.arange(10, len(ipdd['read'])-3):
            if ~np.isnan(ipdd['forwardM'][pos]) and ipdd['read'][pos] == 'T':
                forcontext = reverse_complement(ipdd['read'][pos-10:pos+4])
                for offs in checkbases:
                    contextmat[ic, contextColDic[offs][forcontext[3+offs]]] = True
                for i in range(4):
                    percsmat[ic,np.arange(len(usepercs)) + (i * len(usepercs))] = ipdd['percentiles'][bases[i]][usepercs]
                ic += 1
    
    if not os.path.exists(dataPBase + '/processed/forNN/'):
        os.makedirs(dataPBase + '/processed/forNN/')
        
    print("saving")

    np.savez(os.path.join(dataPBase, 'processed','forNN', os.path.basename(infile).replace('_full.pickle','_forNN.npz')), 
             contextmat = contextmat,
             ipdArr = ipdArr,
             contextColDic = contextColDic, 
             zmwArr = zmwArr, 
             zmwPos = zmwPos,
             sampleArr = sampleArr,
             percsmat = percsmat,
             numSubArr = numSubArr)
    


def main():
    args=parse_args()

    ## read sampleRef
    sampleRef = pd.read_csv(args.referenceFile,sep=',',index_col='index')

    ## set output directory
    dataPBase=args.output_directory

    ## for each sample
    for samp in args.samples:
        ## for each block of pickles
        parts = sorted(glob(os.path.join(dataPBase,'processed','full','{}_{}*_full.pickle'.format(sampleRef['cell'][samp],+sampleRef['sampleName'][samp]))))
        for part in parts:
            createNNinputs(infile=part,samp=samp, dataPBase=dataPBase)

if __name__=='__main__':
    main()