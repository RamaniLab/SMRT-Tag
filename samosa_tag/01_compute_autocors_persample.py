#!/usr/bin/env python3
'''
01_compute_autocors_persample.py
Vijay Ramani

This script takes processed HMM pickles produced by the SAMOSA-ChAAT computational pipeline and computes
per molecule autocorrelograms out to 500 nucleotides (filter out molecules < 1 kb). 

'''

import os,sys,re
import numpy as np
import scipy as sp
import pickle
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Computes per molecule autocorrelograms out to 500 nucleotides (filter out molecules < 1 kb).')
    parser.add_argument('hmm_pickles', nargs='+', help='HMM pickle files')
    parser.add_argument('-o','--output-dir',nargs=1,dest='output_directory',help='Output directory',required=True)
    parser.add_argument('-p','--project',nargs=1,dest='project',help='Project',required=True)

    args = parser.parse_args()
    return args

def eat_pickle_binary(pick):
    with open(pick, 'rb') as fout:
        tipds = pickle.load(fout, encoding="latin1")    
    return (tipds)

def ccf(x, y):
    result = np.correlate(y - np.mean(y), x - np.mean(x), 'same') / (np.std(y) * np.std(x) * len(y))
    length = (len(result) - 1) // 2
    return result[length:]

def xcor(a,b, maxlen = 1000, lengths=False):
    foo = []
    for i in range(len(a)):
        res = ccf(a[i],b[i])
        foo.append(res)
    return foo

def process_xcors(tipds, min_len=1000, lengths=False):
    cors = []
    zmws = []

    for zmw in tipds:
        read = np.nan_to_num(tipds[zmw]) >= 0.5
        if len(read) < min_len: continue
        cor = ccf(read, read)[:500]
        cors.append(cor)
        zmws.append(zmw)
    return (np.vstack(cors), zmws)

def main():
    args = parse_args()
    for hmm_pickle in args.hmm_pickles:
        #hmm_pickle=_hmm_pickle.replace('.pickle','')
        print(hmm_pickle)

        binary = eat_pickle_binary(hmm_pickle)
        cors, zmws = process_xcors(binary)

        np.save('{}/{}.autocors.npy'.format(args.output_directory,hmm_pickle), cors)
        print("File written.")

        with open('{}/{}.zmw_ids.txt'.format(args.output_directory,hmm_pickle), 'w') as fho:
            for zmw in zmws:
                print(zmw, file = fho)
            print('ZMWs written.')


if __name__ == "__main__":
    main()
