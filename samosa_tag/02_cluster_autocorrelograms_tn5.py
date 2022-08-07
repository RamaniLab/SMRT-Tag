#!/usr/bin/env python3
'''
02_cluster_autocorrelograms_tn5.py
Vijay Ramani

This script takes a series of npy autocorrelograms calculated from the pickle HMM files produced by the SAMOSA-ChAAT pipeline
and clusters them using leiden clustering (res = 0.4, n_neighbors = 10). The script also filters out all molecules in clusters
that altogether account for < 10% of the data.
'''

import argparse
import os,sys,re
import numpy as np
import pandas as pd
import pickle
import scanpy

def parse_args():
    parser = argparse.ArgumentParser(description='Computes per molecule autocorrelograms out to 500 nucleotides (filter out molecules < 1 kb). Takes processed HMM pickles.')
    parser.add_argument('hmm_pickles', nargs='+', help='HMM pickle files')
    parser.add_argument('-a','--autocor-dir',nargs=1,dest='autocor_directory',help='Autocorrelation directory, (output-dir from 01_compute_autocors_persample.py)',required=True)
    parser.add_argument('-o','--output-dir',nargs=1,dest='output_directory',help='Output directory',required=True)
    parser.add_argument('-p','--project',nargs=1,dest='project',help='Project',required=True)

    args = parser.parse_args()
    return args

def process_samps(tipds, min_len=1000):
    zmws = []
    for zmw in tipds:
        read = np.nan_to_num(tipds[zmw]) >= 0.5
        if len(read) < min_len: continue
        zmws.append(read[:1000])
    return(np.vstack(zmws))

def eat_pickle_binary(pick):
    with open(pick, 'rb') as fout:
        tipds = pickle.load(fout, encoding="latin1")    
    return (tipds)

def zmw_filter(tipds, zmw_list,samp_label):
    '''We are only going to cluster ZMWs that we calculated autocorrelograms for.'''
    list_open = open(zmw_list)
    signal = []
    zmws = []
    for line in list_open:
        split = line.split()
        zmw = int(split[0])
        zmws.append("%s_%s" % (samp_label, zmw))
        read = np.nan_to_num(tipds[zmw]) >= 0.5
        signal.append(read[:1000])
    signal = np.vstack(signal)
    list_open.close()
    zmws = np.array(zmws)
    return (signal, zmws)

    
def cluster_autocorrelograms(auto_cor):
    np.nan_to_num(auto_cor, copy=False)
    auto_cor = scanpy.AnnData(X=auto_cor)
    scanpy.tl.pca(auto_cor)
    scanpy.pp.neighbors(auto_cor, metric='correlation', n_neighbors=10)
    scanpy.tl.leiden(auto_cor, resolution=0.4)
    clusters = np.array(auto_cor.obs['leiden'])
    return(clusters)

def cluster_cutoff(clusters, cut=0.95):
    tot = len(clusters)
    cum = 0
    for i in range(len(np.unique(clusters))):
        frac = len(clusters[clusters == str(i)]) / tot
        cum += frac
        if cum >= cut:
            cutoff = i
            break
    return cutoff
        


def main():
    args=parse_args()

    filelist_signal_npy=[]
    filelist_npy=[]
    filelist_txt=[]

    for hmm_pickle in args.hmm_pickles:
        fname = os.path.basename(hmm_pickle)
        filelist_signal_npy.append(hmm_pickle)
        filelist_npy.append("{}/{}.autocors.npy".format(args.autocor_directory,fname))
        filelist_txt.append("{}/{}.zmw_ids.txt".format(args.autocor_directory,fname))

    #file_labels = ['Rep1','Rep2','Rep3','Rep4','Rep5','Rep6','Rep7','Rep8'...]
    file_labels = ['Rep{}'.format(i) for i in range(len(args.hmm_pickles))]
    
    file_label_lookup = zip(filelist_npy, filelist_signal_npy, filelist_txt, file_labels)

    zmw_tot = []
    samps_tot = []
    arrs = []
    arrs_sig = []

    for a,b,c,d in file_label_lookup:
        file = a #pwd_hmm + a
        zmws = c #pwd_hmm + c
        signal_file = eat_pickle_binary(b)
        autocor = np.load(file)
        arrs.append(autocor)
        signal, zmw_list = zmw_filter(signal_file, zmws, d)
        zmw_tot.append(zmw_list)
        arrs_sig.append(signal)

    arrs = np.vstack(arrs)
    arrs_sig = np.vstack(arrs_sig)
    zmw_tot = np.concatenate(zmw_tot)

    print(len(arrs))
    print(len(arrs_sig))

    clusters = cluster_autocorrelograms(arrs)
    #filter molecules (we'll say all clusters that collectively make up <5% of the data)
    cutoff = cluster_cutoff(clusters,cut=0.9)

    fho1 = open('{}/{}_autocor_clusters.data'.format(args.output_directory,args.project), 'w')
    fho2 = open('{}/{}_autocor_averages.data'.format(args.output_directory,args.project), 'w')
    fho3 = open('{}/{}_autocor_signal_averages.data'.format(args.output_directory,args.project), 'w')
    for i in range(len(zmw_tot)):
        if int(clusters[i]) > cutoff: continue #use cutoff to filter out molecules from one of the filtered clusters
        print("%s\t%s" % (zmw_tot[i], clusters[i]), file = fho1)

    for i in np.unique(clusters):
        if int(i) > cutoff: continue #use cutoff to filter out molecules from one of the filtered clusters
        avg = np.nanmean(arrs[clusters == i],axis=0)
        avg_sig = np.nanmean(arrs_sig[clusters == i], axis=0)
        for idx in range(len(avg)):
            print("%s\t%s\t%s" % (idx, avg[idx], i), file = fho2)
        for idx in range(len(avg_sig)):
            print("%s\t%s\t%s" % (idx, avg_sig[idx], i), file = fho3)

    fho1.close()
    fho2.close()
    fho3.close()
    
if __name__ == "__main__":
    main()
    

