'''
primrose_samosa_ctcf_integration_OS.py
Vijay Ramani

This script takes CpG run length encodings provided by primrose that have been converted into a pickle format by 05_cpg2pickle.py
and computes CpG signal enrichment at CTCF sites as is done with m6dA.
'''

import os,sys,re
import pandas as pd
import pickle
import numpy as np
import scipy as sp


def parse_args():
    parser = argparse.ArgumentParser(description='Computes CpG signal enrichment at CTCF sites as is done with m6dA')
    parser.add_argument('-c','--cpg-pickle', nargs=1,dest='cpg_pickle', help='Path to CpG pickle file produced by 05_cpg2pickle.py')
    parser.add_argument('-s','--site-list',nargs=1,dest='site_list',help='TSV-formatted site list ',required=True)
    parser.add_argument('-o','--output-dir',nargs=1,dest='output_directory',help='Output directory',required=True)
    parser.add_argument('-p','--project',nargs=1,dest='project',help='Name of the project',required=True)

    args = parser.parse_args()
    return args


def eat_pickle_binary(pick):
    with open(pick, 'rb') as fout:
        tipds = pickle.load(fout, encoding="latin1")    
    return (tipds)

def distill_tipds_flat_density(tipds, sitelist, win_len, label, subsample=0):
    sites = open(sitelist,'r')
    hole_nos = {}
    for i in sites:
        split = i.split()
        #check to make sure the feature is in the middle of the read 
        #and enough on both sides
        if int(split[2]) <= int(split[4]) <= int(split[3]):
            index1 = int(split[4]) - int(split[2])
            index2 = int(split[3]) - int(split[4])
            strand = split[-1]
            r_strand = split[-2]
            if index1 > (win_len / 2) and index2 > (win_len / 2):
                hole_nos[split[0]] = (index1, strand, r_strand, (split[1],split[4]))
    sites.close()
    lengths = []
    labs = []
    reads = []
    densities = []
    fracs = []
    sites = []
    mno = 0
    for hole in hole_nos:
        prefix = 'm64182_220430_012143/' ## OS152 cell prefix
        suffix = '/ccs'
        key = prefix + hole + suffix
        if key not in tipds: 
            print("This happens\t%s" % hole)
            continue
        if len(tipds[key]) == 0: continue
        vecs = tipds[key][0]
        length = tipds[key][1]
        num_cpgs = len(vecs)
        read = np.empty(length)
        read[:] = np.nan
        for idx, val in vecs:
            read[idx] = val        
        index,strand,r_strand, site_info = hole_nos[hole]
        if r_strand == '-':
            read = read[::-1]
        if strand == "+":
            extract = read[int(index - (win_len / 2)): int(index + (win_len / 2))]
            if len(extract) != win_len: continue
            reads.append(extract)
            labs.append(hole)
            lengths.append(length)
            sites.append(site_info)
        else:
            extract = read[::-1][int(index - (win_len / 2)):int(index + (win_len / 2))]
            if len(extract) != win_len: continue
            reads.append(extract)
            labs.append(hole)
            lengths.append(length)
            sites.append(site_info)
        if subsample != 0:
            mno += 1
            if mno == subsample: break
    new_mat = np.vstack(reads)
    return (new_mat, np.array(lengths), np.array(labs), sites)

def main():

    args=parse_args()

    ## a list of ZMWs that overlap with CTCF sites 
    sitelist = args.site_list
    meth_pick=args.cpg_pickle
    
    length = 750
    #factors = ['OS152']    
    factors=[args.project]

    mats = []
    labels_tot = []
    densities_tot = []
    fracs_tot = []
    sites_tot = []    
    rep = eat_pickle_binary(meth_pick)
    mat, lengths, labels, sites = distill_tipds_flat_density(rep, sitelist, length, factors[0])
    print(sites)
    #mat = mat > 0.5
    chrs_tot = []
    sites_tot = []
    for chrid, site in sites:
        chrs_tot.append(chrid)
        sites_tot.append(site)
    print(len(labels))
    print(len(chrs_tot))
    print(len(sites_tot))
    #Save densities, fractions, and labels as dataframe
    mdata = pd.DataFrame({'labels_agg':labels, 
                          'chrid':chrs_tot, 'sites':sites_tot})

    mdata.to_csv('{}/{}_CTCF_meth_ctcf_final.data'.format(args.output_directory,args.project),sep=',')
    #mdata_labels = pd.DataFrame(columns=['sample','bio_rep','tech_rep','factor'])
    #mdata_labels[['sample','bio_rep','tech_rep','factor']] = mdata['labels_agg'].str.split('_',expand=True,n=3)
    #mdata_labels.to_csv('OS_CTCF_access_ctcf_sep_labels.data',sep=',')
    #Save mols matrix for plotting, clustering, etc.
    np.save('{}/{}_Ctcf_methmols_{}bp.npy'.format(args.output_directory,args.project,length), mat)
    
if __name__ == "__main__":
    main()


