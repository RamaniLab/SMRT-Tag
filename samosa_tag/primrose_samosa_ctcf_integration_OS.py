import os,sys,re
import pandas as pd
import pickle
import numpy as np
import scipy as sp

#This script takes CpG run length encodings provided by primrose that have been converted into a pickle format
#and computes signal enrichment at CTCF sites as with m6dA


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
        prefix = 'm64182_220430_012143/'
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
    sitelist = '/avicenna/vramani/analyses/pacbio/OS_SMRT_tag/symlink_for_aligned_files/U2OS_total_ZMWs.5kb.zmws'
    meth_pick = '/avicenna/vramani/analyses/pacbio/Tn5_OS_total_data.pickle'
    
    length = 750
    factors = ['OS152']    
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
    mdata.to_csv('OS_CTCF_meth_ctcf_final.data',sep=',')
    #mdata_labels = pd.DataFrame(columns=['sample','bio_rep','tech_rep','factor'])
    #mdata_labels[['sample','bio_rep','tech_rep','factor']] = mdata['labels_agg'].str.split('_',expand=True,n=3)
    #mdata_labels.to_csv('OS_CTCF_access_ctcf_sep_labels.data',sep=',')
    #Save mols matrix for plotting, clustering, etc.
    np.save('OS_Ctcf_methmols_%sbp.npy' % length, mat)
    
if __name__ == "__main__":
    main()


