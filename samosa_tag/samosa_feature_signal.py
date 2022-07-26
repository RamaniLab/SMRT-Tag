import os,sys,re
import pandas as pd
import pickle
import numpy as np
import scipy as sp

#Compute signal enrichment at a specific feature (in this case CTCF sites)
#using output from zmw_selector.py. In this case, filepaths are hardcoded
#but can easily be changed.


def eat_pickle_binary(pick):
    with open(pick, 'rb') as fout:
        tipds = pickle.load(fout, encoding="latin1")    
    return (tipds)

def distill_tipds_flat_density(tipds, sitelist, length, label, subsample=0):
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
            if index1 > (length / 2) and index2 > (length / 2):
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
        if int(hole) not in tipds: 
#           print("This happens\t%s" % hole)
            continue
        read = tipds[int(hole)][12:-12]
        index,strand,r_strand, site_info = hole_nos[hole]
        if r_strand == '-':
            read = read[::-1]
        if strand == "+":
            extract = tipds[int(hole)][int(index - (length / 2)): int(index + (length / 2))]
            if len(extract) != length: continue
            reads.append(extract)
            labs.append(hole)
            lengths.append(len(tipds[int(hole)]))
            sites.append(site_info)
        else:
            extract = tipds[int(hole)][::-1][int(index - (length / 2)):int(index + (length / 2))]
            if len(extract) != length: continue
            reads.append(extract)
            labs.append(hole)
            lengths.append(len(tipds[int(hole)]))
            sites.append(site_info)
        if subsample != 0:
            mno += 1
            if mno == subsample: break
    new_mat = np.vstack(reads)
    return (new_mat, np.array(lengths), np.array(labs), sites)

def main():
    pwd_hmm = '/avicenna/vramani/analyses/pacbio/OS_SMRT_tag/symlink_for_hmm_files/'
    pwd_aligned = '/avicenna/vramani/analyses/pacbio/OS_SMRT_tag/symlink_for_aligned_files/'
    
    
    filelist_npy = ['OS152_OS152_PR27_plusM_NNsingle_HMM-v2.pickle',
                    'OS152_OS152_PR28_plusM_NNsingle_HMM-v2.pickle',
                    'OS152_OS152_PR29_plusM_NNsingle_HMM-v2.pickle',
                    'OS152_OS152_PR30_plusM_NNsingle_HMM-v2.pickle',
                    'OS152_OS152_PR31_plusM_NNsingle_HMM-v2.pickle',
                    'OS152_OS152_PR32_plusM_NNsingle_HMM-v2.pickle',
                    'OS152_OS152_PR33_plusM_NNsingle_HMM-v2.pickle',
                    'OS152_OS152_PR34_plusM_NNsingle_HMM-v2.pickle']

    filelist_txt = ['Tn5.A13.subread.ccs.OS152_PR27_plusM.aln.sorted.bam',
                    'Tn5.A13.subread.ccs.OS152_PR28_plusM.aln.sorted.bam',
                    'Tn5.A13.subread.ccs.OS152_PR29_plusM.aln.sorted.bam',
                    'Tn5.A13.subread.ccs.OS152_PR30_plusM.aln.sorted.bam',
                    'Tn5.A13.subread.ccs.OS152_PR31_plusM.aln.sorted.bam',
                    'Tn5.A13.subread.ccs.OS152_PR32_plusM.aln.sorted.bam',
                    'Tn5.A13.subread.ccs.OS152_PR33_plusM.aln.sorted.bam',
                    'Tn5.A13.subread.ccs.OS152_PR34_plusM.aln.sorted.bam']    
    
    file_labels = ['Rep1','Rep2','Rep3','Rep4','Rep5','Rep6','Rep7','Rep8']
    
    file_label_lookup = zip(filelist_npy, filelist_txt, file_labels)
    length = 750
    factors = ['U2OS']    
    mats = []
    labels_tot = []
    densities_tot = []
    fracs_tot = []
    sites_tot = []
    for signal, aligned, rep_label in file_label_lookup:
        signal_open = "%s/%s" % (pwd_hmm, signal)
        rep = eat_pickle_binary(signal_open)
        for factor in factors:
            site_open = "%s/%s.%s.5kb.zmws" % (pwd_aligned, factor, aligned)
            mat_rep1, lengths, labels, sites = distill_tipds_flat_density(rep, site_open, length, "%s_%s" % (rep_label, factor))
            mats.append(mat_rep1)
            labels_tot.append(labels)
            sites_tot.append(sites)
        rep = 0
    labels_tot = np.concatenate(labels_tot)
    sites_tot = np.concatenate(sites_tot)
    #print(sites_tot)
    mats = np.vstack(mats) > 0.5
    chrs = []
    sites = []
    for chrid, site in sites_tot:
        chrs.append(chrid)
        sites.append(site)
    #Save densities, fractions, and labels as dataframe
    mdata = pd.DataFrame({'labels_agg':labels_tot, 
                          'chrid':chrs, 'sites':sites})
    mdata.to_csv('OS_CTCF_access_ctcf_final_include_hole.data',sep=',')
    #mdata_labels = pd.DataFrame(columns=['sample','bio_rep','tech_rep','factor'])
    #mdata_labels[['sample','bio_rep','tech_rep','factor']] = mdata['labels_agg'].str.split('_',expand=True,n=3)
    #mdata_labels.to_csv('OS_CTCF_access_ctcf_sep_labels.data',sep=',')
    #Save mols matrix for plotting, clustering, etc.
    np.save('OS_Ctcf_mols_%sbp.npy' % length, mats)
    
if __name__ == "__main__":
    main()


