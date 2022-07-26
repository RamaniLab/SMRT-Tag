import os,sys,re
import numpy as np
import scipy as sp
from scipy import stats
import pandas as pd

#Same fishers exact script, but subsets the filtered data on the basis of replicate to test
#reproducibility of odds ratios across each technical replicate

fhi = open('test_cluster_cpg.dataout.filtered')
rep_labels = open('/avicenna/vramani/analyses/pacbio/OS_SMRT_tag/invivo_final_autocor_clusters_Tn5_2.data')
reps = {}
for line in rep_labels:
    split = line.split()
    rep, zmw = split[0].split('_')
    reps[zmw] = rep
clusters = []
meth_level = []
cpgs = []
rep_total = []
for line in fhi:
    split = line.split()
    zmw = split[0].split('/')[1]
    rep_total.append(reps[zmw])
    clusters.append(split[1])
    meth_level.append(float(split[3]))
    cpgs.append(float(split[2]))

cpg_df = pd.DataFrame({'reps':rep_total, 'clusters':clusters,'meth':meth_level,'cpgs':cpgs})
cpg_df['cpg_high'] = cpg_df['cpgs'].values > 10
cpg_df['meth_high'] = cpg_df['meth'].values > 0.5
fho = open('OS_fishers_by_rep.data','w')
for i in np.unique(cpg_df['cpg_high']):
    for j in np.unique(cpg_df['meth_high']):
        for k in np.unique(cpg_df['clusters']):
            for rep in np.unique(cpg_df['reps']):
                subSamp = cpg_df[cpg_df['reps'] == rep]
                num_clust_lab = len(subSamp[(subSamp['clusters'] == k) & (subSamp['cpg_high'] == i) & (subSamp['meth_high'] == j) ].values)
                num_clust_notlab = len(subSamp[(subSamp['clusters'] == k) & ((subSamp['cpg_high'] != i) | (subSamp['meth_high'] != j))].values)
                num_notclust_lab = len(subSamp[(subSamp['clusters'] != k) & (subSamp['cpg_high'] == i) & (subSamp['meth_high'] == j)].values)
                num_notclust_notlab = len(subSamp[(subSamp['clusters'] != k) & ((subSamp['cpg_high'] != i) | (subSamp['meth_high'] != j))].values)
                odds_r, pval = sp.stats.fisher_exact([[num_clust_lab, num_notclust_lab], \
                    [num_clust_notlab, num_notclust_notlab]])
                print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (i, j, k, num_clust_lab, odds_r, pval,rep), file=fho)
fho.close()
