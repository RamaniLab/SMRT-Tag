import os,sys,re
import numpy as np
import scipy as sp
from scipy import stats
import pandas as pd

#Takes the list of filtered clusters and computes Fisher's exact test odds ratios
#and p-values to determine fiber type enrichment and depletion in differentially
#methylated genomic regions.


fhi = open('test_cluster_cpg.dataout.filtered')
clusters = []
meth_level = []
cpgs = []
for line in fhi:
    split = line.split()
    clusters.append(split[1])
    meth_level.append(float(split[3]))
    cpgs.append(float(split[2]))

cpg_df = pd.DataFrame({'clusters':clusters,'meth':meth_level,'cpgs':cpgs})
cpg_df['cpg_high'] = cpg_df['cpgs'].values > 10
cpg_df['meth_high'] = cpg_df['meth'].values > 0.5
fho = open('OS_fishers.data','w')
for i in np.unique(cpg_df['cpg_high']):
    for j in np.unique(cpg_df['meth_high']):
        for k in np.unique(cpg_df['clusters']):
            subSamp = cpg_df
            num_clust_lab = len(subSamp[(subSamp['clusters'] == k) & (subSamp['cpg_high'] == i) & (subSamp['meth_high'] == j) ].values)
            num_clust_notlab = len(subSamp[(subSamp['clusters'] == k) & ((subSamp['cpg_high'] != i) | (subSamp['meth_high'] != j))].values)
            num_notclust_lab = len(subSamp[(subSamp['clusters'] != k) & (subSamp['cpg_high'] == i) & (subSamp['meth_high'] == j)].values)
            num_notclust_notlab = len(subSamp[(subSamp['clusters'] != k) & ((subSamp['cpg_high'] != i) | (subSamp['meth_high'] != j))].values)
            odds_r, pval = sp.stats.fisher_exact([[num_clust_lab, num_notclust_lab], \
                [num_clust_notlab, num_notclust_notlab]])
            print("%s\t%s\t%s\t%s\t%s\t%s" % (i, j, k, num_clust_lab, odds_r, pval), file=fho)
fho.close()
