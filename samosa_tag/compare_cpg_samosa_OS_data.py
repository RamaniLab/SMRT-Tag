import os,sys,re
import pickle
import numpy as np

#Link CpG information and m6dA signal between ZMWs

def link_clusters_cpg(cluster_labels, cpg_pickle,fho):
    with open(cpg_pickle, 'rb') as fout:
        tipds = pickle.load(fout, encoding="latin1")
    prefix = 'm64182_220430_012143/'
    suffix = '/ccs'
    fhi = open(cluster_labels)
    for line in fhi:
        split = line.split()
        cluster = split[1]
        zmw = split[0].split('_')[1]
        key = prefix + zmw + suffix
        if key not in tipds: continue
        if len(tipds[key]) == 0: continue
        print('this happens')
        vecs = tipds[key][0]
        length = tipds[key][1]
        num_cpgs = len(vecs)
        cpg_content = num_cpgs / length * 1000
        vals = []
        for idx, val in vecs:
            vals.append(val)
        cpg_meth = np.mean(vals)
        print("%s\t%s\t%s\t%s" % (key, cluster, cpg_content, cpg_meth), file=fho)

def main():
    cluster_labels = sys.argv[1]
    cpg_pickle = sys.argv[2]
    fho = open(sys.argv[3], 'w')
    link_clusters_cpg(cluster_labels, cpg_pickle, fho)
    fho.close()

if __name__ == '__main__':
    main()
        
        
    
