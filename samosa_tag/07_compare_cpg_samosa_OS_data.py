'''
compare_cpg_samosa_OS_data.py
Vijay Ramani

This script links CpG information and m6dA signal between ZMWs.
'''
import os,sys,re
import pickle
import numpy as npd

def parse_args():
    parser = argparse.ArgumentParser(description='Links CpG information and m6dA signal between ZMW')
    parser.add_argument('-c','--cpg-pickle',nargs=1,dest='cpg_pickle',help='CpG data per read (produced by 04_cpg2pickle.py)',required=True)
    parser.add_argument('-l','--cluster-labels',nargs=1,dest='cluster_labels',help='Cluster assignments per read (produced by 02_cluster_autocorrelograms_tn5.py)',required=True)
    parser.add_argument('-o','--output-dir',nargs=1,dest='output_directory',help='Output directory',required=True)
    parser.add_argument('-p','--project',nargs=1,dest='project',help='Name of the project',required=True)

    args = parser.parse_args()
    return args



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
        key = prefix + zmw + suffix ## specify ZMWids 
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
    args=parse_args()
    #cluster_labels = sys.argv[1] # autocor_clusters_Tn5_2.data
    #cpg_pickle = sys.argv[2] # OS152_plusM_cpg_data.pickle
    cluster_labels=args.cluster_labels
    cpg_pickle=args.cpg_pickle

    fho = open('{}/{}_cluster_cpg.dataout.filtered'.format(args.output_directory,args.project), 'w') ##
    
    link_clusters_cpg(cluster_labels, cpg_pickle, fho)
    fho.close()

if __name__ == '__main__':
    main()
        
        
    
