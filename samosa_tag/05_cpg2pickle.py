'''
05_cpg2pickle.py
Vijay Ramani

This script converts run-length-encoded primrose outputs in BAM files into continuous value vectors (per read), and aggregates them into a pickle file.
'''
import pysam
import re
import pickle
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Converts run-length-encoded primrose outputs in BAM files into a pickle file')
    parser.add_argument('primrose_bams', nargs='+', help='Primrose-annotated BAM files')
    parser.add_argument('-o','--output-dir',nargs=1,dest='output_directory',help='Output directory',required=True)
    parser.add_argument('-p','--project',nargs=1,dest='project',help='Name of the project',required=True)

    args = parser.parse_args()
    return args


def main():
    cnt = 0
    run_ids = {}

    args=parse_args()
    
    for _samfile in args.primrose_bams:
        samfile=pysam.AlignmentFile(_samfile,'rb',check_sq=False)
        for record in samfile:
            title = record.query_name
            if title not in run_ids:
                run_ids[title] = []
            seq = record.get_forward_sequence()
            c_loc = [i.start() for i in re.finditer('C',seq)]
            try:
                c_ind = record.get_tag('MM').split(',')
            except KeyError:
                continue
            if len(c_ind) <= 2: continue
            start = c_loc[int(c_ind[1])]
            adds = []
            meth_vals = []
            for char in c_ind[2:]:
                adds.append(int(char.split(';')[0]))
            c_meth = record.get_tag('ML')
            for val in c_meth[1:]:
                meth_vals.append(val / 255)
            c_meth_start = c_meth[0] / 255
            meth_vec = []
            meth_vec.append((start, c_meth_start))
            zipped = zip(adds, meth_vals)
            adder = 0
            for add, meth_val in zipped:
                adder += add + 1
                idx = c_loc[int(c_ind[1]) + adder]
                meth_vec.append((idx,meth_val))
            run_ids[title] = (meth_vec, len(seq))
        
        ## Close bam file
        samfile.close()

    ## Save aggregated CpG data
    with open('{}/{}_cpg_data.pickle'.format(args.output_dir,args.project), 'wb') as handle:
        pickle.dump(run_ids, handle, protocol=pickle.HIGHEST_PROTOCOL)

if __name__=='__main__':
    main()
