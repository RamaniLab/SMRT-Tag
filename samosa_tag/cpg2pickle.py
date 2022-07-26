import pysam
import re
import pickle

#Convert run-length-encoded primrose output into a pickle file (a la m6dA storage)

#methylation analysis
cnt = 0
run_ids = {}

samfile= pysam.AlignmentFile("/avicenna/vramani/data/pacbio/tn5/B02/demux/Tn5.A13.subread.OS152_plusM.ccs_kinetics.primrose.align.sorted.bam", "rb", check_sq=False)

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

with open('Tn5_OS_total_data.pickle', 'wb') as handle:
    pickle.dump(run_ids, handle, protocol=pickle.HIGHEST_PROTOCOL)
