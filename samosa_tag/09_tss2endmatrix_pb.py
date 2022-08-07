'''
tss2endmatrix_pb.py
Vijay Ramani

This script takes an aligned, sorted PacBio BAM file, as well as a list of genomic sites formatted
as chrid\tcoordinate\tstrand and tabulates all read ends falling in a 5kb window centered on the site. The script
returns a numpy array where rows are individual sites and counts at each distance from center are columns.
'''

import os,sys,re
import scipy
import numpy as np
import pysam 
from collections import Counter

def parse_args():
    parser = argparse.ArgumentParser(description='Takes an aligned, sorted PacBio BAM file, as well as a list of genomic sites formatted as chrid\tcoordinate\tstrand and tabulates all read ends falling in a 5kb window centered on the site.')
    parser.add_argument('aligned_bam_files', nargs='+', help='Aligned PacBio BAM files.')
    parser.add_argument('-g','--genomic-sites',nargs=1,dest='genomic_sites',help='TSV file containing genomic sites to be examined for read end enrichment',required=True)
    parser.add_argument('-n','--genomic-sites-name',nargs=1,dest='genomic_sites_name',help='Label / Name to describe set of sites examined for enrichment',required=True)
    parser.add_argument('-c','--valid-chroms',nargs='1',dest='valid_chroms',help='Comma separated string of chromsomes to examine for read end enrichment',required=True)
    parser.add_argument('-m','--simulation-mode',nargs=1,dest='simulation_mode',help='Output directory',required=True,)

    parser.add_argument('-o','--output-dir',nargs=1,dest='output_directory',help='Enable simulation mode',required=True,action='store_false',default=True)
    parser.add_argument('-p','--project',nargs=1,dest='project',help='Project',required=True)

    args = parser.parse_args()
    return args

def readIterator(filenames,chrom, start, end):
  for bamfile in filenames:
    if os.path.exists(bamfile) and (os.path.exists(bamfile.replace(".bam",".bai")) or os.path.exists(bamfile+".bai")):
      input_file = pysam.Samfile( bamfile, "rb" )
      for read in input_file.fetch(chrom, start, end):
        yield read
      input_file.close()

def calculatePerBase(filenames, tss, mode, sim_type, valid_chroms):
  mat_total = []

  for line in tss:
    #print("this happens")
    split = line.split()
    chrom = split[0]
    if chrom not in valid_chroms: continue
    t_start, strand = int(split[1]), split[2]
    start = t_start - 2500
    if start <= 0: continue
    end = t_start + 2500
    end_mat = np.zeros(5001)
    for read in readIterator(filenames, chrom, start, end):
        rstart = read.reference_start
        rend = read.reference_end
        rlength = rend - rstart 
        rmid = (int(rend) + int(rstart)) / 2
        if strand == "-":
            register1 = t_start - rend + 2500
            register2 = t_start - rstart + 2500
        elif strand == '+':
            register1 = rstart - t_start + 2500
            register2 = rend - t_start + 2500
        if register1 < 5001 and register1 >= 0:
            end_mat[register1] += 1
        if register2 < 5001 and register2 >= 0:
            end_mat[register2] += 1
    mat_total.append(end_mat)
  mat_total = np.vstack(mat_total)
  return mat_total

def main():
  args=parse_args()

  with open(args.genomic_sites) as fho:
    tss=fho.readlines()

  valid = {}
  valid_chroms=args.valid_chroms.split(',')
  for v in valid_chroms:
    valid[v]=True

  sim_type=args.simulation_mode

  for filename in args.aligned_bam_files:
      filenames = [filename]
      vplots = calculatePerBase(filenames, tss, "long", sim_type, valid)
      np.save("{}/{}_mat.ends.{}".format(args.output_directory,args.genomic_sites_name,os.path.basename(filename)), vplots)

if __name__ == "__main__":
  main()
