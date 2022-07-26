import os,sys,re
import scipy
import numpy as np
import pysam 
from collections import Counter
#Author: VR
#Goal of this script is to take an aligned, sorted PacBio BAM file, as well as a list of genomic sites formatted
#as chrid\tcoordinate\tstrand and tabulate all read ends falling in a 5kb window centered on the site. The script
#returns a numpy array where rows are individual sites and counts at each distance from center are columns.

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
  valid = {}
  tfile = open(sys.argv[1])
  valid_chroms = open(sys.argv[2])
  for line in valid_chroms:
    split = line.split()
    valid[split[0]] = True
  tss = tfile.readlines()
  filename = sys.argv[3]
  label = sys.argv[4]
  sim_type = sys.argv[5]
  if sim_type == "True":
    sim_type = True
  else:
    sim_type = False
# sim_type = True
  filenames = [filename]
  vplots = calculatePerBase(filenames, tss, "long", sim_type, valid)
  np.save("TSS_mat.ends.%s" % label, vplots)
  tfile.close()
  valid_chroms.close()

if __name__ == "__main__":
  main()
