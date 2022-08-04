'''
06_annotate_bams.py
Scott Nanda
2022/08/03

Annotate individual reads with MM and ML tags defined by accessibility profiles predicted by the SAMOSA-ChAAT pipeline.
'''
import sys
import os
import numpy as np ## v1.23.0
import pandas as pd
import pickle
from tqdm import tqdm
from numba import njit ## v0.53.1
from numba.typed import Dict
from numba.core import types
import glob
from collections import defaultdict
import pysam


def parse_args():
    parser = argparse.ArgumentParser(description="Annotate SAMOSA reads with accessibility footprints")
    parser.add_argument('samples', nargs='+',help='An integer index into the reference file')
    parser.add_argument('-m','--merged-bam', dest='merged_bam',nargs=1, required=True, help='The sample reference file from which to load sample information')
    parser.add_argument('-r','--referenceFile', nargs=1, required=True, help='The sample reference file from which to load sample information')
    parser.add_argument('-o', '--outputlocation',dest='output_directory',help='Location to save the outputs')

    args = parser.parse_args()
    return args


## Define integer dict for type conversion
int_dict=Dict.empty(
    key_type=types.unicode_type,
    value_type=types.int64
)
for key in range(50000):
    int_dict[str(key)]=key


def compress_state_vector(state_vectors):
    '''
    Function to compress state vectors into state strings, run length encoded representations of binarized accessibility calls
    Requires as input a dictionary with ZMW IDs as keys, and state vectors as values, as is produced by HMMv2, or
    multiple HMMv1 calls. 
        
    Parameters:
        state_vectors : dict
            Dictionary with ZMW IDs as keys and binarized acccesibility calls in a np.array as values. 
        
    Returns:
        state_strings : dict
            Dictionary with ZMW IDs as keys and state strings, run length encoded repesentations of binarized accessibilty calls,
            as values. 
    '''
    
    state_strings=dict()
    
    for k,v in tqdm(state_vectors.items()):
        state_strings[k] = _compress_state_vector_to_state_string(v)
        
    return state_strings


@njit
def _compress_state_vector_to_state_string(state_vector):
    '''
    Function to compress a state vector (binarized accessibility numpy array) into a "state-string", a
    compact run-length encoded representation of accessibility calls for a single molecule
    
    Parameters:
        state_vector : np.array
            Binarized accessibiltiy stored as numpy array with 1/0, and continous values at changepoints. np.nan is used to indicate values that could not be inferred (molecule ends).
    
    Returns:
        compressed_string : str
            A run-length encoded representation of the state vector, using "N,C,I,A" to represent "NaN", "changepoints",
            "accessible", "inaccessible" calls. Accessible and Inacessible calls are (1,0) respectively, and Changepoint 
            calls represent the region of uncertainty between accessible and inaccessible calls that is represeted by
            linear interpolation between the values i.e {0.2, 0.4, 0.6...et}.
    
    '''
    
    compressed_string=""

    ## Start at state 0
    state = state_vector[0]

    ## Depending on the state, pick a starting character. 
    if state==0:
        prev_char="I"
    elif state==1:
        ## 1 is accessible
        prev_char="A"
    elif np.isnan(state):
        ## N is NaN state 
        prev_char="N"
    else:
        ## if value is float, it is part of a changepoint
        prev_char="C"

    ## Increase the processed count
    count=1

    
    ## Set up a variable to store the state character
    current_char=""
    ## For each state in the state_vector
    for state in state_vector[1:]:
        ## simple switch statement 
        if state==0:
            ## 0 is inaccessible
            current_char="I"
        elif state==1:
            ## 1 is accessible
            current_char="A"
        elif np.isnan(state):
            ## N is NaN state 
            current_char="N"
        else:
            ## if value is float, it is part of a changepoint
            current_char="C"

        ## if switching characters, save the result 
        if current_char!=prev_char:
            compressed_string+=str(count)+prev_char
            ## reset count
            count=1
            prev_char=current_char

        else:
            ## otherwise increment the count
            count+=1
    
    ## one final round of addition for the final character
    if current_char==prev_char:
        compressed_string+=str(count)+prev_char

    return compressed_string



@njit
def convert_state_string_to_MM_ML(state_string,int_dict):
    '''
    Function to convert state strings into associated total-molecule MM and ML tags i.e using each base
    
    Parameters:
        state_string : str
            A run-length encoded representation of a state vector. 

        int_dict : numba.Dict
            Dictionary containing int values as keys, and the character representation of those integer values as values. 
            Allows for conversion between int and char, which is not supported in numba natively. 
    
    Returns:
        MM : numba.List
            A list containing modification positional encodings as specified in the SAMv4.2 format, where every base
            is considered as being potentially modified
        
        ML : numba.List
            A list containing modification likelihoods as specified in the SAMv4.2 format (uint8), where every base
            is considered as being potentially modified. 
    
    '''
    
    
    ## parse the string on the fly and assemble the vector instead of
    ## using a pre-allocated vector or two passes.
    cur_count=""
    cur_count_int=0
    ops=["N","C","I","A"]

    append_dict={"I":0,"A":1}

    MM=[]
    ML=[]
    previous_skip=0
    prev_op=""

    ## first loop
    # for chr in state_string:
    for i in range(len(state_string)):
        char=state_string[i]
        
        ## for each character in the state string
        if char in ops:
            cur_count_int=int_dict[cur_count]#int(cur_count)

            if char=="N":
                ## increment the previous_skip, but dont add it 
                previous_skip+=cur_count_int
            else:
                ## always happens 
                MM.append(previous_skip)

                if char=="C":                
                    append_prob = append_dict[prev_op] ## 0 / 1 
                    numerator=1-2*append_prob
                    slope=numerator/(cur_count_int+1)
                    
                    ML.append(int(255*(slope+append_prob)))

                    for j in range(1,cur_count_int):
                        MM.append(0)
                        ML.append(int(255*(slope*(j+1)+append_prob)))

                else: ## I/A
                    append_prob = int(255*(append_dict[char])) ## 0 / 1 
                    ML.append(append_prob)

                    for j in range(1,cur_count_int):
                        MM.append(0)
                        ML.append(append_prob)
    
                previous_skip=0
                
            ## store the op as previous op
            prev_op=char
            ## reset cur_count
            cur_count=""

        else:
            cur_count+=char
        
    return MM, ML


def main():
    args=parse_args()

    ## read sampleRef
    sampleRef = pd.read_csv(args.referenceFile,sep=',',index_col='index')


    ## 1) load HMM results
    hmmDict={}
    for sample in samples:
        with open('{}/processed/binarized/HMMout/{}_{}_NNsingle.pickle'.format(args.output_directory,sampleRef['cell'][samp],sampleRef['sampleName'][samp]),'rb') as f:
            hmmDict.update(pickle.load(f))

    state_strings = compress_state_vector(hmmDict)

    ## 2) Load the ccs file associateed with the reads
    infile = pysam.AlignmentFile(args.merged_bam,check_sq=False)
    ## Set up an output file with the same header
    outfile = pysam.AlignmentFile('{}/processed/annot/{}'.format(args.output_directory,os.path.basename(args.merged_bam).replace('.bam','.samosa.bam')),'wb',check_sq=False,template=infile,threads=5)

    ## 3) Write out reads + MM and ML information
    ## for each read in the bamfile
    for read in tqdm(infile):
        ## get the ZMW
        zmw = read.get_tag('zmw')
        
        ## if the ZMW is in the set of molecules solved by HMM, then process
        ## otherwise skip
        if zmw in state_strings:
            ## generate MM and ML tags
            mm,ml= convert_state_string_to_MM_ML(state_strings[zmw],int_dict)
            
            ## Overwrite the tags present and add the MM and ML tags
            ## Use pacbio style of "Mm" vs "MM" and "Ml" vs "ML"
            ## use the "N+a" spec to indicate these are N=any base, and a=m6dA (not strictly possible)
            ## but for our purposes, records that the modification is NOT anything other than m6dA
            read.set_tags(
                read.get_tags()+
                [
                ("Mm","N+a,{};".format(",".join(map(str,mm))),"Z"),
                ("Ml",ml,"C")
                ]
            )

            
            outfile.write(read)

    infile.close()
    outfile.close()


if __name__=='__main__':
    main()