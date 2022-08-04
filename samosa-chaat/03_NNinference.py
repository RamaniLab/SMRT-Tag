'''
03_NNinference.py
Colin McNally
2021/07/13

For each sample, load the formatted data, run through network inference, and produce results that are inputs for HMM binarization.
'''

import os
import pickle
from glob import glob
from tqdm import tqdm
import pandas as pd
import argparse
import keras ## version


def parse_args():
    parser = argparse.ArgumentParser(description="Run NN inference on formatted examples")
    parser.add_argument('samples', nargs='+',help='Either an integer index into the reference file, or a string identifier in the format [cell].[samplename]')
    parser.add_argument('-r','--referenceFile', nargs=1, required=True, help='The sample reference file from which to load sample information')
    parser.add_argument('-o', '--outputlocation',dest='output_directory',help='Location to save the outputs')
    parser.add_argument('-t', '--threshold', type=float, help='Threshold for defining IPD residuals as methylated. [default] 0.42',default=0.42)
    parser.add_argument('-m','--modeldir',nargs=1,dest='model_directory',help='Directory where IPD models are located')
    args = parser.parse_args()
    return args

def load_models(model_directory):
    ipdModel = keras.models.load_model('{}/ipdModel'.format(model_directory) )
    posResidModel = keras.models.load_model('{}/posResidualModel'.format(model_directory))
    posProbModel = keras.models.load_model('{}/posProbModel'.format(model_directory) )
    negProbModel = keras.models.load_model('{}/negProbModel'.format(model_directory) )

    return ipdModel,posResidModel,posProbModel,negProbModel
    
def makeHMMinput(samp,dataPBase,sampleRef,model_directory,threshold=0.42):

    ipdModel,posResidModel,posProbModel,negProbModel = load_models(model_directory)

    with open(dataPBase + '/processed/full/{}_{}_full_zmwinfo.pickle'.format(sampleRef['cell'][samp],sampleRef['sampleName'][samp]), 'rb') as fin:
        zmwinfo = pickle.load(fin, encoding="latin1")

    ## Load the formatted IPDs
    parts = sorted(glob(os.path.join(dataPBase,'processed','forNN',
                                    sampleRef['cell'][samp]+'_'+sampleRef['sampleName'][samp] + '*_forNN.npz')))
   
   
    pieceN = 0

    for inputPart in tqdm(parts, position=0, desc=sampleRef['sampleName'][samp]):
        
        ## takes a minute to load the data for 10K reads
        with np.load(inputPart) as data:
            ipdDat = data['ipdArr']
            contextDat = data['contextmat']
            percsDat = data['percsmat'][:,9:]
            zmwDat = data['zmwArr']
            zmwPosDat = data['zmwPos']

    
        predIPD = ipdModel.predict([contextDat, percsDat], batch_size=131072)
        resid = ipdDat - predIPD
        methPred = resid > threshold
        resPred = posResidModel.predict([contextDat], batch_size=131072)
        posProb = posProbModel.predict([contextDat], batch_size=131072)
        negProb = negProbModel.predict([contextDat], batch_size=131072)
        usePred = np.nonzero(resPred > 0.6)[0] 

        if not os.path.exists(dataPBase + '/processed/forHMM/'):
            os.makedirs(dataPBase + '/processed/forHMM/')


        hmmInput = {}
        ia = 0
        while ia < resPred.shape[0]:
            zmw = zmwDat[ia,0]
            istart = ia
            while ia < resPred.shape[0] and zmwDat[ia,0] == zmw:
                ia += 1
            iend = ia - 1

            zmwInd = np.arange(istart, iend + 1)
            goodInd = zmwInd[resPred[zmwInd,0] > 0.6]
            goodDat = pd.DataFrame({'pos':zmwPosDat[goodInd, 0], 'methPred':methPred[goodInd, 0], 'posProb':posProb[goodInd, 0],
                          'negProb':negProb[goodInd, 0]})
            goodDat.sort_values(axis=0, by='pos', inplace=True)
            goodDat.reset_index(inplace=True)

            if goodDat.shape[0] >= 2:
                hmmInput[zmw] = {'inDat':goodDat,
                                 'cclen':zmwinfo['cclen'][zmwinfo['zmw'] == zmw].iloc[0]}

            if len(list(hmmInput.keys())) >= 1000 or ia == resPred.shape[0]:
                
                with open(dataPBase + '/processed/forHMM/{0}_{1}_forHMM_resid-{2}_piece{3:05d}.pickle'.format(
                                                                                                   sampleRef['cell'][samp],
                                                                                                   sampleRef['sampleName'][samp],
                                                                                                   threshold,
                                                                                                   pieceN), 'wb') as fout:
                    pickle.dump(hmmInput, fout, protocol=4)
                pieceN += 1
                hmmInput = {}
                
def main():
    args=parse_args()

    ## read sampleRef
    sampleRef = pd.read_csv(args.referenceFile,sep=',',index_col='index')

    ## run NN inference
    for samp in args.samples:
        makeHMMinput(samp,args.output_directory,sampleRef,args.model_directory,threshold=args.threshold)

if __name__=='__main__':
    main()