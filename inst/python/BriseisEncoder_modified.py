##########################################################################################################################
#BriseisEncoder 1.0.0: A deep learning encoder for CDR3 sequences in the TCR space
##########################################################################################################################
#BriseisEncoder is capable of transforming CDR3 peptide sequences from productive
#TCR-beta chains into 15-digit informative numerical vectors.
#01312019 Ze Zhang <Ze.Zhang@UTsouthwestern.edu>
##########################################################################################################################
#Dependencies: 
#Python 3.6.4 (preferred)
#numpy (>=1.15.4), pandas (>=0.23.4), keras (>=2.2.4)
##########################################################################################################################
#Parameters:
#tcr_dir: a .csv file containing columns 'contig_id' and 'cdr3'.
#   'contig_id' contains non-repeated ID strings
#   'cdr3' contains valid TRB CDR3 peptide sequences.

#aa_dict_dir: a .csv file storing Atchley's amino acid embedding, details in https://www.ncbi.nlm.nih.gov/pubmed/15851683.

#model_dir: a .h5 file containing trained encoding models from CDR3 seqs called from bulk-tumor RNA-seq data.
##########################################################################################################################

#Import dependencies
import numpy as np
import pandas as pd
import os
import csv
from keras.models import load_model
from keras.models import Model

# Main function
#def BriseisEncoder(tcr_dir, aa_dict_dir, model_dir):
tcr_dir = 'Atchley_factor_tcr_only_tmp.csv'
sysfiles = pd.read_csv('temp_sysfiles.csv')['sysfiles']
aa_dict_dir = sysfiles[0]
model_dir = sysfiles[1]
encode_dim=80
#Define helper functions
def preprocess(filedir):
    #Preprocess TCR files
    print('Processing: '+filedir)
    if not os.path.exists(filedir):
        print('Invalid file path: ' + filedir)
        return 0
    dataset = pd.read_csv(filedir, header=0)
    if dataset.isnull().values.any():
        print('Input data contains NAs.')
        dataset=dataset.dropna()
    data_new=pd.DataFrame({
        'contig_id':dataset['contig_id'],
        'cdr3':dataset['cdr3']
    })
    data_new.index=range(0,dataset.shape[0])
    return data_new

def aamapping(peptideSeq, aa_dict, encode_dim):
    #Transform aa seqs to Atchley's factors.
    peptideArray = []
    if len(peptideSeq)>encode_dim:
        print('Length: '+str(len(peptideSeq))+' over bound!')
        peptideSeq=peptideSeq[0:encode_dim]
    for aa_single in peptideSeq:
        try:
            peptideArray.append(aa_dict[aa_single])
        except KeyError:
            print('Not proper aaSeqs: '+peptideSeq)
            peptideArray.append(np.zeros(5,dtype='float64'))
    for i in range(0,encode_dim-len(peptideSeq)):
        peptideArray.append(np.zeros(5,dtype='float64'))
    return np.asarray(peptideArray)

def datasetMap(dataset, aa_dict, encode_dim):
    #Wrapper of aamapping
    TCR_dict=dict()
    for i in range(0,len(dataset['cdr3'])):
        TCR_key=dataset['contig_id'][i]
        if TCR_key in TCR_dict.keys():
            TCR_key=TCR_key+'_'+str(i)
        TCR_dictarray=aamapping(dataset['cdr3'][i],aa_dict,encode_dim)
        TCR_dict[TCR_key]=TCR_dictarray
    return TCR_dict

# Perform encoding
print('BriseisEncoder: Mission loading.')
tcr=preprocess(tcr_dir)
aa_dict=dict()
with open(aa_dict_dir,'r') as aa:
    aa_reader=csv.reader(aa)
    next(aa_reader, None)
    for rows in aa_reader:
        aa_name=rows[0]
        aa_factor=rows[1:len(rows)]
        aa_dict[aa_name]=np.asarray(aa_factor,dtype='float')
TCR_dict=datasetMap(tcr, aa_dict, encode_dim)
TCR_contigs=np.stack(TCR_dict.values())
TCR_contigs=TCR_contigs.reshape(-1,encode_dim,5,1)
#Model prediction
TCRencoder=load_model(model_dir)
encoder=Model(TCRencoder.input,TCRencoder.layers[-12].output)
encoded_mat=encoder.predict(TCR_contigs)
encoded_mat=pd.DataFrame(encoded_mat,index=tcr['contig_id'])
encoded_mat.to_csv('temp_atchley_factors_encoded.csv', sep = ',')
print('BriseisEncoder: Mission Accomplished.\n')
    
#    return encoded_mat
