# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 20:07:08 2018

@author: Wajid Abbasi
"""

# -*- coding: utf-8 -*-
import random,os
from Bio.Data import IUPACData
from Bio.SubsMat.MatrixInfo import blosum62
import numpy as np
from propy.PyPro import GetProDes
#from sklearn.externals import joblib
#from sklearn import svm, cross_validation,ensemble
window_size = 21
half_window_size = window_size / 2

def plot_binding_site(BS_results,seq_id):
    if not os.path.exists('binding_site_plots'):
        os.makedirs('binding_site_plots')
    seq_id="".join(x for x in seq_id if x.isalnum())
    import matplotlib.pyplot as plt
    import numpy as np
    BS_plot='binding_site_plots/'+seq_id+'_BS_plot.jpg'
    max_score,max_position,scores=BS_results[0],BS_results[1],BS_results[2]
    fig=plt.figure()
    plt.ioff()# prevent plot to display
    plt.plot(range(1,scores.shape[0]+1),scores,'b-')
    plt.plot(max_position,max_score,'ro')
    plt.xlabel('Position')
    plt.ylabel('Score')
    plt.title('Max score is %f and max position is %d' %(max_score,max_position))
    plt.grid()
    fig.savefig(BS_plot)
    
def return_blosum62_dict():
    prot_dic = dict((k, 0) for k in IUPACData.protein_letters)
    A=[]
    Dictblosum={}
    for aa in prot_dic:
        for bb in prot_dic:
            if (aa,bb) in blosum62.keys():
                A.append(blosum62[(aa,bb)])
            else:
                A.append(blosum62[(bb,aa)])
        Dictblosum[aa]=A
        A=[]
    Dictblosum['X']=[0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2,0,0,-2,-1,-1]
    Dictblosum['B']=[-2,-1,3,4,-3,0,1,-1,0,-3,-4,0,-3,-3,-2,0,-1,-4,-3,-3]
    Dictblosum['Z']=[-1,0,0,1,-3,3,4,-2,0,-3,-3,1,-1,-3,-1,0,-1,-3,-2,-2]
    return Dictblosum
def blosum62_protein_features(seq):#Simply Compute Amino Acid Composition
    Dictblosum=return_blosum62_dict()
    features_list=[]
    feature=[]
    AA=["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I","L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    for i in range((len(seq)-window_size)+1):
        wenseq= str(seq[i:window_size+i])
        wenseq=wenseq.replace('U',AA[random.randint(0,len(AA)-1)])
        for j in range(len(wenseq)):
            feature.extend(Dictblosum[wenseq[j]])
        #feature=np.mean(feature,axis=0)
        feature=feature/np.linalg.norm(feature)
        features_list.append(feature)
        feature=[]
    return features_list
def mean_varrianace_normalization_single(examples):
    arab_propy_feats_mean=np.load(os.path.join('','arab_propy_feats_mean.npy'))
    arab_propy_feats_std=np.load(os.path.join('','arab_propy_feats_std.npy'))
    examples=(examples-arab_propy_feats_mean)/(arab_propy_feats_std)

    examples=(examples/np.linalg.norm(examples))

    return examples
def get_features_proPy_list_Full_des(seq):
        """
        This function takes fasta record and save path as input and generate features using
        proPy packege. IN this a search is implemented to get dictionary of all
        terms to get same feature lenght. Normalized fatures will
        be saved at given path with prot id.
        """
        seq=str(seq)
        AA=["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I","L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
        keyList=np.load(os.path.join('','key_list_proPy.npy'))
        feat=[]
        seq=seq.replace('X',AA[random.randint(0,len(AA)-1)])
        seq=seq.replace('B',AA[random.randint(0,len(AA)-1)])
        seq=seq.replace('U',AA[random.randint(0,len(AA)-1)])
        seq=seq.replace('Z',AA[random.randint(0,len(AA)-1)])
        des=GetProDes(seq)
        propy_FD=des.GetALL()
        for key in keyList:
            try:
                feat.append(propy_FD[key])
            except KeyError:
                feat.append(0.0)
        feat=mean_varrianace_normalization_single(feat)
        return feat
def validate(sequence):
    amino_acids = 'ACDEFGHIKLMNPQRSTVWYXBUZ'
    sequence = sequence.upper()
    #if len(set(sequence).difference(amino_acids)):
       # raise ValueError("Input sequence contains non-standard amino acid codes")
    if len(sequence)<window_size or len(set(sequence).difference(amino_acids))>0:
        #raise ValueError("Input sequence is shorter than %d amino acids." %window_size)
        return 0
    else:
        return 1
def window_scores(features):
    """
    Obtain raw window scores for all windows
    """
    filepath = os.path.join('', 'blosum_weight_vector.npy')
    W_Blosum=np.load(filepath)
    scores =  np.dot(features,W_Blosum)
    max_position = np.argmax(scores)
    max_score = np.round(scores[max_position],decimals=2)
    scores=np.round([np.max(scores[max(0,i-20):min(i+1,len(scores))]) for i in range(len(scores))],decimals=2)
    print ('Position of central residue of maximum scoring window:', max_position+11)
    print ('Maximum Score:', max_score)
    print ('Scores for  central residue of all windows',scores)
    return (max_score,max_position+11,scores)
def predict_BS(sequence):
    if not(validate(sequence)):
        print('Input sequence is shorter than 21 amino acids')
        return
    feats=blosum62_protein_features(sequence)
    Scores=window_scores(feats)
    return Scores
def Gaussian(x,z,gamma):
    return np.exp(-(gamma*(np.linalg.norm(x-z)**2)))
def Gaussian_matrix(X,Y):
    G_matrix=[]
    for i in range(len(X)):
        G_matrix.append(Gaussian(X[i],Y,2.0))
    return np.asarray(G_matrix)
def CaM_interaction(feats):
    intercept=np.load('clf_intercept.npy')
    coeff=np.load('clf_coefficients.npy')
    S_vectors=np.load('clf_S_vectors.npy')
    Guassian=Gaussian_matrix(S_vectors,feats)
    score=np.dot(coeff,Guassian)+intercept
    A = -9.16386985458
    B = -4.3849677898
    score = 100.0/(1+np.exp(score*A+B))
    return round(score,2)
def predict_CaM_interaction(sequence):
    if not(validate(sequence)):
        print('Input sequence is shorter than 21 amino acids')
        return
    feats=get_features_proPy_list_Full_des(sequence)
    feats=feats.reshape(1, -1)
    score=CaM_interaction(feats)
    return score
