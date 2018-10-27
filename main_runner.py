# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 20:30:04 2018

@author: Wajid Abbasi
"""
import os
from Bio import SeqIO
from CaMELS import *
import base64

def plot_binding_site(BS_results,seq_id):
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
    
    
all_seqs_recs=list(SeqIO.parse("example.fasta", "fasta"))#Give your fasta file instead of example.fasta

for rec in all_seqs_recs:
    print (rec.id)
    if validate(rec.seq)==0:
        print('Input sequence is not valid')
    print 'Interaction Prediction Score:',predict_CaM_interaction(rec.seq)
    print 'Binding Site Prediction:'
    binding_site_results=predict_BS(rec.seq)
    plot_binding_site(binding_site_results, rec.id)
