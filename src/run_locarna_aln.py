# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 20:05:50 2020

@author: willi
"""

import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import sys
import subprocess
import os
import time
import pickle

targets = pd.read_csv("./cnn_hits.csv") 

utr_dir = './combined/'
rs_dir = './rs_fa/'
rs_dir_fs = './rs_fa_fs/'
save_dir = './utr_aln/'

def make_rs_fasta(name, seq, dot):
  with open('rs_fixed.fasta','w') as f:
    f.write(">" + name + "\n" +seq + "\n")
    f.write(dot + " #FS" + "\n")  

def make_utr_fasta(name, seq):
  with open('utr.fasta','w') as f:
    f.write(">" + name + "\n" +seq + "\n")

def make_non_fs_rs_fasta(name, seq):
  with open('rs.fasta','w') as f:
    f.write(">" + name + "\n" +seq + "\n")
    

def make_utr_sub_fasta(name, seq):
  with open('utr_sub.fasta','w') as f:
    f.write(">" + name + "\n" +seq + "\n")   
   
   
def align(file1,file2, e):
    aln = subprocess.check_output( ('locarna ' + file1 + ' ' + file2  +' --struct-local=True'), shell=True,stderr= e)
    aln_score = int(aln.decode('utf-8').split('\n')[0].split(' ')[1])
    return aln_score





    
utr_files = os.listdir(utr_dir)
rs_files = os.listdir(rs_dir)
rs_fs_files = os.listdir(rs_dir_fs)

utr_names = [x.split('.')[0].split('_')[-1] for x in utr_files]

for i in range(len(utr_files)):
    utr = utr_files[i]
    utr_name = utr_names[i]
    e = open('err.txt','w')
    
    
    print('-------------------')
    print('aligning: %s'%utr)
    print('starting at: %s'% time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))
    
    
    alignment_scores = np.zeros((len(rs_files),3)) 
    alignment_names = []
    alignment_scores_fs = np.zeros((len(rs_files),3)) 
    alignment_names_fs = []    
    for j in range( len(rs_files)):
        rs = rs_files[j]
        ref_score = align(rs,rs)
        utr_score = align(utr,rs )
        
        alignment_scores[i,:] =  [ref_score,  utr_score,   float(utr_score)/ref_score]
        alignment_names.append( (utr,rs)  )
    

    for j in range( len(rs_files)):
        rs = rs_files[j]
        ref_score = align(rs,rs)
        utr_score = align(utr,rs )
        alignment_scores_fs[i,:] =  [ref_score,  utr_score,   float(utr_score)/ref_score]
        alignment_names_fs.append( (utr,rs)  )        
        
    e.close()
    
    print('finished at: %s'% time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))
    
    np.save((save_dir +  utr_name ), alignment_scores)
    np.save((save_dir +  utr_name + '_fs' ), alignment_scores_fs)
    
    with open((utr_name + '_names.p'), 'wb') as handle:
        pickle.dump(alignment_names, handle)
    with open((utr_name + '_names_fs.p'), 'wb') as handle:
        pickle.dump(alignment_names_fs, handle)        
