# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 20:05:50 2020

@author: willi
"""

import sys
import subprocess
import os
import time

from os import path
from csv import reader

utr_dir         = '../data/utr_negative_aln/lstm100/'

rs_dir          = '../data/utr_negative_aln/eukaryotic_RS_fastas/'
save_dir        = '../data/utr_negative_aln/utr_aln_combined/'

rs_self_file    = rs_dir + 'self_compare_list.txt'
rs_list         = {}

def align(file1, file2, e):
    aln = subprocess.check_output((f'locarna --struct-local=1 {file1} {file2}'), shell=True, stderr=e)
    aln_score = int(aln.decode('utf-8').split('\n')[0].split(' ')[1])
    return aln_score


def load_rs_values(rs_file):
    with open(rd_dir + rs_self_file, 'r') as fh:
        csv_reader = reader(fh)
        for line in csv_reader:
            rs_list[line[0]] = line[1]
    return rs_list


if os.is_file(rd_dir + rs_self_file):
    rs_list = load_rs_values(rd_dir + rs_self_file)
else: 
    rs_files = os.listdir(rs_dir)
    e = open('err.txt', 'w')
    with open(rd_dir + rs_self_file, 'w') as fh:
    for rs_file in rs_files:
        token = rs_file.split('_')[1]
        ref_score = align(rs_dir + rs_file,rs_dir + rs_file, e) # can we move this to outside the j loop and into the i loop,
        fh.write(f'{token},{ref_score}')
    rs_list = load_rs_values(rd_dir + rs_self_file)


utr_names = [x.split('.')[0].split('_')[-1] for x in utr_files]

for i in range(len(utr_files)):
    print('i is', i)
    if i >= 1:
        break

    print('There are', len(utr_files), 'utr files (parent loop)')
    print('There are', len(rs_files), 'rs files (child 1 loop)')

    utr = utr_files[i]
    utr_name = utr_names[i]
    e = open('err.txt', 'a')

    print('-------------------')
    print('aligning: %s' % utr)
    print('starting at: %s' % time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))

    st = time.time()

    alignment_scores = np.zeros((len(rs_list), 3))
    alignment_names = []
    start = time.time()

    for j in range(len(rs_list)):

        if j % 5 == 0:
            diff = time.time() - start
            print('I have completed', j, 'alns, it has taken', diff, 'to do alns', (j - 5), 'thru', j)
            start = time.time()

        rs = rs_list[j]

        # ref_score = align(rs_dir + rs,rs_dir+rs,e) # can we move this to outside the j loop and into the i loop,
        # or completely separate program that calculates and then fetch in loop

        utr_score = align(utr_dir + utr, rs_dir + rs, e)

        alignment_scores[i, :] = [ref_score, utr_score, float(utr_score) / ref_score]
        alignment_names.append((utr, rs))
        # print('The score is:',float(utr_score)/ref_score)

    e.close()
    print('total time: %s' % time.time() - st)
    print('finished at: %s' % time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))

    np.save((save_dir + utr_name), alignment_scores)

    with open((utr_name + '_names.p'), 'wb') as handle:
        pickle.dump(alignment_names, handle)

