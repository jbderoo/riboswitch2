# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 20:05:50 2020

@author: willi
"""

import pandas as pd
import numpy as np

import sys
import subprocess
import os
import time
import pickle

# mp stuff

from multiprocessing import Pool, TimeoutError, freeze_support
import multiprocessing

import itertools
from contextlib import contextmanager


def test_pool(args):
    return do_func(*args)


def do_func(f1, f2):
    return '{} & {}'.format(f1, f2)


def unpack_fnames(args):
    return align_readout(*args)


# def make_rs_fasta(name, seq, dot):
#   with open('rs_fixed.fasta','w') as f:
#     f.write(">" + name + "\n" +seq + "\n")
#     f.write(dot + " #FS" + "\n")

# def make_utr_fasta(name, seq):
#   with open('utr.fasta','w') as f:
#     f.write(">" + name + "\n" +seq + "\n")

# def make_non_fs_rs_fasta(name, seq):
#   with open('rs.fasta','w') as f:
#     f.write(">" + name + "\n" +seq + "\n")


# def make_utr_sub_fasta(name, seq):
#   with open('utr_sub.fasta','w') as f:
#     f.write(">" + name + "\n" +seq + "\n")

def align(file1, file2):
    aln = subprocess.check_output(('locarna ' + file1 + ' ' + file2 + ' --struct-local=1'), shell=True)
    aln_score = int(aln.decode('utf-8').split('\n')[0].split(' ')[1])
    aln_string = aln.decode('utf-8')
    return aln_score, aln_string


def align_readout(utr_file, rs_file):
    #e = errortext
    rs_score, a = align(rs_file, rs_file)
    aln_score, aln_report = align(utr_file, rs_file)

    return '{},{}|U:{},R:{}@{}'.format(utr_file, rs_file, str(aln_score), str(rs_score), aln_report)


def write_score_report(results):
    csv = []
    for i in range(len(results)):
        a = results[i]
        aln_report = a.split('@')[-1]
        a = a.split('@')[0]
        r = [x.split(',') for x in a.split('|')]
        utr_file = r[0][0]
        rs_file = r[0][1]

        utr_score = float(r[1][0].split(':')[1])
        rs_score = float(r[1][1].split(':')[1])
        per_score = utr_score / rs_score
        utr_file = utr_file.split('/')[-1]
        rs_file = rs_file.split('/')[-1]
        aln_report = aln_report.replace('\n','|')
        csv.append([utr_file, rs_file, utr_score, rs_score, per_score, aln_report])

    df = pd.DataFrame(csv, columns=['utr', 'rs', 'utr_score', 'rs_score', 'per_score', 'aln_report'])
    df.to_csv((save_dir + result_name + '.csv'))


@contextmanager
def poolcontext(*args, **kwargs):
    pool = multiprocessing.Pool(*args, **kwargs)
    yield pool
    pool.terminate()


if __name__ == '__main__':
    # start 4 worker processes
    utr_dir = sys.argv[1]
    rs_dir = sys.argv[2]
    save_dir = sys.argv[3]

    utr_files = os.listdir(utr_dir)
    rs_files = os.listdir(rs_dir)

    to_aln = ()
    for i in range(len(utr_files)):
        for j in range(len(rs_files)):
            to_aln = (((utr_dir+utr_files[i]), (rs_dir+rs_files[j]),),) + to_aln

    utr_names = [x.split('.')[0].split('_')[-1] for x in utr_files]

    result_name = 'test1'
    errortext = './Err/err.txt'

    print('Aligning %i Combinations' % len(to_aln))
    print('aligning {} to {}'.format(utr_dir, rs_dir))
    print('starting pool....')
    st = time.time()
    with poolcontext() as pool:

        r = pool.starmap(align_readout, to_aln)
    #r = align_readout( to_aln[0][0] , to_aln[0][1] )
    results = r
    #print(results)

    print('pool time:')
    print(time.time() - st)
    print('writing report...')

    write_score_report(results)



