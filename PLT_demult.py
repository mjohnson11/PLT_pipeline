"""
A program to demultiplex reads for the PLT
Milo Johnson
Modifiying from PLT_parse in March 2022
"""

import time
import gzip
import pandas as pd
import re
import csv
import numpy as np
import os
import subprocess
from PLT_tools import FourLineFastq

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('primer_map_file', help='file containing primer demultiplex info')
parser.add_argument('out_base', help='basic output directory')
parser.add_argument('R1_files', help='comma separated list of R1 files to parse - '
                                     '"R1" will be changed to "R2" to get the paired read files')
parser.add_argument('run_name', help='name for the run')
args = parser.parse_args()

# POSITIONAL ARGS - THESE SHOULD ALWAYS BE THE SAME FOR A RUN
primer_map_file = args.primer_map_file
R1_files = args.R1_files.split(',')
R2_files = [i.replace('R1', 'R2') for i in R1_files]
out_base = args.out_base
run_name = args.run_name
output_dir = out_base + run_name + '/'
reads_out_base = output_dir + 'read_files/'

for odir in [output_dir, reads_out_base]:
    if not os.path.isdir(odir):
        print('Making directory:', odir)
        subprocess.call(['mkdir', odir])



otime = time.time()

dem_df = pd.read_csv(primer_map_file)
R1_index_lens = list(set([len(ind) for ind in dem_df['R1_index']]))
R2_index_lens = list(set([len(ind) for ind in dem_df['R2_index']]))
if len(R1_index_lens) != 1:
    raise Exception('R1 indices are not all the same length')
if len(R2_index_lens) != 1:
    raise Exception('R2 indices are not all the same length')
R1_ind_len = R1_index_lens[0]
R2_ind_len = R2_index_lens[0]

libraries = list(dem_df['Library'])
num_libraries = len(libraries)
index_list = np.array(dem_df[['R1_index', 'R2_index']])
index_dict = {index_list[i][0]+index_list[i][1]: i for i in range(len(index_list))}
R1_bp_to_bc = list(dem_df['R1_bp_to_BC'])
R2_bp_to_bc = list(dem_df['R2_bp_to_BC'])

# this dict has keys for each library and values are the # of reads for that library
lib_counts = {l: 0 for l in libraries}
other_reads = 0  # reads not in any of the libraries

if len(R1_files) != len(R2_files):
    raise Exception("R1 and R2 file lists must be the same length")

output_files = dict()
for lib in libraries:
    output_files[lib] = [gzip.open(reads_out_base+run_name+'_'+lib+rep_suffix+'.fastq.gz', 'wt') for rep_suffix in ['_R1', '_R2']]
for r in range(len(R1_files)):
    R1in = R1_files[r]
    R2in = R2_files[r]
    if R1in[len(R1in)-3:] == '.gz':
        reads1 = gzip.open(R1in, 'rt')
    else:
        reads1 = open(R1in)
    if R2in[len(R2in)-3:] == '.gz':
        reads2 = gzip.open(R2in, 'rt')
    else:
        reads2 = open(R2in)
    rc = 0
    R2_iterator = FourLineFastq(reads2)
    for R1_title, R1_seq, R1_qual in FourLineFastq(reads1):
        R2_title, R2_seq, R2_qual = next(R2_iterator)
        rc += 1
        indices = R1_seq[8:8+R1_ind_len] + R2_seq[8:8+R2_ind_len]
        if indices in index_dict:   # demultiplexing
            tmp_lib_ind = index_dict[indices]
            tmp_lib = libraries[tmp_lib_ind]
            lib_counts[tmp_lib] += 1
            output_files[tmp_lib][0].write('\n'.join([R1_title, R1_seq, R1_qual, '+'])+'\n')
            output_files[tmp_lib][1].write('\n'.join([R2_title, R2_seq, R2_qual, '+'])+'\n')
        else:
            other_reads += 1  # increment counter for non-demultiplexed reads

reads1.close()
reads2.close()

print(rc, 'total reads')

for lib in libraries:
    print(lib, lib_counts[lib])
    output_files[lib][0].close()
    output_files[lib][1].close()

print('not demultiplexed', other_reads)

print('Time:', time.time()-otime)