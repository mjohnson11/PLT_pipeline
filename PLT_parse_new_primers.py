"""
A program to demultiplex reads and count barcodes for the PLT
editing on 3/20/18 for the new primers and single-end sequencing
"""

import time
import gzip
import pandas as pd
import re
import csv
import numpy as np
import os
import subprocess
from PLT_tools import FourLineFastq, reverse_transcribe
from glob import glob


import argparse
parser = argparse.ArgumentParser()
parser.add_argument('primer_map_file', help='file containing primer demultiplex info')
parser.add_argument('out_base', help='basic output directory')
parser.add_argument('input_base', help='directory with fastq.gz files')
parser.add_argument('run_name', help='name for the run')
parser.add_argument("-quality_cutoff", type=int, default=30, help='quality threshold for bc region')
args = parser.parse_args()

# POSITIONAL ARGS - THESE SHOULD ALWAYS BE THE SAME FOR A RUN
primer_map_file = args.primer_map_file
input_base = args.input_base
out_base = args.out_base
run_name = args.run_name
QUALITY_CUTOFF = args.quality_cutoff
output_dir = out_base + run_name + '/'
output_file = output_dir + run_name + '_bc_counts.csv'
stats_out_base = output_dir + 'run_statistics/'
lib_stats_out = stats_out_base + run_name + '_library_statistics.csv'

if not os.path.isdir(output_dir):
    print('Making main output directory:', output_dir)
    subprocess.call(['mkdir', output_dir])
if not os.path.isdir(stats_out_base):
    print('Making stats output directory:', stats_out_base)
    subprocess.call(['mkdir', stats_out_base])

MY_REGEX = re.compile('\D*?(GTACC|GGACC|GGTCC|G.TACC|GG.ACC|GGT.CC|GGTA.C|GGTAC.)(\D{24,28})(.TAACT|A.AACT|AT.ACT|ATA.CT|ATAA.T|ATAAC|AAACT|ATACT|ATAAT)\D*')
EBC_REGEX = re.compile('\D*?(ATTAT|AGTAT|AGTTT|GTTAT|A.TTAT|AG.TAT|AGT.AT|AGTT.T|AGTTA.)(\D{24,28})(.GTACC|G.TACC|GG.ACC|GGT.CC|GGTA.C|GGACC|GGTCC|GGTAC)\D*')


class BcCounter:

    """
    BcCounter counts barcodes and corrects counts based on the unique molecular indices (UMIs) they are paired with.
    """

    def __init__(self, demultiplex_file):
        # this bc_dict has entries like: [bc div, bc env, total counts, [counts_in_lib_1, counts_in_lib_2, etc.]]
        self.bc_dict = dict()

        self.dem_df = pd.read_csv(demultiplex_file)

        self.libraries = list(self.dem_df['Library'])
        self.num_libraries = len(self.libraries)
        self.indices = list(self.dem_df['Inline_index'])
        self.R1_bp_to_bc = list(self.dem_df['R1_bp_to_BC'])

        # this dict has keys for each library, and entries like
        # [total reads, failed on quality, failed on regex, primer dimer likely, failed on inline index]
        self.lib_stats = {l: [0, 0, 0, 0, 0] for l in self.libraries}
        self.other_reads = 0  # reads not in any of the libraries

    def add_count(self, bc_div, bc_env, library_index):
        bc = bc_div + bc_env
        if bc in self.bc_dict:
            tmp_entry = self.bc_dict[bc]
            tmp_entry[2] += 1
            tmp_entry[3][library_index] += 1
        else:
            tmp_entry = self.bc_dict[bc] = [bc_div, bc_env, 1, [0]*self.num_libraries]
            tmp_entry[3][library_index] += 1

    def read_file(self, library_index):
        infile_list = glob(input_base + self.libraries[library_index] + '_S*_R1_001.fastq.gz')
        assert len(infile_list) == 1
        R1in = infile_list[0]
        print('Reading', R1in)
        index_expected = self.indices[library_index]
        r1_dbc_start = self.R1_bp_to_bc[library_index]
        tmp_lib_stats = self.lib_stats[self.libraries[library_index]]
        if R1in[len(R1in)-3:] == '.gz':
            reads1 = gzip.open(R1in, 'rt')
        else:
            reads1 = open(R1in)
        rc = 0
        for R1_title, R1_seq, R1_qual in FourLineFastq(reads1):
            rc += 1
            index = R1_seq[:8]
            if index == index_expected:   # demultiplexing
                r1_ebc_start = r1_dbc_start + 60
                tmp_lib_stats[0] += 1   # count for total reads
                # quality check
                if np.mean([ord(c)-33 for c in R1_qual[r1_dbc_start:r1_dbc_start+26] + R1_qual[r1_ebc_start:r1_ebc_start+26]]) < QUALITY_CUTOFF:
                    tmp_lib_stats[1] += 1
                else:
                    # regex check
                    reghit1 = MY_REGEX.match(R1_seq[r1_dbc_start-10:r1_dbc_start+36])
                    reghit2 = EBC_REGEX.match(R1_seq[r1_ebc_start-10:r1_ebc_start+36])
                    if (not reghit1) or (not reghit2):
                        # regex failed
                        tmp_lib_stats[2] += 1
                        # checks if regex failed because this was a primer-dimer fragment
                        if 'CTGTCTCTT' in R1_seq:
                            tmp_lib_stats[3] += 1
                    else:
                        bc_div = reghit1.group(2)
                        bc_env = reverse_transcribe(reghit2.group(2))
                        self.add_count(bc_div, bc_env, library_index)
            else:
                tmp_lib_stats[4] += 1  # increment counter for non-demultiplexed reads

        print('Read', rc, 'reads,', rc - tmp_lib_stats[4], 'matched library inline index.')

        reads1.close()

    def write_output(self, fout):

        with open(fout, 'w') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(['Diverse.BC', 'Environment.BC', 'Total.Counts'] + self.libraries)
            sorted_bcs = sorted(self.bc_dict, key=lambda x: self.bc_dict[x][2], reverse=True)
            for bc in sorted_bcs:
                entry = self.bc_dict[bc]
                writer.writerow(entry[:3]+entry[3])

    def write_lib_stats(self, fout):

        with open(fout, 'w') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(['Library', 'Total.Reads', 'Quality.Failed', 'Regex.Failed', 'Primer.Dimer',
                             'Inline.Index.Wrong', 'Usable.Reads'])
            for lib in self.libraries:
                entry = self.lib_stats[lib]
                writer.writerow([lib] + entry + [entry[0]-entry[1]-entry[2]-entry[4]])


otime = time.time()
bcc = BcCounter(primer_map_file)
for lib_index in range(len(bcc.libraries)):
    bcc.read_file(lib_index)

bcc.write_lib_stats(lib_stats_out)

bcc.write_output(output_file)

print('Time:', time.time()-otime)