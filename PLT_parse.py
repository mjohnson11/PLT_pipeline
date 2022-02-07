"""
A program to demultiplex reads and count barcodes for the PLT
Milo Johnson
Started 3_6_16
made simple 2_6_17
"""

import time
import gzip
import pandas as pd
import re
import csv
import numpy as np
import os
import subprocess
from PLT_htmlReport import make_html_report
from PLT_tools import FourLineFastq


import argparse
parser = argparse.ArgumentParser()
parser.add_argument('primer_map_file', help='file containing primer demultiplex info')
parser.add_argument('out_base', help='basic output directory')
parser.add_argument('R1_files', help='comma separated list of R1 files to parse - '
                                     '"R1" will be changed to "R2" to get the paired read files')
parser.add_argument('run_name', help='name for the run')
parser.add_argument("-quality_cutoff", type=int, default=30, help='quality threshold for bc region')
args = parser.parse_args()

# POSITIONAL ARGS - THESE SHOULD ALWAYS BE THE SAME FOR A RUN
primer_map_file = args.primer_map_file
R1_files = args.R1_files.split(',')
R2_files = [i.replace('R1', 'R2') for i in R1_files]
out_base = args.out_base
run_name = args.run_name
QUALITY_CUTOFF = args.quality_cutoff
output_dir = out_base + run_name + '/'
output_file = output_dir + run_name + '_bc_counts.csv'
stats_out_base = output_dir + 'run_statistics/'
lib_stats_out = stats_out_base + run_name + '_library_statistics.csv'
umi_fam_size_out = stats_out_base + run_name + '_umi_family_sizes.csv'
html_report_out = stats_out_base + run_name + '_html_report'

if not os.path.isdir(output_dir):
    print('Making main output directory:', output_dir)
    subprocess.call(['mkdir', output_dir])
if not os.path.isdir(stats_out_base):
    print('Making stats output directory:', stats_out_base)
    subprocess.call(['mkdir', stats_out_base])

MY_REGEX = re.compile('\D*?(GTACC|GGACC|GGTCC|G.TACC|GG.ACC|GGT.CC|GGTA.C|GGTAC.)(\D{24,28})(.TAACT|A.AACT|AT.ACT|ATA.CT|ATAA.T|ATAAC|AAACT|ATACT|ATAAT)\D*')


class BcCounter:

    """
    BcCounter counts barcodes and corrects counts based on the unique molecular indices (UMIs) they are paired with.
    """

    def __init__(self, demultiplex_file):
        # this bc_dict has entries like: [bc div, bc env, total counts, [counts_in_lib_1, counts_in_lib_2, etc.]]
        self.bc_dict = dict()

        self.dem_df = pd.read_csv(demultiplex_file)
        R1_index_lens = list(set([len(ind) for ind in self.dem_df['R1_index']]))
        R2_index_lens = list(set([len(ind) for ind in self.dem_df['R2_index']]))
        if len(R1_index_lens) != 1:
            raise Exception('R1 indices are not all the same length')
        if len(R2_index_lens) != 1:
            raise Exception('R2 indices are not all the same length')
        self.R1_ind_len = R1_index_lens[0]
        self.R2_ind_len = R2_index_lens[0]

        self.libraries = list(self.dem_df['Library'])
        self.num_libraries = len(self.libraries)
        index_list = self.dem_df.as_matrix(['R1_index', 'R2_index'])
        self.index_dict = {index_list[i][0]+index_list[i][1]: i for i in range(len(index_list))}
        self.R1_bp_to_bc = list(self.dem_df['R1_bp_to_BC'])
        self.R2_bp_to_bc = list(self.dem_df['R2_bp_to_BC'])

        # this dict has keys for each library, and entries like
        # [total reads, failed on quality, failed on regex, failed on UMI, primer dimer likely]
        self.lib_stats = {l: [0, 0, 0, 0, 0] for l in self.libraries}
        self.other_reads = 0  # reads not in any of the libraries
        # this dict has keys that are the libraries and values that are dictionaries
        # these dictionaries have concatenated UMIs and entries that are counts in that UMI family
        self.umi_dict = {l: dict() for l in self.libraries}

    def add_count(self, bc_div, bc_env, library_index):
        bc = bc_div + bc_env
        if bc in self.bc_dict:
            tmp_entry = self.bc_dict[bc]
            tmp_entry[2] += 1
            tmp_entry[3][library_index] += 1
        else:
            tmp_entry = self.bc_dict[bc] = [bc_div, bc_env, 1, [0]*self.num_libraries]
            tmp_entry[3][library_index] += 1

    def read_files(self, R1_files, R2_files):
        if len(R1_files) != len(R2_files):
            raise Exception("R1 and R2 file lists must be the same length")
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
                indices = R1_seq[8:8+self.R1_ind_len] + R2_seq[8:8+self.R2_ind_len]
                if indices in self.index_dict:   # demultiplexing
                    tmp_lib_ind = self.index_dict[indices]
                    r1_bc_start = self.R1_bp_to_bc[tmp_lib_ind]
                    r2_bc_start = self.R2_bp_to_bc[tmp_lib_ind]
                    tmp_lib = self.libraries[tmp_lib_ind]
                    tmp_lib_stats = self.lib_stats[self.libraries[tmp_lib_ind]]
                    tmp_lib_stats[0] += 1   # count for total reads
                    # quality check
                    if np.mean([ord(c)-33 for c in (R1_qual[r1_bc_start:r1_bc_start+26] +
                                                    R2_qual[r2_bc_start:r2_bc_start+26])]) < QUALITY_CUTOFF:
                        tmp_lib_stats[1] += 1
                    else:
                        # regex check
                        reghit1 = MY_REGEX.match(R1_seq[r1_bc_start-10:r1_bc_start+36])
                        reghit2 = MY_REGEX.match(R2_seq[r2_bc_start-10:r2_bc_start+36])
                        if (not reghit1) or (not reghit2):
                            # regex failed
                            tmp_lib_stats[2] += 1
                            # checks if regex failed because this was a primer-dimer fragment
                            if 'TTGAATTCGA' in R1_seq:
                                tmp_lib_stats[4] += 1
                        else:
                            # umi check
                            umi_combined = R1_seq[:8] + R2_seq[:8]
                            if umi_combined in self.umi_dict[tmp_lib]:
                                self.umi_dict[tmp_lib][umi_combined] += 1
                                tmp_lib_stats[3] += 1
                            else:
                                self.umi_dict[tmp_lib][umi_combined] = 1
                                # all checks passed, group 2 in the regex match is the barcode region
                                bc_div = reghit1.group(2)
                                bc_env = reghit2.group(2)
                                self.add_count(bc_div, bc_env, tmp_lib_ind)
                else:
                    self.other_reads += 1  # increment counter for non-demultiplexed reads

        print('Read', rc, 'reads,', rc - self.other_reads, 'matched library indices.')

        reads1.close()
        reads2.close()

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
            writer.writerow(['Library', 'Total.Reads', 'Quality.Failed', 'Regex.Failed', 'UMI.Repeats', 'Primer.Dimer',
                             'Usable.Reads'])
            for lib in self.libraries:
                entry = self.lib_stats[lib]
                writer.writerow([lib] + entry + [entry[0]-entry[1]-entry[2]-entry[3]])

    def write_umi_fam_sizes(self, fout):
        # outputs rows with umi family size distributions for all libraries in the run
        biggest_fam = max([max(self.umi_dict[l].values()) for l in self.libraries if len(self.umi_dict[l]) > 0])
        with open(fout, 'w') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(['Library'] + [str(i) for i in range(1, biggest_fam + 1)])
            for lib in self.libraries:
                tmp_row = [lib]
                lib_fam_sizes = list(self.umi_dict[lib].values())
                tmp_row += [lib_fam_sizes.count(i) for i in range(1, biggest_fam + 1)]
                writer.writerow(tmp_row)


otime = time.time()
bcc = BcCounter(primer_map_file)

bcc.read_files(R1_files, R2_files)

bcc.write_lib_stats(lib_stats_out)

bcc.write_output(output_file)

bcc.write_umi_fam_sizes(umi_fam_size_out)

make_html_report(lib_stats_out, umi_fam_size_out, html_report_out, run_name)

print('Time:', time.time()-otime)