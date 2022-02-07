import regex
import csv
import gzip
from collections import Counter

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('sample_index', help='file with samples listed on each line, these should just be missing the _R1_001.fastq.gz ending')
parser.add_argument('sample_directory', help='directory where fastq.gz files are found')
parser.add_argument('suffix_type', help='if the suffixes are like _R1_001.fastq.gz, input "underscore", if they are like .R1.fastq.gz input "dot"')
parser.add_argument('output_file_base', help='output file base directory')
parser.add_argument('job_index', help='for job array, tells which row of the file_list_file to start on')
parser.add_argument('job_index_step', help='how many files to parse in one job')
args = parser.parse_args()

### Reading in positional arguments
sample_index_file = args.sample_index
base_dir = args.sample_directory
output_file_base = args.output_file_base
job_index = int(args.job_index)
job_index_step = int(args.job_index_step)

if args.suffix_type == 'underscore':
    suffixes = ['_R1_001.fastq.gz', '_R2_001.fastq.gz']
elif args.suffix_type == 'dot':
    suffixes = ['.R1.fastq.gz', '.R2.fastq.gz']
else:
    print('Invalid input for suffix_type:', args.suffix_type, '. Must be "underscore" or "dot"')
    exit()

### Utility functions

def reverse_transcribe(seq):
    watson_crick = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    ls = len(seq)
    snew = ''
    for i in range(ls):
        if not seq[ls-i-1] in watson_crick.keys():
            snew += seq[ls-i-1]
        else:
            snew += watson_crick[seq[ls-i-1]]
    return snew

def FourLineFastq(handle):
    """
    Reads 4 lines from the file and returns 1, 2, and 4 with newlines stripped
    The only check for fastq format is that line 3 starts with '+'
    """
    while True:
        line = handle.readline()     
        if not line:
            # end of file
            break       
        title = line.rstrip()
        seq = handle.readline().rstrip()
        jnk_line = handle.readline()
        if jnk_line[0] != "+":
            print(title, seq, jnk_line)
            raise ValueError("Looks like this isn’t a strictly 4-line fastq file")
        qual = handle.readline().rstrip()
        yield (title, seq, qual)

### MAIN FUNCTIONS

def get_bcs(r):
    # Looks for barcodes in a read by searching the regex for the two barcodes on both the read and its reverse complement
    rt = reverse_transcribe(r)
    dmatch = re_dbc.search(r)
    dmatch_rt = re_dbc.search(rt)
    ematch = re_ebc.search(r)
    ematch_rt = re_ebc.search(rt)
    return dmatch, dmatch_rt, ematch, ematch_rt


def run_sample(samp, fout):
    print('parsing', samp)
    R1 = base_dir + '/' + samp + suffixes[0]
    R2 = base_dir + '/' + samp + suffixes[1]
    dbcs, ebcs = Counter(), Counter()
    r2in = gzip.open(R2, 'rt')
    r2_generator = FourLineFastq(r2in)
    r1in = gzip.open(R1, 'rt')
    c = 0
    for r1_title, r1_seq, r1_qual in FourLineFastq(r1in):
        r2_title, r2_seq, r2_qual = next(r2_generator)
        dbcs_this_read, ebcs_this_read = [], []  # we simply keep a list of all the dbcs and ebcs we find in this read pair
        for read in [r1_seq, r2_seq]:
            # search for barcodes by regex, these are: dbc regex, dbc regex on reverse comp, ebc regex, ebc regex on reverse comp
            dm, rc_dm, em, rc_em = get_bcs(read)
            if dm:
                dbcs_this_read.append(dm.group(2))
            elif rc_dm:
                dbcs_this_read.append(rc_dm.group(2))
            if em:
                ebcs_this_read.append(em.group(2))
            elif rc_em:
                ebcs_this_read.append(rc_em.group(2))
        # this avoids double counting if barcodes are read on the forward and reverse read, just taking the first one we find
        # will do error correction later when matching them up with the known barcodes
        if len(dbcs_this_read) > 0:
            dbcs[dbcs_this_read[0]] += 1
        if len(ebcs_this_read) > 0:
            ebcs[ebcs_this_read[0]] += 1
        c += 1

    # This is just a file we keep appending to with a table of filenames and the number of reads, to reference later
    with open('file_reads.txt', 'a') as outfile:
        outfile.write(R1 + '\t' + str(c) + '\n')
    
    # Writing output
    with open(fout, 'w') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['BC', 'Type', 'Reads'])
        for td, name in [(dbcs, 'DBC'), (ebcs, 'EBC')]:
            for bc in sorted([i for i in td], key=lambda x: td[x], reverse=True):
                writer.writerow([bc, name, td[bc]])
                  

## BEGIN MAIN CALLS ##
# reading sample_index file to get a list of sample names
samples_to_do = []
with open(sample_index_file, 'r') as infile:
    for line in infile:
        samples_to_do.append(line.strip('\n'))
        
re_dbc = regex.compile('(ATCGGTACC){e<=1}(\D{5}AA\D{5}AA\D{5}TT\D{5})(ATAACTTCG){e<=1}')
re_ebc = regex.compile('(ATCGGTACC){e<=1}(\D{5}AA\D{5}TT\D{5}TT\D{5})(ATAACTTCG){e<=1}')
    
for sample in samples_to_do[job_index*job_index_step:(job_index+1)*job_index_step]:
    run_sample(sample, output_file_base  + sample + '.csv')
