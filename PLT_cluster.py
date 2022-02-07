"""
For clustering PLT BC counts
"""

import time
import csv
from DeletionErrorCorrector import DeletionClusterator
import pandas as pd

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('out_base', help='output file directory, should be the same as in PLT_parse, this is where'
                                     'we will be reading the unclustered bc counts file from')
parser.add_argument('run_name', help='run name, should be the same as in PLT_parse')
parser.add_argument("-edit_dist_thresh", type=int, default=3, help='threshold for edit distance clustering within UMI families, default=3')
parser.add_argument("-centroid_count_thresh", type=int, default=10, help='threshold of counts for a centroid to be considered real, default=10')

args = parser.parse_args()

# POSITIONAL ARGS
out_base = args.out_base
run_name = args.run_name

output_dir = out_base + run_name + '/'
# input file: should have columns with headers that say "Diverse.BC", "Environment.BC", and "Total.Counts"
bc_csv = output_dir + run_name + '_bc_counts.csv'
# cluster output file
out_csv = output_dir + run_name + '_clusters.csv'
# output base for cluster info
out_base = output_dir + run_name

# OPTIONAL ARGS

# max edits for errors to cluster
EDITDIST_THRESH = args.edit_dist_thresh

# potential centroids with less than or equal to this number of total reads are not included
CENTROID_COUNT_THRESH = args.centroid_count_thresh


class ClusterParser:
    """
    ClusterParser is a class used to parse a cluster/error-correction output file and output various bc count files
    The input files are
    1) (f_in) the cluster file with rows like
    centroid index, index that clustered to that centroid, another index that clustered to that centroid, etc.
    2) (bc_in) the original barcode counts file used for error correction, with columns for counts and columns
    for Diverse.BC, Environment.BC, Total.Counts
    """

    def __init__(self, f_in, bc_in):

        print('\nReading barcode count data file:', bc_in)
        self.bcdata = pd.read_csv(bc_in)
        self.cols_to_add = [i for i in self.bcdata.select_dtypes(include=['int']).keys()]
        self.info_cols = [i for i in self.bcdata.keys() if i not in self.cols_to_add]
        self.info_matrix = self.bcdata.as_matrix(self.info_cols)
        self.int_matrix = self.bcdata.as_matrix(self.cols_to_add)
        self.bcdata['Index'] = range(len(self.bcdata))  # make index column
        self.total_counts_col = list(self.bcdata['Total.Counts'])
        print('Read', len(self.bcdata), 'rows.')
        self.total_initial_barcodes = len(self.bcdata)
        self.num_low_count_centroids = 0
        self.low_count_centroids_reads = 0
        print('\nReading cluster file:', f_in)
        self.cluster_rows = []
        self.centroids_dict = dict()
        self.error_cor_dict = dict()
        self.all_clustered = set()
        with open(f_in, 'r') as infile:
            reader = csv.reader(infile)
            for row in reader:
                row = [int(i) for i in row]
                centroid_reads = self.total_counts_col[row[0]]
                # final check that the diverse bc + environment bc combo for a centroid has more reads than the threshold
                if centroid_reads <= CENTROID_COUNT_THRESH:
                    self.num_low_count_centroids += 1
                    self.low_count_centroids_reads += centroid_reads
                else:
                    self.cluster_rows.append(row)
                    self.centroids_dict[row[0]] = set(row)
                    self.all_clustered.update(row)
                    for index in row:
                        # note that every centroid corrects to itself
                        self.error_cor_dict[index] = row[0]

        # I will keep track of when I see "error" barcodes at > 1% the frequency of the real barcode
        # this is possible for mono-nucleotide run errors, but in other cases it could be an error correction mistake
        self.error_percent_too_high = 0
        self.ep_thresh = 0.01  # threshold for calling an error weird

    def output_centroids_without_errors(self, f_out):
        print('\nOutputting centroid barcode rows (no errors) to:', f_out)
        self.bcdata.loc[self.bcdata['Index'].isin(self.centroids_dict)].to_csv(f_out, index=False)
        print('Wrote', len(self.bcdata.loc[self.bcdata['Index'].isin(self.centroids_dict)]), 'rows.')

    def output_clusters(self, f_out):
        print('\nOutputting clustered lineages to:', f_out)

        titlerow = self.info_cols + self.cols_to_add
        total_reads_index = self.cols_to_add.index('Total.Counts')
        total_reads_clust = 0
        total_bc_clust = 0
        lcta = len(self.cols_to_add)
        with open(f_out, 'w') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(titlerow)
            for centroid in sorted(self.centroids_dict.keys()):
                tmp_row = self.int_matrix[centroid]
                for index in self.centroids_dict[centroid]:
                    if not index == centroid:
                        # add read counts to centroid row
                        tmp_counts = self.int_matrix[index]
                        for i in range(lcta):
                            tmp_row[i] += tmp_counts[i]

                writer.writerow(list(self.info_matrix[centroid]) + list(tmp_row))
                total_reads_clust += tmp_row[total_reads_index]
                total_bc_clust += len(self.centroids_dict[centroid])

        print('Wrote', len(self.centroids_dict), 'rows.')
        return total_bc_clust, total_reads_clust

    def output_cluster_stats(self, f_out):
        print('\nOutputting stats on clusters to:', f_out)
        with open(f_out, 'w') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(self.info_cols + ['Number above ' + str(self.ep_thresh*100) + ' % of centroid reads'])
            for centroid in sorted(self.centroids_dict.keys()):
                centroid_reads = self.total_counts_col[centroid]
                ep_too_high = 0  # error percent too high counts
                for index in self.centroids_dict[centroid]:
                    if not index == centroid:
                        if (self.total_counts_col[index] / centroid_reads) > self.ep_thresh:
                            ep_too_high += 1
                tmp_row = list(self.info_matrix[centroid]) + [ep_too_high]
                writer.writerow(tmp_row)
                self.error_percent_too_high += ep_too_high

        print('Found', self.error_percent_too_high, 'clustering events in which the clustered barcode had more than',
              self.ep_thresh*100, '% of the total reads of the centroid barcode.')

    def output_unclustered(self, f_out):
        print('\nOutputting unclustered barcodes to:', f_out)
        unclust_data = self.bcdata.loc[~self.bcdata['Index'].isin(self.all_clustered)]
        unclust_data.to_csv(f_out, index=False)
        num_unclust_rows = len(unclust_data)
        unclust_reads = sum(unclust_data['Total.Counts'])
        print('Wrote', num_unclust_rows, 'rows, comprising', unclust_reads, 'reads.')
        return num_unclust_rows, unclust_reads


# Main function

def div_bc_clustering(ec, dat):
    div_bcs_list = dat.as_matrix(['Diverse.BC', 'Total.Counts'])
    print('Correcting Div BCs...')
    div_corrector = DeletionClusterator(div_bcs_list)
    div_corrector.get_deletion_clusters(True)
    div_corrector.cluster_within_delnets(EDITDIST_THRESH, CENTROID_COUNT_THRESH, False)
    full_bc_list = dat.as_matrix(['Diverse.BC', 'Environment.BC', 'Total.Counts'])
    clusters_by_index = dict()
    odd_fail = 0  # counts the odd case where the most common thing that corrects to a ebc dbc pair is not that actual pair
    for r in range(len(full_bc_list)):
        row = full_bc_list[r]
        ebc = row[1]
        dbc = row[0]
        if ebc in ec.corrector:
            real_ebc = ec.corrector[ebc]
            if dbc in div_corrector.corrector:
                real_dbc = div_corrector.corrector[dbc]
                bc_combo = real_dbc + real_ebc
                if bc_combo in clusters_by_index:
                    clusters_by_index[bc_combo].append(r)
                else:
                    if ebc == real_ebc and dbc == real_dbc:
                        clusters_by_index[bc_combo] = [r]
                    else:
                        odd_fail += 1

    total_clusters = len(clusters_by_index)
    print('Made', total_clusters, 'total clusters.', odd_fail, 'failed because the top hit did not match the true bcs.')

    # PART 2: going through the original file and correcting, outputting a file with row indices like:
    # bc index, bc index, bc index, etc.
    """
    This is the main error-correction/clustering output file
    It has rows like:
    centroid index, index that clustered to that centroid, another index that clustered to that centroid, etc.
    """
    print('Writing error correction cluster output to:', out_csv)
    with open(out_csv, 'w') as outfile:
        writer = csv.writer(outfile)
        for bc_combo in sorted(clusters_by_index.keys(), key=lambda x: clusters_by_index[x][0]):
            writer.writerow(clusters_by_index[bc_combo])

    print('Parsing and outputting clusters now')

    bc_out = out_base + '_bc_counts_clustered.csv'
    centroids_only_out = out_base + '_centroids_only.csv'
    cluster_details = out_base + '_cluster_info.csv'
    unclustered_out = out_base + '_unclustered.csv'
    my_clus_pars = ClusterParser(out_csv, bc_csv)
    unclustered_bcs, unclustered_reads = my_clus_pars.output_unclustered(unclustered_out)
    my_clus_pars.output_centroids_without_errors(centroids_only_out)
    clustered_bcs, clustered_reads = my_clus_pars.output_clusters(bc_out)
    my_clus_pars.output_cluster_stats(cluster_details)

    print("Looked at", my_clus_pars.total_initial_barcodes, "barcodes")
    print("Made", len(my_clus_pars.centroids_dict), "clusters from", clustered_bcs, "bcs, totaling",
          clustered_reads, 'reads.')
    print("Could not match", unclustered_bcs, "bcs, totaling", unclustered_reads, 'reads.')

    with open(out_base + '_error_correct_stats.csv', 'w') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['Run', 'Raw.BCs', 'Raw.Reads', 'Centroids', 'Clustered.BCs', 'Clustered.Reads',
                         'Unclustered.BCs', 'Unclustered.Reads', 'Centroids.Excluded.On.Thresh',
                         'Reads.Excluded.On.Thresh'])
        writer.writerow([run_name, len(dat), sum(dat['Total.Counts']), len(my_clus_pars.centroids_dict),
                         clustered_bcs, clustered_reads, unclustered_bcs, unclustered_reads,
                         my_clus_pars.num_low_count_centroids, my_clus_pars.low_count_centroids_reads])



#
# MAIN CALLS
#
otime = time.time()

bcdata = pd.read_csv(bc_csv)
env_bcs_list = bcdata.as_matrix(['Environment.BC', 'Total.Counts'])
print('Correcting Env BCs...')
env_corrector = DeletionClusterator(env_bcs_list)
env_corrector.get_deletion_clusters(True)
env_corrector.cluster_within_delnets(EDITDIST_THRESH, CENTROID_COUNT_THRESH, False)
# div_bc_clustering_by_env_bc(env_corrector, bcdata)
div_bc_clustering(env_corrector, bcdata)
print('Done, time elapsed:', time.time()-otime)

