#!/bin/bash
#SBATCH -J hBFA2_merging  #job name for array
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-00:10              # Runtime in D-HH:MM
#SBATCH -p serial_requeue       # Partition to submit to
#SBATCH --mem=3500               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o hBFA2_merging.out      # File to which STDOUT will be written
#SBATCH -e hBFA2_merging.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent

source activate MILO_ENV

python ../PLT_merge_and_remove_chimeras.py ../../BFA_data/hBFA2_Harvard_1/hBFA2_Harvard_1_bc_counts_clustered.csv,../../BFA_data/hBFA2_Harvard_2/hBFA2_Harvard_2_bc_counts_clustered.csv,../../BFA_data/hBFA2_Stanford/hBFA2_Stanford_bc_counts_clustered.csv ../../BFA_data/Combined_Counts/hBFA2_counts.csv ../../BFA_data/Combined_Counts/Chimeric/hBFA2_chimeras.csv
