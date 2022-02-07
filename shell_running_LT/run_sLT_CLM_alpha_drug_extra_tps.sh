#!/bin/bash
#SBATCH -J sLT_CLM_alpha-R1_redo  #job name for array
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-30:00              # Runtime in D-HH:MM
#SBATCH -p general       # Partition to submit to
#SBATCH --mem=60000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o sLT_CLM_alpha-R1_redo.out      # File to which STDOUT will be written
#SBATCH -e sLT_CLM_alpha-R1_redo.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent

source activate MILO_ENV

python3 ../PLT_parse.py ../demult_maps/Lineage_Tracking/Stanford_CLM_alpha-R1_redos_indices.csv ../../LT_data/ /n/regal/desai_lab/mjohnson/PLT/LT/2N_LT_171222_Stanford_drugs_R1.fastq.gz CLM_alpha_R1_redo

python3 ../PLT_cluster.py ../../LT_data/ CLM_alpha_R1_redo
