#!/bin/bash
#SBATCH -J dBFA2_Harvard_1  #job name for array
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-30:00              # Runtime in D-HH:MM
#SBATCH -p general       # Partition to submit to
#SBATCH --mem=60000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o dBFA2_Harvard_1.out      # File to which STDOUT will be written
#SBATCH -e dBFA2_Harvard_1.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent

source activate MILO_ENV

python3 ../PLT_parse.py ../demult_maps/BFAs/dBFA2_Harvard_1_indices.csv ../../BFA_data/ /n/regal/desai_lab/mjohnson/PLT/BFAs/BFA_raw_reads/dBFA2_Harvard1/dBFA2_Harvard1_R1.fastq.gz dBFA2_Harvard_1

python3 ../PLT_cluster.py ../../BFA_data/ dBFA2_Harvard_1

