#!/bin/bash
#SBATCH -J demult_S_hap  #job name for array
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-60:00              # Runtime in D-HH:MM
#SBATCH -p shared       # Partition to submit to
#SBATCH --mem=60000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o output/demult_S_hap.out      # File to which STDOUT will be written
#SBATCH -e output/demult_S_hap.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent


python3 ../PLT_demult.py ../demult_maps/Lineage_Tracking/Stanford_hap_all.csv /n/holyscratch01/desai_lab/mjohnson/PLT/demult_reads/LT/ /n/holyscratch01/desai_lab/mjohnson/PLT/raw_reads/LT/alpha_LT_Stanford_R1.fastq.gz Stanford_hap_LT

