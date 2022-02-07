#!/bin/bash
#SBATCH -J BC_ext_FLC_2N_Plate3  #job name for array
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-10:00              # Runtime in D-HH:MM
#SBATCH -p serial_requeue       # Partition to submit to
#SBATCH --mem=3500               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ./shell_outs/bc_ext_%A_%a.out      # File to which STDOUT will be written
#SBATCH -e ./shell_outs/bc_ext_%A_%a.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent

module load Anaconda3/5.0.1-fasrc01
source activate milo_simple_conda5

python extract_bcs_V4.py ../FLC_2N_Plate3/sample_index.txt ../FLC_2N_Plate3 underscore bc_counts/FLC_2N_Plate3_ "${SLURM_ARRAY_TASK_ID}" 5
