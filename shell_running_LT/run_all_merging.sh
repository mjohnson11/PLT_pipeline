#!/bin/bash
#SBATCH -J LT_merging  #job name for array
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-10:00              # Runtime in D-HH:MM
#SBATCH -p serial_requeue       # Partition to submit to
#SBATCH --mem=10000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o LT_merging.out      # File to which STDOUT will be written
#SBATCH -e LT_merging.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent

source activate MILO_ENV

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/GlyEtOH_2N_R1/GlyEtOH_2N_R1_bc_counts_clustered.csv ../../LT_data/Combined_Counts/GlyEtOH_2N_R1_counts.csv ../../LT_data/Combined_Counts/Chimeric/GlyEtOH_2N_R1_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/GlyEtOH_alpha_R1/GlyEtOH_alpha_R1_bc_counts_clustered.csv ../../LT_data/Combined_Counts/GlyEtOH_alpha_R1_counts.csv ../../LT_data/Combined_Counts/Chimeric/GlyEtOH_alpha_R1_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/02M_NaCl_2N_R1/02M_NaCl_2N_R1_bc_counts_clustered.csv,../../LT_data/02M_NaCl_2N_R1_restart/02M_NaCl_2N_R1_restart_bc_counts_clustered.csv ../../LT_data/Combined_Counts/02M_NaCl_2N_R1_counts.csv ../../LT_data/Combined_Counts/Chimeric/02M_NaCl_2N_R1_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/02M_NaCl_alpha_R1/02M_NaCl_alpha_R1_bc_counts_clustered.csv ../../LT_data/Combined_Counts/02M_NaCl_alpha_R1_counts.csv ../../LT_data/Combined_Counts/Chimeric/02M_NaCl_alpha_R1_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/08M_NaCl_2N_R1/08M_NaCl_2N_R1_bc_counts_clustered.csv ../../LT_data/Combined_Counts/08M_NaCl_2N_R1_counts.csv ../../LT_data/Combined_Counts/Chimeric/08M_NaCl_2N_R1_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/08M_NaCl_alpha_R1/08M_NaCl_alpha_R1_bc_counts_clustered.csv ../../LT_data/Combined_Counts/08M_NaCl_alpha_R1_counts.csv ../../LT_data/Combined_Counts/Chimeric/08M_NaCl_alpha_R1_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/CLM_2N_R1/CLM_2N_R1_bc_counts_clustered.csv,../../LT_data/CLM_2N_R1_redo/CLM_2N_R1_redo_bc_counts_clustered.csv ../../LT_data/Combined_Counts/CLM_2N_R1_counts.csv ../../LT_data/Combined_Counts/Chimeric/CLM_2N_R1_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/CLM_alpha_R1/CLM_alpha_R1_bc_counts_clustered.csv,../../LT_data/CLM_alpha_R1_redo/CLM_alpha_R1_redo_bc_counts_clustered.csv ../../LT_data/Combined_Counts/CLM_alpha_R1_counts.csv ../../LT_data/Combined_Counts/Chimeric/CLM_alpha_R1_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/FLC4_2N_R1/FLC4_2N_R1_bc_counts_clustered.csv,../../LT_data/FLC4_2N_R1_redo/FLC4_2N_R1_redo_bc_counts_clustered.csv ../../LT_data/Combined_Counts/FLC4_2N_R1_counts.csv ../../LT_data/Combined_Counts/Chimeric/FLC4_2N_R1_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/FLC4_alpha_R1/FLC4_alpha_R1_bc_counts_clustered.csv,../../LT_data/FLC4_alpha_R1_redo/FLC4_alpha_R1_redo_bc_counts_clustered.csv ../../LT_data/Combined_Counts/FLC4_alpha_R1_counts.csv ../../LT_data/Combined_Counts/Chimeric/FLC4_alpha_R1_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/FLC32_2N_R1/FLC32_2N_R1_bc_counts_clustered.csv ../../LT_data/Combined_Counts/FLC32_2N_R1_counts.csv ../../LT_data/Combined_Counts/Chimeric/FLC32_2N_R1_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/FLC32_alpha_R1/FLC32_alpha_R1_bc_counts_clustered.csv ../../LT_data/Combined_Counts/FLC32_alpha_R1_counts.csv ../../LT_data/Combined_Counts/Chimeric/FLC32_alpha_R1_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/M3_2N_R1/M3_2N_R1_bc_counts_clustered.csv ../../LT_data/Combined_Counts/M3_2N_R1_counts.csv ../../LT_data/Combined_Counts/Chimeric/M3_2N_R1_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/M3_alpha_R1/M3_alpha_R1_bc_counts_clustered.csv ../../LT_data/Combined_Counts/M3_alpha_R1_counts.csv ../../LT_data/Combined_Counts/Chimeric/M3_alpha_R1_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/21C_2N_R1/21C_2N_R1_bc_counts_clustered.csv,../../LT_data/21C_2N_R1_restart/21C_2N_R1_restart_bc_counts_clustered.csv ../../LT_data/Combined_Counts/21C_2N_R1_counts.csv ../../LT_data/Combined_Counts/Chimeric/21C_2N_R1_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/21C_alpha_R1/21C_alpha_R1_bc_counts_clustered.csv ../../LT_data/Combined_Counts/21C_alpha_R1_counts.csv ../../LT_data/Combined_Counts/Chimeric/21C_alpha_R1_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/37C_2N_R1/37C_2N_R1_bc_counts_clustered.csv ../../LT_data/Combined_Counts/37C_2N_R1_counts.csv ../../LT_data/Combined_Counts/Chimeric/37C_2N_R1_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/37C_2N_R2/37C_2N_R2_bc_counts_clustered.csv ../../LT_data/Combined_Counts/37C_2N_R2_counts.csv ../../LT_data/Combined_Counts/Chimeric/37C_2N_R2_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/YPD_2N_R1/YPD_2N_R1_bc_counts_clustered.csv ../../LT_data/Combined_Counts/YPD_2N_R1_counts.csv ../../LT_data/Combined_Counts/Chimeric/YPD_2N_R1_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/YPD_2N_R2/YPD_2N_R2_bc_counts_clustered.csv ../../LT_data/Combined_Counts/YPD_2N_R2_counts.csv ../../LT_data/Combined_Counts/Chimeric/YPD_2N_R2_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/SC_2N_R1/SC_2N_R1_bc_counts_clustered.csv ../../LT_data/Combined_Counts/SC_2N_R1_counts.csv ../../LT_data/Combined_Counts/Chimeric/SC_2N_R1_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/SC_2N_R2/SC_2N_R2_bc_counts_clustered.csv ../../LT_data/Combined_Counts/SC_2N_R2_counts.csv ../../LT_data/Combined_Counts/Chimeric/SC_2N_R2_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/48Hr_2N_R1/48Hr_2N_R1_bc_counts_clustered.csv ../../LT_data/Combined_Counts/48Hr_2N_R1_counts.csv ../../LT_data/Combined_Counts/Chimeric/48Hr_2N_R1_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/48Hr_2N_R2/48Hr_2N_R2_bc_counts_clustered.csv ../../LT_data/Combined_Counts/48Hr_2N_R2_counts.csv ../../LT_data/Combined_Counts/Chimeric/48Hr_2N_R2_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/pH7_3_2N_R1/pH7_3_2N_R1_bc_counts_clustered.csv ../../LT_data/Combined_Counts/pH7_3_2N_R1_counts.csv ../../LT_data/Combined_Counts/Chimeric/pH7_3_2N_R1_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/pH7_3_2N_R2/pH7_3_2N_R2_bc_counts_clustered.csv ../../LT_data/Combined_Counts/pH7_3_2N_R2_counts.csv ../../LT_data/Combined_Counts/Chimeric/pH7_3_2N_R2_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/pH3_8_2N_R1/pH3_8_2N_R1_bc_counts_clustered.csv ../../LT_data/Combined_Counts/pH3_8_2N_R1_counts.csv ../../LT_data/Combined_Counts/Chimeric/pH3_8_2N_R1_chimeras.csv -merge_type outer

python ../PLT_merge_and_remove_chimeras.py ../../LT_data/pH3_8_2N_R2/pH3_8_2N_R2_bc_counts_clustered.csv ../../LT_data/Combined_Counts/pH3_8_2N_R2_counts.csv ../../LT_data/Combined_Counts/Chimeric/pH3_8_2N_R2_chimeras.csv -merge_type outer
