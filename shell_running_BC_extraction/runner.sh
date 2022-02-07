#sbatch --array=0-14 run_CLM_2N_Batch2.sh

sbatch --array=0-38 run_CLM_clones_1N.sh

sbatch --array=0-33 run_FLC_R1_clones_1N.sh

sbatch --array=0-17 run_FLC_2N_Plate3.sh

sbatch --array=0-38 run_FLC_R2_clones_1N.sh

sbatch --array=0-19 run_GlyEtOH_2N_20170616.sh

sbatch --array=0-23 run_GlyEtOH_1N.sh

sbatch --array=0-8 run_CLM_2N.sh

sbatch --array=0-18 run_FLC_2N_Plate2.sh

sbatch --array=0-18 run_FLC_2N_Plate1.sh

sbatch --array=0-19 run_CLM_R2_clones_1N.sh

sbatch --array=0-14 run_Dip_clone_Env.sh

