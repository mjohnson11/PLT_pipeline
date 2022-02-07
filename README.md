# Scripts used to process data from bulk barcode-based fitness assays in the pleiotropy lineage tracking project

## 1) Extracting and error-correcting barcodes, merging and cleaning counts, adding home environment info

* PLT_parse.py, then PLT_cluster.py, then PLT_merge_and_remove_chimeras.py, and PLT_BFA_env_calling do this

* see shell scripts in shell_running_BFAs, and shell_running_BFA_subpools for the specific jobs submitted

* Note: for the BFA subpools, we use PLT_parse_new_primers.py instead (bc of new offset primer system)

## 2) Extracting barcodes from whole genome sequencing data

* extract_bcs_V4.py and compile_bc_extractions.py do this

* see shell scripts in shell_running_BC_extraction

## 3) Final Analyses

* S_measure.ipynb infers fitnesses from count data

* Mutation_parsing.ipynb parses mutation data and combines it with fitness data

* AT_bias_demonstration.ipynb shows the methods to create the supplemental figure showing AT/GC-content biases

* Figs.ipynb has some preliminary figures

