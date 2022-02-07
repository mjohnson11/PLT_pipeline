import pandas as pd
import numpy as np

# Function to call the environment based on subpool sequencing, criteria: 
# 95% of the total subpool counts have to be from one environment
# This has to be at least 3 reads
def call_env(row, env_list):
    env_sums = {env: sum([row[c] for c in row.keys() if env in c.split('-')]) for env in env_list}
    max_counts = max(list(env_sums.values()))
    if max_counts/row['Total.Counts'] > 0.95 and max_counts > 2:
        max_env = [env for env in env_list if env_sums[env] == max_counts]
        assert len(max_env) == 1
        return max_env[0]
    else:
        return 'none'
    
# Simple function to determine if the bc is present in the R1 subpool, the R2 subpool, or both
def call_env_reps(row):
    if row['Subpool.Environment'] == 'none':
        return 'none'
    else:
        env = row['Subpool.Environment'] 
        env_reps = [c for c in row.keys() if env in c.split('-') if row[c]/row['Total.Counts'] >= 0.05]
        return '.'.join([e[e.index(env)+len(env):] for e in env_reps])


bfas = ['hBFA1', 'hBFA2', 'dBFA2']

# Reading count files
dats = dict()
for b in bfas:
    dats[b] = pd.read_csv('../../BFA_data/Combined_Counts/' + b + '_counts.csv')
    
# Reading this new subpool sequencing
subp = dict()
for b in bfas:
    subp[b] = pd.read_csv('../../BFA_Subpool_data/' + b + '_subpools/' + b + '_subpools_bc_counts_clustered.csv')
    subp[b]['Full.BC'] = subp[b]['Diverse.BC'] + subp[b]['Environment.BC']


# For dBFA2, many of the subpools were sequenced and processed along with the BFA, so I will add those in:
dbfa2_extra = pd.read_csv('../../BFA_data/dBFA2_Harvard_1/dBFA2_Harvard_1_bc_counts_clustered.csv')
renamer = {'Ancestor_YPD_2N_R1_1': 'dBFA2-Ancestor_YPD_2N-R1-1', 'Ancestor_YPD_2N_R1_2': 'dBFA2-Ancestor_YPD_2N-R1-2'}
for c in [i for i in dbfa2_extra.columns if 'Subpool' in i]:
    renamer[c] = 'dBFA2-' + c[8:].replace('_R', '-R') + '-Subpool'
dbfa2_extra.rename(index=str, columns=renamer, inplace=True)
dbfa2_extra['Full.BC'] = dbfa2_extra['Diverse.BC'] + dbfa2_extra['Environment.BC']
cols = ['Diverse.BC', 'Environment.BC', 'Full.BC']
subp['dBFA2'] = subp['dBFA2'].merge(dbfa2_extra[cols + list(renamer.values())], on=cols, how='outer').fillna(0)

# Totalling counts
for b in bfas:
    subp[b]['Total.Counts'] = np.sum(subp[b][[i for i in subp[b].columns if b in i]], axis=1)
# Getting rid of the bcs that have 0 counts in subpools but be in dbfa2_extra
subp['dBFA2'] = subp['dBFA2'].loc[subp['dBFA2']['Total.Counts'] > 0]

# Assigning environment calls based on subpool sequencing
for b in bfas:
    envs = set([i.split('-')[1] for i in subp[b].columns if b in i])
    subp[b]['Subpool.Environment'] = subp[b].apply(lambda row: call_env(row, envs), axis=1)
    subp[b]['Which.Subpools'] = subp[b].apply(lambda row: call_env_reps(row), axis=1)


# Adding to data files
for b in bfas:
    dats[b] = dats[b].merge(subp[b][['Full.BC', 'Subpool.Environment', 'Which.Subpools']], on='Full.BC', how='left').fillna('not_read')
    dats[b].to_csv('../../BFA_data/Combined_Counts/' + b + '_counts_with_env_info.csv', index=False)