import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input_files', help='comma separated bc count files to combine')
parser.add_argument('output_file', help='output csv')
parser.add_argument('chimera_output', help='chimera output csv')
parser.add_argument("-merge_type", type=str, default='inner', help='include only shared bcs (default) or all ("outer")')
args = parser.parse_args()

input_files = args.input_files.split(',')
output_file = args.output_file
chimera_output = args.chimera_output
merge_type = args.merge_type

print('Processing files:', '\n'.join(input_files))

dats = []
for f in input_files:
    td = pd.read_csv(f)
    td['Full.BC'] = td['Diverse.BC'] + td['Environment.BC']
    dats.append(td)

bc_cols = ['Full.BC', 'Diverse.BC', 'Environment.BC']
time_cols = [c for c in dats[0].columns if 'Time' in c]
running_dat = dats[0][bc_cols + time_cols]
for d in dats[1:]:
    tmp_time_cols = [c for c in d.columns if 'Time' in c]
    repeat_time_cols = [c for c in tmp_time_cols if c in time_cols]
    new_time_cols = [c for c in tmp_time_cols if c not in time_cols]
    running_dat = running_dat.merge(d, on=bc_cols, how=merge_type, suffixes=('_old1', '_new1')).fillna(0)
    for t in repeat_time_cols:
        running_dat[t] = running_dat[t + '_old1'] + running_dat[t + '_new1']
    time_cols += new_time_cols
    running_dat = running_dat[bc_cols + time_cols]

running_dat['Total.Counts'] = running_dat.apply(lambda drow: sum([drow[i] for i in time_cols]), axis=1)
use_cols = bc_cols + ['Total.Counts'] + sorted(time_cols,
                                               key=lambda c: c[:c.index('Time')+4] + c[c.index('Time')+4:].zfill(3))

running_dat = running_dat.sort('Total.Counts', ascending=False)
# excluding chimeras
dd = dict(running_dat['Diverse.BC'].value_counts())
dbc_with_mult_hits = set([d for d in running_dat['Diverse.BC'] if dd[d] > 1])
fullbc_chimeras = []
dbc_top_counts = dict()
for entry in running_dat.as_matrix(['Full.BC', 'Diverse.BC', 'Total.Counts']):
    if entry[1] in dbc_with_mult_hits:
        if entry[1] in dbc_top_counts:
            if (entry[2] / dbc_top_counts[entry[1]]) < 0.01:
                fullbc_chimeras.append(entry[0])
        else:
            dbc_top_counts[entry[1]] = entry[2]
running_dat['is_chimeric'] = running_dat.apply(lambda drow: drow['Full.BC'] in fullbc_chimeras, axis=1)
d_no_chimeras = running_dat.loc[~running_dat['is_chimeric']]
d_yes_chimeras = running_dat.loc[running_dat['is_chimeric']]

d_yes_chimeras[use_cols].to_csv(chimera_output, index=False)

print(len(d_no_chimeras), 'bcs after combining and excluding chimeras.')

d_no_chimeras[use_cols].to_csv(output_file, index=False)
