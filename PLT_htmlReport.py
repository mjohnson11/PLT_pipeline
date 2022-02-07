"""
A program to demultiplex reads and count barcodes for the PLT
Milo Johnson
Started 3_6_16
made simple 2_6_17
"""

import pandas as pd
import matplotlib
from io import StringIO
import os
import subprocess

matplotlib.use('agg')  # to make matplotlib work on a remote machine
import matplotlib.pyplot as pl

def add_qual_perc(row):
    return row['Quality.Failed'] / row['Total.Reads']


def add_reg_perc(row):
    return row['Regex.Failed'] / (row['Total.Reads']-row['Quality.Failed'])


def add_umi_perc(row):
    return row['UMI.Repeats'] / (row['Total.Reads']-row['Quality.Failed']-row['Regex.Failed'])


def add_usable_perc(row):
    return row['Usable.Reads'] / row['Total.Reads']


def make_4_subplots(tmp_rows, out_dir, graph_counter):
    fig, subps = pl.subplots(2, 2, figsize=(12,9))
    d = 0
    for subarr in subps:
        for sub in subarr:
            if len(tmp_rows) > d:
                sub.semilogy([i for i in range(1, len(tmp_rows[d]))], tmp_rows[d][1:], color='black')
                if d > 1:
                    sub.set_xlabel('UMI family size')
                if d % 2 == 0:
                    sub.set_ylabel('# of families')
                sub.set_title(tmp_rows[d][0])
                d += 1

    img_out = out_dir + '/graphs/UMI_graphs_' + str(graph_counter) + '.png'
    fig.savefig(img_out)
    pl.close(fig)


def make_html_report(stats_file, umi_fam_size_file, output_path, run_name):
    stats = pd.read_csv(stats_file)
    umi_sizes = pd.read_csv(umi_fam_size_file)

    stats['Percent.Quality.Failed'] = stats.apply(lambda row: add_qual_perc(row), axis=1)
    stats['Percent.Regex.Failed'] = stats.apply(lambda row: add_reg_perc(row), axis=1)
    stats['Percent.UMI.Repeats'] = stats.apply(lambda row: add_umi_perc(row), axis=1)
    stats['Percent.Usable.Reads'] = stats.apply(lambda row: add_usable_perc(row), axis=1)

    out_dir = output_path
    if not os.path.isdir(out_dir):
        print('making output directory:', out_dir)
        subprocess.call(['mkdir', out_dir])
        subprocess.call(['mkdir', out_dir + '/graphs'])
    elif not os.path.isdir(out_dir + '/graphs'):
        subprocess.call(['mkdir', out_dir + '/graphs'])
    tmp_rows = []
    graph_counter = 1
    for row in umi_sizes.as_matrix():
        tmp_rows.append(row)
        if len(tmp_rows) == 4:
            make_4_subplots(tmp_rows, out_dir, graph_counter)
            graph_counter += 1
            tmp_rows = []
    # final leftover graphs case:
    if len(tmp_rows) > 0:
        make_4_subplots(tmp_rows, out_dir, graph_counter)

    favorite_columns = ['Library', 'Total.Reads', 'Usable.Reads', 'Percent.Quality.Failed', 'Percent.Regex.Failed',
                        'Percent.UMI.Repeats', 'Percent.Usable.Reads']

    out_css_filename = out_dir + '/library_report.css'
    out_filename = out_dir + '/index.html'
    buf = StringIO()
    buf.write('<head>\n')
    buf.write('<title>'+run_name+'</title>\n')
    buf.write("<link rel='stylesheet' href='library_report.css' />\n</head>\n")
    buf.write('<div id="content_container">\n')
    buf.write('<h1>'+run_name+' Library Statistics</h1>\n')
    # stat_dats[r].to_html(buf, index=False)
    tmp_rows = []
    graph_counter = 1
    for row in umi_sizes.as_matrix():
        tmp_rows.append(row)
        if len(tmp_rows) == 4:
            libs = set([r[0] for r in tmp_rows])
            stats.loc[stats['Library'].isin(libs)][favorite_columns].to_html(buf, index=False)
            img_out = 'graphs/UMI_graphs_' + str(graph_counter) + '.png'
            buf.write('<img class="umi_graph" src=' + img_out + ' ></img>\n')
            graph_counter += 1
            tmp_rows = []
    # final leftover graphs case:
    if len(tmp_rows) > 0:
        libs = set([r[0] for r in tmp_rows])
        stats.loc[stats['Library'].isin(libs)][favorite_columns].to_html(buf, index=False)
        img_out = 'graphs/UMI_graphs_' + str(graph_counter) + '.png'
        buf.write('\n<img class="umi_graph" src=' + img_out + ' ></img>\n')

    buf.write('</div>\n')
    with open(out_filename, 'w') as f:
        f.write(buf.getvalue())
    with open(out_css_filename, 'w') as out_css:
        out_css.write('@charset "UTF-8";\n\n')
        out_css.write('#content_container {\ntext-align:center;\nwidth:1200px;\n')
        out_css.write('position: relative;\nmargin: 0px auto;\n}\n\n')
        out_css.write('.umi_graph {\nposition: relative;\nmargin: 0px auto;\n')
        out_css.write('width: 600px;\n}\n\n')
        out_css.write('table {\nposition: relative;\nmargin: 0px auto;\nborder-collapse: collapse;\n}\n\n')
        out_css.write('th, td {\ntext-align: left;\npadding: 8px;\n}\n\n')
        out_css.write('tr:nth-child(even){background-color: #f2f2f2}\n\n')