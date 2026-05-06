#!/usr/bin/env python3
"""
gen_table.py — Generate a Markdown table of the top-k fastest select variants per dataset.

Reads the same CSV files as plot_results.py and outputs a Markdown table
split into two halves: Elias-Fano (sux Rust/C++) and la_vector.

Usage:
    python gen_table.py [--results-dir results/] [--top 4]
"""

import argparse
import os
import sys

import pandas as pd

GROUPS = {
    'EF (C++)': [
        ('sux_cpp_ef_9_3_time_select',    'sux_cpp_ef_9_3_bpk',    'C++(9)'),
        ('sux_cpp_ef_10_3_time_select',   'sux_cpp_ef_10_3_bpk',   'C++(10)'),
        ('sux_cpp_ef_11_3_time_select',   'sux_cpp_ef_11_3_bpk',   'C++(11)'),
    ],
    'EF (Rust)': [
        ('sux_ef_9_time_select',          'sux_ef_9_bpk',          'Rust(9)'),
        ('sux_ef_10_time_select',         'sux_ef_10_bpk',         'Rust(10)'),
        ('sux_ef_11_time_select',         'sux_ef_11_bpk',         'Rust(11)'),
        ('sux_ef_12_time_select',         'sux_ef_12_bpk',         'Rust(12)'),
        ('sux_ef_13_time_select',         'sux_ef_13_bpk',         'Rust(13)'),
    ],
    'LA': [
        ('la_vector_6_time_select',       'la_vector_6_bpk',       'la_vector<6>'),
        ('la_vector_7_time_select',       'la_vector_7_bpk',       'la_vector<7>'),
        ('la_vector_8_time_select',       'la_vector_8_bpk',       'la_vector<8>'),
        ('la_vector_9_time_select',       'la_vector_9_bpk',       'la_vector<9>'),
        ('la_vector_10_time_select',      'la_vector_10_bpk',      'la_vector<10>'),
        ('la_vector_11_time_select',      'la_vector_11_bpk',      'la_vector<11>'),
        ('la_vector_12_time_select',      'la_vector_12_bpk',      'la_vector<12>'),
        ('la_vector_13_time_select',      'la_vector_13_bpk',      'la_vector<13>'),
        ('la_vector_14_time_select',      'la_vector_14_bpk',      'la_vector<14>'),
    ],
}

DATASET_DISPLAY = {
    'GOV2_AVG_100K-1M': 'GOV2 avg 100K–1M',
    'GOV2_AVG_1M-10M':  'GOV2 avg 1M–10M',
    'GOV2_AVG_10M-':     'GOV2 avg 10M+',
}


def load_csv(path):
    if not os.path.isfile(path):
        return None
    df = pd.read_csv(path)
    df['filename'] = df['filename'].str.strip("'")
    meta = {'filename', 'n', 'u', 'ratio', 'mean_gap', 'stddev_gap',
            'mean_run', 'stddev_run'}
    for col in df.columns:
        if col not in meta:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    return df


def merge_sux(df, sux_path):
    sux = load_csv(sux_path)
    if sux is None:
        return df
    sux_cols = [c for c in sux.columns if c.startswith('sux_')]
    return df.merge(sux[['filename'] + sux_cols], on='filename', how='left')


def top_k(row, variants, k):
    results = []
    for sel_col, bpk_col, name in variants:
        if sel_col not in row.index or bpk_col not in row.index:
            continue
        t = row[sel_col]
        bpk = row[bpk_col]
        if pd.notna(t) and pd.notna(bpk) and t > 0 and bpk > 0:
            results.append((t, bpk, name))
    results.sort()
    return results[:k]


def fmt_entry(t, bpk, name):
    return f'{name} {t:.1f} ns ({bpk:.2f})'


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results-dir', default='results',
                    help='directory containing the CSV files (default: results/)')
    ap.add_argument('--top', type=int, default=4,
                    help='number of fastest selects per group (default: 4)')
    args = ap.parse_args()

    csv_path = os.path.join(args.results_dir, 'comparison.csv')
    sux_path = os.path.join(args.results_dir, 'sux_ef_comparison.csv')
    sux_cpp_path = os.path.join(args.results_dir, 'sux_cpp_ef_comparison.csv')

    df = load_csv(csv_path)
    if df is None:
        print(f'ERROR: {csv_path} not found.', file=sys.stderr)
        sys.exit(1)

    if os.path.isfile(sux_path):
        df = merge_sux(df, sux_path)
    if os.path.isfile(sux_cpp_path):
        df = merge_sux(df, sux_cpp_path)

    k = args.top
    group_names = list(GROUPS.keys())

    hdrs = []
    for gname in group_names:
        for i in range(k):
            hdrs.append(f'{gname} #{i+1}')
    header = '| Dataset | ' + ' | '.join(hdrs) + ' |'
    sep = '|---' * (1 + len(hdrs)) + '|'

    print(header)
    print(sep)

    for _, row in df.iterrows():
        name = row['filename']
        display = DATASET_DISPLAY.get(name, name)

        cells = []
        for gname in group_names:
            tops = top_k(row, GROUPS[gname], k)
            gcells = [fmt_entry(*e) for e in tops]
            while len(gcells) < k:
                gcells.append('')
            cells.extend(gcells)

        print(f'| {display} | ' + ' | '.join(cells) + ' |')


if __name__ == '__main__':
    main()
