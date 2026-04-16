#!/usr/bin/env python3
"""
plot_results.py — Space-time tradeoff plots for rank/select dictionaries.

Reads results/comparison.csv (DNA/5GRAM/URL + GOV2 averages) and optionally
results/sux_ef_comparison.csv, then produces PDF figures: bits-per-element (x)
vs. nanoseconds-per-query (y) for select and rank, grouped by dataset family.

Usage:
    python plot_results.py [--results-dir results/] [--output results/figures/]
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# ── Style palette ─────────────────────────────────────────────────────────────
# Each entry: color, marker, linestyle, linewidth, markersize, display label
STYLES = {
    'Array':             dict(color='#555555', marker='X',  ls='--',    lw=1.2, ms=7,  label='Array'),
    'EF (SDSL)':         dict(color='#1f77b4', marker='o',  ls='-',     lw=1.5, ms=7,  label='EF (SDSL)'),
    'EF (sux)':          dict(color='#17becf', marker='D',  ls='-',     lw=1.5, ms=6,  label='EF (sux)'),
    'RRR':               dict(color='#ff7f0e', marker='s',  ls='-',     lw=1.5, ms=6,  label='RRR'),
    'RLE':               dict(color='#2ca02c', marker='^',  ls='-',     lw=1.5, ms=6,  label='RLE'),
    'la_vector':         dict(color='#d62728', marker='*',  ls='-',     lw=2.2, ms=10, label='la_vector'),
    'la_vector_opt':     dict(color='#8c0000', marker='P',  ls='none',  lw=1.0, ms=10, label='la_vector (opt)'),
    'enc (delta)':       dict(color='#9467bd', marker='v',  ls='-',     lw=1.5, ms=6,  label='enc_vector (δ)'),
    'enc (gamma)':       dict(color='#c5b0d5', marker='^',  ls='-',     lw=1.5, ms=6,  label='enc_vector (γ)'),
    'hyb (uniform)':     dict(color='#8c564b', marker='p',  ls='-',     lw=1.5, ms=7,  label='hyb (uniform)'),
    'hyb (partitioned)': dict(color='#e377c2', marker='h',  ls='-',     lw=1.5, ms=7,  label='hyb (partitioned)'),
    's18':               dict(color='#bcbd22', marker='d',  ls='-',     lw=1.5, ms=6,  label='s18'),
}

# ── Column definitions ────────────────────────────────────────────────────────
# Each structure: list of (select_col, bpk_col, rank_col) per variant (block size / bpc).
# Variants within a group are connected by a line (space-time curve).
STRUCTURES = {
    'Array': [
        ('array_time_select', 'array_bpk', 'array_time_rank'),
    ],
    'EF (SDSL)': [
        ('ef_sd_time_select', 'ef_sd_bpk', 'ef_sd_time_rank'),
    ],
    'RRR': [
        ('rrr_15_time_select',  'rrr_15_bpk',  'rrr_15_time_rank'),
        ('rrr_31_time_select',  'rrr_31_bpk',  'rrr_31_time_rank'),
        ('rrr_63_time_select',  'rrr_63_bpk',  'rrr_63_time_rank'),
        ('rrr_127_time_select', 'rrr_127_bpk', 'rrr_127_time_rank'),
    ],
    'RLE': [
        ('rle_32_time_select',  'rle_32_bpk',  'rle_32_time_rank'),
        ('rle_64_time_select',  'rle_64_bpk',  'rle_64_time_rank'),
        ('rle_96_time_select',  'rle_96_bpk',  'rle_96_time_rank'),
        ('rle_128_time_select', 'rle_128_bpk', 'rle_128_time_rank'),
        ('rle_160_time_select', 'rle_160_bpk', 'rle_160_time_rank'),
    ],
    'la_vector': [
        ('la_vector_6_time_select',  'la_vector_6_bpk',  'la_vector_6_time_rank'),
        ('la_vector_7_time_select',  'la_vector_7_bpk',  'la_vector_7_time_rank'),
        ('la_vector_8_time_select',  'la_vector_8_bpk',  'la_vector_8_time_rank'),
        ('la_vector_9_time_select',  'la_vector_9_bpk',  'la_vector_9_time_rank'),
        ('la_vector_10_time_select', 'la_vector_10_bpk', 'la_vector_10_time_rank'),
        ('la_vector_11_time_select', 'la_vector_11_bpk', 'la_vector_11_time_rank'),
        ('la_vector_12_time_select', 'la_vector_12_bpk', 'la_vector_12_time_rank'),
        ('la_vector_13_time_select', 'la_vector_13_bpk', 'la_vector_13_time_rank'),
        ('la_vector_14_time_select', 'la_vector_14_bpk', 'la_vector_14_time_rank'),
    ],
    'la_vector_opt': [
        ('la_vector_opt_time_select', 'la_vector_opt_bpk', 'la_vector_opt_time_rank'),
    ],
    'enc (delta)': [
        ('elias_delta_4_time_select', 'elias_delta_4_bpk', 'elias_delta_4_time_rank'),
        ('elias_delta_5_time_select', 'elias_delta_5_bpk', 'elias_delta_5_time_rank'),
        ('elias_delta_6_time_select', 'elias_delta_6_bpk', 'elias_delta_6_time_rank'),
        ('elias_delta_7_time_select', 'elias_delta_7_bpk', 'elias_delta_7_time_rank'),
    ],
    'enc (gamma)': [
        ('elias_gamma_4_time_select', 'elias_gamma_4_bpk', 'elias_gamma_4_time_rank'),
        ('elias_gamma_5_time_select', 'elias_gamma_5_bpk', 'elias_gamma_5_time_rank'),
        ('elias_gamma_6_time_select', 'elias_gamma_6_bpk', 'elias_gamma_6_time_rank'),
        ('elias_gamma_7_time_select', 'elias_gamma_7_bpk', 'elias_gamma_7_time_rank'),
    ],
    'hyb (uniform)': [
        ('uniform_hyb0_select', 'uniform_hyb0_sequence_bpk', 'uniform_hyb0_rank'),
        ('uniform_hyb1_select', 'uniform_hyb1_sequence_bpk', 'uniform_hyb1_rank'),
    ],
    'hyb (partitioned)': [
        ('partitioned_hyb0_select', 'partitioned_hyb0_bpk', 'partitioned_hyb0_rank'),
        ('partitioned_hyb1_select', 'partitioned_hyb1_bpk', 'partitioned_hyb1_rank'),
        ('partitioned_hyb2_select', 'partitioned_hyb2_bpk', 'partitioned_hyb2_rank'),
    ],
    's18': [
        ('s18_1_select', 's18_1_bpk', 's18_1_rank'),
        ('s18_2_select', 's18_2_bpk', 's18_2_rank'),
        ('s18_3_select', 's18_3_bpk', 's18_3_rank'),
        ('s18_4_select', 's18_4_bpk', 's18_4_rank'),
        ('s18_5_select', 's18_5_bpk', 's18_5_rank'),
    ],
}

# sux-bench adds a separate row for each dataset; we attach it via a join.
SUX_STRUCT = 'EF (sux)'
SUX_COLS = ('sux_ef_time_select', 'sux_ef_bpk', 'sux_ef_time_rank')

# ── Dataset groups ────────────────────────────────────────────────────────────
# (display_name, [list of filename strings as stored in the CSV, without quotes])
DATASET_GROUPS = [
    ('DNA',             ['DNA_1', 'DNA_2', 'DNA_3']),
    ('5-GRAM',          ['5GRAM_1', '5GRAM_2', '5GRAM_3']),
    ('URL',             ['URL_1', 'URL_2', 'URL_3']),
    ('GOV2 (10M+)',     ['GOV2_AVG_10M-']),
    ('GOV2 (1M–10M)',   ['GOV2_AVG_1M-10M']),
    ('GOV2 (100K–1M)',  ['GOV2_AVG_100K-1M']),
]


# ── I/O helpers ───────────────────────────────────────────────────────────────

def load_csv(path: str) -> pd.DataFrame | None:
    """Load a benchmark CSV, stripping the surrounding quotes from filenames."""
    if not os.path.isfile(path):
        return None
    df = pd.read_csv(path)
    # Strip surrounding single-quotes from the filename column
    df['filename'] = df['filename'].str.strip("'")
    # Coerce all non-metadata columns to float (some may be read as object)
    meta = {'filename', 'n', 'u', 'ratio', 'mean_gap', 'stddev_gap',
            'mean_run', 'stddev_run'}
    for col in df.columns:
        if col not in meta:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    return df


def merge_sux(df: pd.DataFrame, sux_path: str) -> pd.DataFrame:
    """Left-join sux-bench columns onto the main dataframe."""
    sux = load_csv(sux_path)
    if sux is None:
        return df
    cols = ['filename'] + [c for c in sux.columns
                           if c.startswith('sux_ef')]
    return df.merge(sux[cols], on='filename', how='left')


def group_avg(df: pd.DataFrame, filenames: list) -> pd.Series | None:
    """Return the row mean across the given filenames (numeric cols only)."""
    rows = df[df['filename'].isin(filenames)]
    if rows.empty:
        return None
    numeric = rows.select_dtypes(include='number')
    return numeric.mean()


# ── Plotting helpers ──────────────────────────────────────────────────────────

def _get_xy(row: pd.Series, tuples: list, op: str) -> tuple[list, list]:
    """
    Extract (bpk, time) pairs for one structure/op from a data row.
    op is 'select' (index 0 of tuple) or 'rank' (index 2).
    Returns (xs, ys) lists, skipping any NaN entries.
    """
    t_idx = 0 if op == 'select' else 2
    xs, ys = [], []
    for tup in tuples:
        bpk_col = tup[1]
        t_col   = tup[t_idx]
        if bpk_col in row.index and t_col in row.index:
            bpk = row[bpk_col]
            t   = row[t_col]
            if pd.notna(bpk) and pd.notna(t) and bpk > 0 and t > 0:
                xs.append(bpk)
                ys.append(t)
    return xs, ys


def plot_tradeoff(ax, row: pd.Series, op: str, structures: dict,
                  sux_present: bool, draw_legend: bool = False):
    """
    Draw a space-time tradeoff scatter/curve on *ax* for one dataset row.

    op:  'select' or 'rank'
    """
    handles = []

    # --- main structures from STRUCTURES dict ---
    for name, tuples in structures.items():
        if name == 'EF (sux)':
            continue  # handled separately below
        st = STYLES[name]
        xs, ys = _get_xy(row, tuples, op)
        if not xs:
            continue
        # Sort by bpk so the connecting line makes sense
        pairs = sorted(zip(xs, ys))
        xs = [p[0] for p in pairs]
        ys = [p[1] for p in pairs]

        line, = ax.plot(xs, ys,
                        color=st['color'],
                        marker=st['marker'],
                        ls=st['ls'],
                        lw=st['lw'],
                        markersize=st['ms'],
                        zorder=3 if name == 'la_vector' else 2)
        handles.append(mlines.Line2D([], [],
                                     color=st['color'],
                                     marker=st['marker'],
                                     ls=st['ls'],
                                     lw=st['lw'],
                                     markersize=st['ms'],
                                     label=st['label']))

    # --- sux EF ---
    if sux_present and SUX_COLS[1] in row.index:
        bpk = row[SUX_COLS[1]]
        t_col = SUX_COLS[0] if op == 'select' else SUX_COLS[2]
        t = row.get(t_col, float('nan'))
        if pd.notna(bpk) and pd.notna(t) and bpk > 0 and t > 0:
            st = STYLES[SUX_STRUCT]
            ax.plot([bpk], [t],
                    color=st['color'],
                    marker=st['marker'],
                    ls='none',
                    markersize=st['ms'],
                    zorder=3)
            handles.append(mlines.Line2D([], [],
                                         color=st['color'],
                                         marker=st['marker'],
                                         ls='none',
                                         markersize=st['ms'],
                                         label=st['label']))

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, which='both', ls=':', lw=0.5, alpha=0.6)
    ax.set_xlabel('bits per element', fontsize=9)
    ax.set_ylabel('ns / query', fontsize=9)
    ax.tick_params(labelsize=8)

    if draw_legend:
        ax.legend(handles=handles, fontsize=7, ncol=2,
                  loc='upper left', framealpha=0.85)

    return handles


# ── Main figure generation ────────────────────────────────────────────────────

def make_figure(df: pd.DataFrame, groups: list, sux_present: bool,
                title: str, out_path: str):
    """
    Create one figure: rows = dataset groups, cols = (select, rank).
    """
    n_groups = len(groups)
    fig, axes = plt.subplots(n_groups, 2,
                             figsize=(10, 3.2 * n_groups),
                             squeeze=False)
    fig.suptitle(title, fontsize=11, y=1.002)

    all_handles = None
    for row_idx, (gname, filenames) in enumerate(groups):
        avg = group_avg(df, filenames)
        if avg is None:
            for col_idx in range(2):
                axes[row_idx][col_idx].set_visible(False)
            continue

        for col_idx, op in enumerate(('select', 'rank')):
            ax = axes[row_idx][col_idx]
            draw_leg = (row_idx == 0 and col_idx == 1)
            handles = plot_tradeoff(ax, avg, op, STRUCTURES,
                                    sux_present, draw_legend=draw_leg)
            if draw_leg:
                all_handles = handles
            ax.set_title(f'{gname} — {op}', fontsize=9)

    fig.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f'  wrote {out_path}')


def make_per_dataset_figures(df: pd.DataFrame, groups: list, sux_present: bool,
                             out_dir: str):
    """
    One figure per dataset group, showing select (left) and rank (right).
    """
    os.makedirs(out_dir, exist_ok=True)
    for gname, filenames in groups:
        avg = group_avg(df, filenames)
        if avg is None:
            print(f'  skip {gname}: no data')
            continue

        fig, (ax_sel, ax_rnk) = plt.subplots(1, 2, figsize=(11, 4.5))
        fig.suptitle(gname, fontsize=11)

        handles_sel = plot_tradeoff(ax_sel, avg, 'select', STRUCTURES,
                                    sux_present)
        handles_rnk = plot_tradeoff(ax_rnk, avg, 'rank',   STRUCTURES,
                                    sux_present)
        ax_sel.set_title('Select', fontsize=10)
        ax_rnk.set_title('Rank',   fontsize=10)

        # Shared legend below the subplots
        handles = handles_sel  # same structures in both
        fig.legend(handles=handles,
                   loc='lower center',
                   ncol=4,
                   fontsize=8,
                   bbox_to_anchor=(0.5, -0.12),
                   framealpha=0.9)
        fig.tight_layout()
        safe = gname.replace(' ', '_').replace('–', '-').replace('+', 'plus')
        out_path = os.path.join(out_dir, f'{safe}.pdf')
        fig.savefig(out_path, bbox_inches='tight', dpi=150)
        plt.close(fig)
        print(f'  wrote {out_path}')


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results-dir', default='results',
                    help='directory containing the CSV files (default: results/)')
    ap.add_argument('--output', default='results/figures',
                    help='output directory for PDF figures (default: results/figures/)')
    ap.add_argument('--combined', action='store_true',
                    help='also write a single combined figure (all groups in one PDF)')
    args = ap.parse_args()

    csv_path = os.path.join(args.results_dir, 'comparison.csv')
    sux_path = os.path.join(args.results_dir, 'sux_ef_comparison.csv')

    df = load_csv(csv_path)
    if df is None:
        print(f'ERROR: {csv_path} not found. Run benchmarks first.', file=sys.stderr)
        sys.exit(1)

    sux_present = os.path.isfile(sux_path)
    if sux_present:
        df = merge_sux(df, sux_path)
        print(f'Loaded sux-bench results from {sux_path}')
    else:
        print(f'Note: {sux_path} not found — EF (sux) will be omitted.')

    print(f'Loaded {len(df)} dataset rows: {list(df["filename"])}')

    # Filter groups to those actually present in the data
    present = set(df['filename'])
    active_groups = []
    for gname, fnames in DATASET_GROUPS:
        if any(f in present for f in fnames):
            active_groups.append((gname, fnames))
        else:
            print(f'  skip group "{gname}": none of {fnames} found in CSV')

    if not active_groups:
        print('No matching dataset groups found. Exiting.', file=sys.stderr)
        sys.exit(1)

    print(f'\nGenerating per-dataset figures → {args.output}/')
    make_per_dataset_figures(df, active_groups, sux_present, args.output)

    if args.combined:
        out_pdf = os.path.join(args.output, 'all_datasets.pdf')
        print(f'\nGenerating combined figure → {out_pdf}')
        make_figure(df, active_groups, sux_present,
                    'Space-time tradeoffs — all datasets', out_pdf)

    print('\nDone.')


if __name__ == '__main__':
    main()
