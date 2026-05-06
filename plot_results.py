#!/usr/bin/env python3
"""
plot_results.py — Space-time tradeoff plots for rank/select dictionaries.

Reads results/comparison.csv (DNA/5GRAM/URL + GOV2 averages) and optionally
results/sux_ef_comparison.csv, then produces three PDF figures:
  - results/figures/select.pdf  (space vs select time per query)
  - results/figures/rank.pdf    (space vs rank time per query)
  - results/figures/build.pdf   (space vs build time per key, log-scale X)

Each figure is a 4×3 grid: rows = dataset families (GOV2, URL, 5GRAM, DNA),
columns = sizes sparse → dense within each family. Y = space (bpk), clipped
at 16 bpk. Select/rank X-axes are linear ns/query; build X-axis is log ns/key
(la_vector_opt dominates, so log scale is needed).

Usage:
    python plot_results.py [--results-dir results/] [--output-dir results/figures/]
"""

import argparse
import os
import sys

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# ── Style palette ─────────────────────────────────────────────────────────────
STYLES = {
    'Array':             dict(color='#555555', marker='X',  ls='--',   lw=1.2, ms=6,  label='Array'),
    'EF (SDSL)':         dict(color='#1f77b4', marker='o',  ls='-',    lw=1.5, ms=7,  zorder=4, label='EF (SDSL)'),
    'EF (sux)':          dict(color='#17becf', marker='D',  ls='-',    lw=1.5, ms=4,  zorder=2, label='EF (sux)'),
    'EF (sux C++ w=2)':  dict(color='#0e8a7a', marker='d',  ls='-',    lw=1.5, ms=5,  zorder=3, label='EF (sux C++ w=2)'),
    'EF (sux C++ w=3)':  dict(color='#6a0dad', marker='d',  ls='-',    lw=1.5, ms=5,  zorder=3, label='EF (sux C++ w=3)'),
    'RRR':               dict(color='#ff7f0e', marker='s',  ls='-',    lw=1.5, ms=5,  label='RRR'),
    'RLE':               dict(color='#2ca02c', marker='^',  ls='-',    lw=1.5, ms=5,  label='RLE'),
    'la_vector':         dict(color='#d62728', marker='*',  ls='-',    lw=2.0, ms=9,  zorder=3, label='la_vector'),
    'la_vector_opt':     dict(color='#8c0000', marker='P',  ls='none', lw=1.0, ms=9,  label='la_vector (opt)'),
    'enc (delta)':       dict(color='#9467bd', marker='v',  ls='-',    lw=1.5, ms=5,  label='enc_vector (δ)'),
    'enc (gamma)':       dict(color='#c5b0d5', marker='^',  ls='-',    lw=1.5, ms=5,  label='enc_vector (γ)'),
    'hyb (uniform)':     dict(color='#8c564b', marker='p',  ls='-',    lw=1.5, ms=6,  label='hyb (uniform)'),
    'hyb (partitioned)': dict(color='#e377c2', marker='h',  ls='-',    lw=1.5, ms=6,  label='hyb (partitioned)'),
    's18':               dict(color='#bcbd22', marker='d',  ls='-',    lw=1.5, ms=5,  label='s18'),
}

# ── Column definitions ────────────────────────────────────────────────────────
# Each structure: list of (select_col, bpk_col, rank_col, build_col) per variant.
# Variants are connected by a line (space-time tradeoff curve).
STRUCTURES = {
    'Array': [
        ('array_time_select', 'array_bpk', 'array_time_rank', 'array_time_build'),
    ],
    'EF (SDSL)': [
        ('ef_sd_time_select', 'ef_sd_bpk', 'ef_sd_time_rank', 'ef_sd_time_build'),
    ],
    'RRR': [
        ('rrr_15_time_select',  'rrr_15_bpk',  'rrr_15_time_rank',  'rrr_15_time_build'),
        ('rrr_31_time_select',  'rrr_31_bpk',  'rrr_31_time_rank',  'rrr_31_time_build'),
        ('rrr_63_time_select',  'rrr_63_bpk',  'rrr_63_time_rank',  'rrr_63_time_build'),
        ('rrr_127_time_select', 'rrr_127_bpk', 'rrr_127_time_rank', 'rrr_127_time_build'),
    ],
    'RLE': [
        ('rle_32_time_select',  'rle_32_bpk',  'rle_32_time_rank',  'rle_32_time_build'),
        ('rle_64_time_select',  'rle_64_bpk',  'rle_64_time_rank',  'rle_64_time_build'),
        ('rle_96_time_select',  'rle_96_bpk',  'rle_96_time_rank',  'rle_96_time_build'),
        ('rle_128_time_select', 'rle_128_bpk', 'rle_128_time_rank', 'rle_128_time_build'),
        ('rle_160_time_select', 'rle_160_bpk', 'rle_160_time_rank', 'rle_160_time_build'),
    ],
    'la_vector': [
        ('la_vector_6_time_select',  'la_vector_6_bpk',  'la_vector_6_time_rank',  'la_vector_6_time_build'),
        ('la_vector_7_time_select',  'la_vector_7_bpk',  'la_vector_7_time_rank',  'la_vector_7_time_build'),
        ('la_vector_8_time_select',  'la_vector_8_bpk',  'la_vector_8_time_rank',  'la_vector_8_time_build'),
        ('la_vector_9_time_select',  'la_vector_9_bpk',  'la_vector_9_time_rank',  'la_vector_9_time_build'),
        ('la_vector_10_time_select', 'la_vector_10_bpk', 'la_vector_10_time_rank', 'la_vector_10_time_build'),
        ('la_vector_11_time_select', 'la_vector_11_bpk', 'la_vector_11_time_rank', 'la_vector_11_time_build'),
        ('la_vector_12_time_select', 'la_vector_12_bpk', 'la_vector_12_time_rank', 'la_vector_12_time_build'),
        ('la_vector_13_time_select', 'la_vector_13_bpk', 'la_vector_13_time_rank', 'la_vector_13_time_build'),
        ('la_vector_14_time_select', 'la_vector_14_bpk', 'la_vector_14_time_rank', 'la_vector_14_time_build'),
    ],
    'la_vector_opt': [
        ('la_vector_opt_time_select', 'la_vector_opt_bpk', 'la_vector_opt_time_rank', 'la_vector_opt_time_build'),
    ],
    'enc (delta)': [
        ('elias_delta_4_time_select', 'elias_delta_4_bpk', 'elias_delta_4_time_rank', 'elias_delta_4_time_build'),
        ('elias_delta_5_time_select', 'elias_delta_5_bpk', 'elias_delta_5_time_rank', 'elias_delta_5_time_build'),
        ('elias_delta_6_time_select', 'elias_delta_6_bpk', 'elias_delta_6_time_rank', 'elias_delta_6_time_build'),
        ('elias_delta_7_time_select', 'elias_delta_7_bpk', 'elias_delta_7_time_rank', 'elias_delta_7_time_build'),
    ],
    'enc (gamma)': [
        ('elias_gamma_4_time_select', 'elias_gamma_4_bpk', 'elias_gamma_4_time_rank', 'elias_gamma_4_time_build'),
        ('elias_gamma_5_time_select', 'elias_gamma_5_bpk', 'elias_gamma_5_time_rank', 'elias_gamma_5_time_build'),
        ('elias_gamma_6_time_select', 'elias_gamma_6_bpk', 'elias_gamma_6_time_rank', 'elias_gamma_6_time_build'),
        ('elias_gamma_7_time_select', 'elias_gamma_7_bpk', 'elias_gamma_7_time_rank', 'elias_gamma_7_time_build'),
    ],
    'hyb (uniform)': [
        ('uniform_hyb0_select', 'uniform_hyb0_sequence_bpk', 'uniform_hyb0_rank', 'uniform_hyb0_time_build'),
        ('uniform_hyb1_select', 'uniform_hyb1_sequence_bpk', 'uniform_hyb1_rank', 'uniform_hyb1_time_build'),
    ],
    'hyb (partitioned)': [
        ('partitioned_hyb0_select', 'partitioned_hyb0_bpk', 'partitioned_hyb0_rank', 'partitioned_hyb0_time_build'),
        ('partitioned_hyb1_select', 'partitioned_hyb1_bpk', 'partitioned_hyb1_rank', 'partitioned_hyb1_time_build'),
        ('partitioned_hyb2_select', 'partitioned_hyb2_bpk', 'partitioned_hyb2_rank', 'partitioned_hyb2_time_build'),
    ],
    's18': [
        ('s18_1_select', 's18_1_bpk', 's18_1_rank', 's18_1_time_build'),
        ('s18_2_select', 's18_2_bpk', 's18_2_rank', 's18_2_time_build'),
        ('s18_3_select', 's18_3_bpk', 's18_3_rank', 's18_3_time_build'),
        ('s18_4_select', 's18_4_bpk', 's18_4_rank', 's18_4_time_build'),
        ('s18_5_select', 's18_5_bpk', 's18_5_rank', 's18_5_time_build'),
    ],
    'EF (sux)': [
        ('sux_ef_9_time_select',  'sux_ef_9_bpk',  'sux_ef_9_time_rank',  'sux_ef_9_time_build'),
        ('sux_ef_10_time_select', 'sux_ef_10_bpk', 'sux_ef_10_time_rank', 'sux_ef_10_time_build'),
        ('sux_ef_11_time_select', 'sux_ef_11_bpk', 'sux_ef_11_time_rank', 'sux_ef_11_time_build'),
        ('sux_ef_12_time_select', 'sux_ef_12_bpk', 'sux_ef_12_time_rank', 'sux_ef_12_time_build'),
        ('sux_ef_13_time_select', 'sux_ef_13_bpk', 'sux_ef_13_time_rank', 'sux_ef_13_time_build'),
    ],
    'EF (sux C++ w=2)': [
        ('sux_cpp_ef_9_2_time_select',  'sux_cpp_ef_9_2_bpk',  'sux_cpp_ef_9_2_time_rank',  'sux_cpp_ef_9_2_time_build'),
        ('sux_cpp_ef_10_2_time_select', 'sux_cpp_ef_10_2_bpk', 'sux_cpp_ef_10_2_time_rank', 'sux_cpp_ef_10_2_time_build'),
        ('sux_cpp_ef_11_2_time_select', 'sux_cpp_ef_11_2_bpk', 'sux_cpp_ef_11_2_time_rank', 'sux_cpp_ef_11_2_time_build'),
    ],
    'EF (sux C++ w=3)': [
        ('sux_cpp_ef_9_3_time_select',  'sux_cpp_ef_9_3_bpk',  'sux_cpp_ef_9_3_time_rank',  'sux_cpp_ef_9_3_time_build'),
        ('sux_cpp_ef_10_3_time_select', 'sux_cpp_ef_10_3_bpk', 'sux_cpp_ef_10_3_time_rank', 'sux_cpp_ef_10_3_time_build'),
        ('sux_cpp_ef_11_3_time_select', 'sux_cpp_ef_11_3_bpk', 'sux_cpp_ef_11_3_time_rank', 'sux_cpp_ef_11_3_time_build'),
    ],
}

BPK_MAX  = 16.0   # filter out structures using more than 16 bits/key (matches paper)
TIME_MAX = 200.0  # clip X-axis (select/rank) at 200 ns/query

# op → tuple index in STRUCTURES variants
OP_INDEX = {'select': 0, 'rank': 2, 'build': 3}

# ── Paper grid layout ─────────────────────────────────────────────────────────
# 4 rows × 3 cols: rows = dataset families (GOV2→DNA top→bottom),
# cols = sparse → dense within each family.
# Each cell: (base_title, [filenames_to_average])
PAPER_GRID = [
    # Row 1: GOV2
    [('GOV2 avg 100K–1M', ['GOV2_AVG_100K-1M']),
     ('GOV2 avg 1M–10M',  ['GOV2_AVG_1M-10M']),
     ('GOV2 avg +10M',    ['GOV2_AVG_10M-'])],
    # Row 2: URL
    [('URL', ['URL_3']), ('URL', ['URL_2']), ('URL', ['URL_1'])],
    # Row 3: 5GRAM
    [('5GRAM', ['5GRAM_3']), ('5GRAM', ['5GRAM_2']), ('5GRAM', ['5GRAM_1'])],
    # Row 4: DNA
    [('DNA', ['DNA_3']),  ('DNA', ['DNA_2']),  ('DNA', ['DNA_1'])],
]


# ── I/O helpers ───────────────────────────────────────────────────────────────

def load_csv(path: str) -> 'pd.DataFrame | None':
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


def merge_sux(df: 'pd.DataFrame', sux_path: str) -> 'pd.DataFrame':
    """Left-join sux-bench columns; backfill GOV2_AVG_* rows by dataset size."""
    sux = load_csv(sux_path)
    if sux is None:
        return df
    sux_cols = [c for c in sux.columns if c.startswith('sux_')]
    merged = df.merge(sux[['filename'] + sux_cols], on='filename', how='left')

    if 'n' in sux.columns and not sux.empty:
        gov2 = sux[sux['filename'].str.match(r'GOV2_\d+')]
        ranges = {
            'GOV2_AVG_10M-':    gov2[gov2['n'] >= 10_000_000],
            'GOV2_AVG_1M-10M':  gov2[(gov2['n'] >= 1_000_000) & (gov2['n'] < 10_000_000)],
            'GOV2_AVG_100K-1M': gov2[(gov2['n'] >= 100_000)   & (gov2['n'] < 1_000_000)],
        }
        for avg_name, subset in ranges.items():
            if subset.empty:
                continue
            avg_vals = subset[sux_cols].mean()
            mask = merged['filename'] == avg_name
            if mask.any():
                for col in sux_cols:
                    merged.loc[mask, col] = avg_vals[col]
    return merged


# ── Plotting helpers ──────────────────────────────────────────────────────────

def _get_time_bpk(row: 'pd.Series', tuples: list, op: str,
                  bpk_max: float = BPK_MAX) -> 'tuple[list, list]':
    """
    Return (time_values, bpk_values) for one structure across its variants.
    Points with bpk > bpk_max are excluded.
    op: 'select', 'rank', or 'build'.
    """
    t_idx = OP_INDEX[op]
    ts, bs = [], []
    for tup in tuples:
        if t_idx >= len(tup):
            continue
        bpk_col, t_col = tup[1], tup[t_idx]
        if bpk_col in row.index and t_col in row.index:
            bpk = row[bpk_col]
            t   = row[t_col]
            if pd.notna(bpk) and pd.notna(t) and bpk > 0 and t > 0:
                if bpk <= bpk_max:
                    ts.append(t)
                    bs.append(bpk)
    return ts, bs


def plot_cell(ax, row: 'pd.Series', op: str) -> list:
    """
    Draw one subplot cell (one dataset × one operation).
    Y = space (bpk), clipped at BPK_MAX.
    X = time per query (select/rank, linear) or per key (build, log).
    Returns legend handles.
    """
    handles = []
    log_x = (op == 'build')

    for name, tuples in STRUCTURES.items():
        st = STYLES[name]
        ts, bs = _get_time_bpk(row, tuples, op)
        if not ts:
            continue
        pairs = sorted(zip(ts, bs))
        ts_s = [p[0] for p in pairs]
        bs_s = [p[1] for p in pairs]
        ax.plot(ts_s, bs_s,
                color=st['color'], marker=st['marker'],
                ls=st['ls'], lw=st['lw'], markersize=st['ms'],
                zorder=st.get('zorder', 2))
        handles.append(mlines.Line2D([], [],
                                     color=st['color'], marker=st['marker'],
                                     ls=st['ls'], lw=st['lw'], markersize=st['ms'],
                                     label=st['label']))

    if log_x:
        ax.set_xscale('log')
    else:
        ax.set_xlim(left=0, right=TIME_MAX)
    ax.set_ylim(top=BPK_MAX)
    ax.grid(True, ls=':', lw=0.4, alpha=0.6, which='both' if log_x else 'major')
    ax.tick_params(labelsize=7)
    return handles


# ── Paper-style figure generation ────────────────────────────────────────────

def make_paper_figure(df: 'pd.DataFrame', op: str,
                      out_path: str) -> None:
    """
    4×3 grid figure matching Figs 7 / 8 of the paper.
    op: 'select' or 'rank'
    """
    present = set(df['filename'])
    x_label = {
        'select': 'Select time (nanoseconds per query)',
        'rank':   'Rank time (nanoseconds per query)',
        'build':  'Build time (nanoseconds per key, log scale)',
    }[op]

    title = {
        'select': 'select queries',
        'rank':   'rank queries',
        'build':  'construction',
    }[op]

    fig, axes = plt.subplots(4, 3, figsize=(11, 13), squeeze=False)
    fig.suptitle(f'Space-time performance of {title}', fontsize=10, y=1.002)

    legend_handles = None

    for row_idx, row_specs in enumerate(PAPER_GRID):
        for col_idx, (base_title, filenames) in enumerate(row_specs):
            ax = axes[row_idx][col_idx]

            avail = [f for f in filenames if f in present]
            if not avail:
                ax.set_visible(False)
                continue

            rows = df[df['filename'].isin(avail)]
            avg  = rows.select_dtypes(include='number').mean()

            # Density label from ratio column
            density = ''
            if 'ratio' in avg.index and pd.notna(avg['ratio']):
                density = f' ({avg["ratio"] * 100:.1f}%)'

            handles = plot_cell(ax, avg, op)
            if legend_handles is None and handles:
                legend_handles = handles

            ax.set_title(f'{base_title}{density}', fontsize=8)

            # Axis labels only on outer edges
            if col_idx == 0:
                ax.set_ylabel('Space (bits per integer)', fontsize=7)
            if row_idx == 3:
                ax.set_xlabel(x_label, fontsize=7)

    fig.tight_layout(rect=[0, 0.06, 1, 1])

    if legend_handles:
        fig.legend(handles=legend_handles,
                   loc='lower center',
                   ncol=4,
                   fontsize=6.5,
                   bbox_to_anchor=(0.5, 0.0),
                   framealpha=0.9)

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f'  wrote {out_path}')


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results-dir', default='results',
                    help='directory containing the CSV files (default: results/)')
    ap.add_argument('--output-dir', default=None,
                    help='output directory for PDF figures (default: <results-dir>/figures/)')
    args = ap.parse_args()

    output_dir = output_dir_dir if output_dir_dir else os.path.join(args.results_dir, 'figures')

    csv_path = os.path.join(args.results_dir, 'comparison.csv')
    sux_path = os.path.join(args.results_dir, 'sux_ef_comparison.csv')
    sux_cpp_path = os.path.join(args.results_dir, 'sux_cpp_ef_comparison.csv')

    df = load_csv(csv_path)
    if df is None:
        print(f'ERROR: {csv_path} not found. Run benchmarks first.', file=sys.stderr)
        sys.exit(1)

    if os.path.isfile(sux_path):
        df = merge_sux(df, sux_path)
        print(f'Loaded sux-bench results from {sux_path}')
    else:
        print(f'Note: {sux_path} not found — EF (sux) will be omitted.')

    if os.path.isfile(sux_cpp_path):
        df = merge_sux(df, sux_cpp_path)
        print(f'Loaded sux C++ results from {sux_cpp_path}')
    else:
        print(f'Note: {sux_cpp_path} not found — EF (sux C++) will be omitted.')

    print(f'Loaded {len(df)} dataset rows: {list(df["filename"])}')

    os.makedirs(output_dir, exist_ok=True)
    for op in ('select', 'rank', 'build'):
        out = os.path.join(output_dir, f'{op}.pdf')
        print(f'\nGenerating {op}.pdf …')
        make_paper_figure(df, op, out)

    print('\nDone.')


if __name__ == '__main__':
    main()
