# Learned-Compressed-Rank-Select-TALG22, Revisited

This repository is a fork of the [original Learned-Compressed-Rank-Select-TALG22
repository](https://github.com/aboffa/Learned-Compressed-Rank-Select-TALG22) for
the paper “[A Learned Approach to Design Compressed Rank/Select Data
Structures](https://dl.acm.org/doi/pdf/10.1145/3524060)”. It contains updated
code, including the C++ implementation (2007) and a new Rust port (2026) of the
Elias–Fano representation of monotone sequences from my paper “[Broadword
Implementation of Rank/Select
Queries](https://link.springer.com/chapter/10.1007/978-3-540-68552-4_12)”.

If you're looking for advice on the choice of a rank/select structure, the
results in the paper are misleading: the main claim of the authors is to have
the “fastest select” (excluding, of course, arrays), but they chose a poor
baseline for the state of the art. The fastest implementation of selection on
bits (and thus, cascading, of an Elias–Fano representation) has been since 2007
the interleaved, byte-aligned two-level inventory described in “Broadword
Implementation of Rank/Select Queries”. Subsequent research has always followed
the same guidelines (up to embedding the selection structure in the bit vector),
trying to improve space usage. The main reference implementation the authors
used (from [SDSL-lite](https://github.com/simongog/sdsl-lite)) violates both
byte alignment and interleaving, and it is indeed very slow.

This repository uses the machinery of the original repository to perform the
same benchmarks, but adding to the mix the 2007 implementation from the
[`sux`](https://github.com/vigna/sux) library used in “Broadword Implementation
of Rank/Select Queries”, and a recent Rust port of the same structure. The
picture is quite different from the one in the paper: the “learned” approach is
not competitive with the state of the art; it is actually slower than the 2007
implementation on select, much slower on rank, and still using in general more
space than the information-theoretical lower bound.

You can see here the results for [select](results/figures/select.pdf),
[rank](results/figures/rank.pdf), and [construction
time](results/figures/build.pdf) (the latter are in logarithmic scale because of
the very large build time of the `la_vector<opt>` variant; timings have been cut
at 200 ns to make the interesting part more understandable). The learned version
is always slower than a state-of-the-art Elias–Fano implementation, in some
cases almost twice as slow, and often uses more space. In fact, apart from a few
data points on the Pareto frontier, the learned version is dominated by
Elias–Fano both in time and space.

The only data point sometimes on the Pareto frontier is `la_vector<opt>`, due
to very good compression. This is, however, more of a proof-of-concept data
structure, as it is optimized exhaustively, and its construction time is 120-150
times slower than an Elias–Fano representation.

The following table shows the fastest variant for the 2007 C++
implementation, the Rust implementation, and the “learned” approach; the number
in parentheses is the space usage in bits per element, and the time is in
nanoseconds per query.

| Dataset          | EF (C++) #1            | EF (Rust) #1            | LA #1                         |
| ---------------- | ---------------------- | ----------------------- | ----------------------------- |
| DNA_1            | C++(9) 47.3 ns (6.66)  | Rust(9) 46.2 ns (6.66)  | la_vector<10> 48.5 ns (10.02) |
| DNA_2            | C++(11) 32.5 ns (6.62) | Rust(11) 31.1 ns (6.62) | la_vector<10> 42.0 ns (10.12) |
| DNA_3            | C++(9) 14.2 ns (10.89) | Rust(9) 13.8 ns (10.89) | la_vector<11> 21.9 ns (11.35) |
| 5GRAM_1          | C++(11) 39.3 ns (5.91) | Rust(11) 37.6 ns (5.93) | la_vector<12> 47.6 ns (12.05) |
| 5GRAM_2          | C++(9) 22.0 ns (10.48) | Rust(10) 21.3 ns (9.05) | la_vector<12> 37.2 ns (12.23) |
| 5GRAM_3          | C++(9) 12.9 ns (11.31) | Rust(9) 13.3 ns (11.33) | la_vector<11> 18.8 ns (11.69) |
| URL_1            | C++(10) 27.3 ns (7.31) | Rust(10) 26.7 ns (7.34) | la_vector<10> 43.4 ns (10.14) |
| URL_2            | C++(9) 13.6 ns (10.67) | Rust(9) 13.2 ns (10.70) | la_vector<10> 22.4 ns (10.58) |
| URL_3            | C++(9) 10.0 ns (12.44) | Rust(9) 9.4 ns (12.49)  | la_vector<6> 14.1 ns (7.49)   |
| GOV2 avg 10M+    | C++(9) 7.8 ns (9.41)   | Rust(9) 7.2 ns (9.41)   | la_vector<8> 14.6 ns (9.41)   |
| GOV2 avg 1M–10M  | C++(9) 12.5 ns (5.78)  | Rust(9) 11.8 ns (5.79)  | la_vector<8> 17.7 ns (8.11)   |
| GOV2 avg 100K–1M | C++(9) 13.0 ns (4.90)  | Rust(9) 12.5 ns (4.90)  | la_vector<6> 22.4 ns (6.10)   |

Note that it is debatable to compare structures with different constant
parameters. If you choose a fixed best structure for each type, however, the
picture does not change. Moreover, all datasets are rather small—the structures
are all running in some level of the cache, so there is no way from these
benchmarks to judge the behavior of the structures when memory access becomes
expensive.

These results were obtained on a Linux server with an Intel i7-12700KF and 64
GiB of RAM using `gcc` 15.2.1 and Rust 1.95. If you want to reproduce them, you
can follow the instructions in the `README-orig.md` file, which contains the
original instructions for the repository, plus instructions for compiling the
Rust variant. You will need:

- a C++ compiler supporting C++20;
- the Rust compiler;
- CMake, git and curl;
- development libraries for Boost, OpenMP and GTest;
- the Python packages `pandas` and `matplotlib`.

There are a few bonuses: the `plot_results.py` script will generate the graphs
above from the data, and the `gen_table.py` script will generate the table above
for the top-_k_ structures. The `run_all.sh` script has been modified to pin
execution to core 2, which is usually a performance core not shared with the OS
scheduler, but you can change the choice (or not pin at all) by modifying the
`PIN` variable in the script.

Caveats:

- You have to download the author data from the link in `README-orig.md` and
  place them in the `data` directory, where some symlinks will make it usable.
  In the form they are distributed, they cannot be used with the code, as the
  code expects data in a slightly different position for the GOV datasets, but
  the authors have ignored my requests to obtain data corresponding exactly to
  the code in their repository. I did my best to make the data usable.

- The `sux` C++ code has been modified from the 2007 version to make it possible
  to try different parameters; it also uses unaligned reads, which in 2007 were
  not supported reliably, and moves the “cold path” of selection to a separate
  function, as it happens in the Rust version.
