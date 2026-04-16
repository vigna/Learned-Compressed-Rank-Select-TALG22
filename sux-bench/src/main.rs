use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::time::Instant;

use mem_dbg::{MemSize, SizeFlags};
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use sux::dict::EliasFanoBuilder;
use sux::rank_sel::{SelectAdaptConst, SelectZeroAdaptConst};
use sux::traits::indexed_dict::{IndexedSeq, Pred};
use sux::traits::TryIntoUnaligned;

/// Read the binary format written by write_data_binary<uint32_t, true>:
/// 8-byte little-endian size_t prefix, then n * 4-byte little-endian u32 values.
fn read_data_binary(path: &str) -> io::Result<Vec<u32>> {
    let mut f = File::open(path)?;
    let mut size_buf = [0u8; 8];
    f.read_exact(&mut size_buf)?;
    let size = usize::from_le_bytes(size_buf);
    let mut bytes = vec![0u8; size * 4];
    f.read_exact(&mut bytes)?;
    Ok(bytes
        .chunks_exact(4)
        .map(|b| u32::from_le_bytes(b.try_into().unwrap()))
        .collect())
}

/// Format a float in C++ std::scientific + setprecision(3) style: "1.234e+02".
/// Rust's {:.3e} omits the + sign and doesn't zero-pad the exponent.
fn fmt_sci(x: f64) -> String {
    let s = format!("{:.3e}", x);
    let e_pos = s.find('e').unwrap();
    let mantissa = &s[..e_pos];
    let exp_str = &s[e_pos + 1..];
    let (sign, digits) = if let Some(rest) = exp_str.strip_prefix('-') {
        ("-", rest)
    } else if let Some(rest) = exp_str.strip_prefix('+') {
        ("+", rest)
    } else {
        ("+", exp_str)
    };
    let exp_num: u32 = digits.parse().unwrap_or(0);
    format!("{mantissa}e{sign}{exp_num:02}")
}

fn main() {
    let files: Vec<String> = std::env::args().skip(1).collect();
    if files.is_empty() {
        eprintln!("Usage: sux-bench <file1> [file2 ...]");
        std::process::exit(1);
    }

    let stdout = io::stdout();
    let mut out = BufWriter::new(stdout.lock());

    writeln!(out, "filename,n,u,ratio,sux_ef_time_select,sux_ef_bpk,sux_ef_time_rank").unwrap();

    for path in &files {
        let data = match read_data_binary(path) {
            Ok(d) => d,
            Err(e) => {
                eprintln!("Error reading {path}: {e}");
                continue;
            }
        };

        let n = data.len();
        if n < 5 {
            eprintln!("Skipping {path}: too few elements ({n})");
            continue;
        }

        // universe = max value + 1 (same as C++)
        let u = *data.last().unwrap() as u64 + 1;
        let ratio = n as f64 / u as f64;

        let filename = std::path::Path::new(path)
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or(path.as_str());

        // n and u are printed as plain integers (C++ prints them as size_t / uint32_t,
        // unaffected by std::scientific); ratio is a float in scientific notation.
        write!(out, "'{filename}',{n},{u},{}", fmt_sci(ratio)).unwrap();

        // 20% of the dataset is queried, matching TIMES_TEST = data.size() / 5 in C++.
        let times_test = n / 5;

        // Select query indices: uniform in [1, n-1] (1-indexed, matching C++).
        // Seed 2323 matches mt1 in my_benchmark.cpp.
        let mut rng1 = SmallRng::seed_from_u64(2323);
        let rands1: Vec<usize> = (0..times_test).map(|_| rng1.gen_range(1..n)).collect();

        // Rank query values: uniform in [data[0], data[n-1]-1], matching C++.
        // Seed 4242 matches mt2 in my_benchmark.cpp.
        // gen_range(a..b) is [a, b), so gen_range(data_min..data_max) = [data_min, data_max-1].
        let data_min = data[0];
        let data_max = data[n - 1]; // exclusive upper bound for gen_range
        let mut rng2 = SmallRng::seed_from_u64(4242);
        let rands2: Vec<u32> = (0..times_test)
            .map(|_| rng2.gen_range(data_min..data_max))
            .collect();

        // Helper: benchmark a built EF structure (select + rank), return (select_ns, bpk, rank_ns).
        macro_rules! bench_ef {
            ($ef:expr) => {{
                let ef = $ef;
                let bpk = ef.mem_size(SizeFlags::default()) as f64 * 8.0 / n as f64;
                let start = Instant::now();
                let mut acc = 0u64;
                for &i in &rands1 {
                    acc = acc.wrapping_add(ef.get(i - 1));
                }
                std::hint::black_box(acc);
                let select_ns = start.elapsed().as_nanos() as f64 / times_test as f64;
                let start = Instant::now();
                let mut acc = 0u64;
                for &x in &rands2 {
                    acc = acc.wrapping_add(ef.rank(x as u64) as u64);
                }
                std::hint::black_box(acc);
                let rank_ns = start.elapsed().as_nanos() as f64 / times_test as f64;
                (select_ns, bpk, rank_ns)
            }};
        }

        // Build base EF (no select support) then wrap with each parameterization.
        let build_base = || {
            let mut efb = EliasFanoBuilder::new(n, u);
            for &v in &data {
                efb.push(v as u64);
            }
            efb.build()
        };

        let ef = unsafe {
            build_base()
                .map_high_bits(SelectAdaptConst::<_, _>::new)
                .map_high_bits(SelectZeroAdaptConst::<_, _>::new)
        }
        .try_into_unaligned()
        .unwrap();
        let (sel, bpk, rank) = bench_ef!(ef);

        writeln!(out, ",{},{},{}", fmt_sci(sel), fmt_sci(bpk), fmt_sci(rank)).unwrap();
    }
}
