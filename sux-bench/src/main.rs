use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::time::Instant;

use mem_dbg::{MemSize, SizeFlags};
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use sux::dict::EliasFanoBuilder;
use sux::traits::indexed_dict::{IndexedSeq, Pred};
use sux::traits::TryIntoUnaligned;

/// Read the binary format written by write_data_binary<uint32_t, true>:
/// 8-byte little-endian size_t prefix, then n * 4-byte little-endian u32 values.
fn read_data_binary(path: &str) -> io::Result<Vec<u32>> {
    let mut f = File::open(path)?;
    let file_len = f.metadata()?.len();
    let mut size_buf = [0u8; 8];
    f.read_exact(&mut size_buf)?;
    let size = u64::from_le_bytes(size_buf);
    let byte_count = size
        .checked_mul(4)
        .and_then(|b| b.checked_add(8))
        .filter(|&total| total <= file_len)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "{path}: declared size {size} would require more bytes than the file's {file_len}"
                ),
            )
        })?
        - 8;
    let mut bytes = vec![0u8; byte_count as usize];
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

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let files: Vec<String> = std::env::args().skip(1).collect();
    if files.is_empty() {
        eprintln!("Usage: sux-bench <file1> [file2 ...]");
        std::process::exit(1);
    }

    let stdout = io::stdout();
    let mut out = BufWriter::new(stdout.lock());

    writeln!(
        out,
        "filename,n,u,ratio,sux_ef_time_build,sux_ef_time_select,sux_ef_bpk,sux_ef_time_rank"
    )?;

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

        // universe = max value + 1
        let u = *data.last().unwrap() as usize + 1;
        let ratio = n as f64 / u as f64;

        let filename = std::path::Path::new(path)
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or(path.as_str());

        // n and u are printed as plain integers (C++ prints them as size_t / uint32_t,
        // unaffected by std::scientific); ratio is a float in scientific notation.
        write!(out, "'{filename}',{n},{u},{}", fmt_sci(ratio))?;

        // 20% of the dataset is queried, matching TIMES_TEST = data.size() / 5 in C++.
        let times_test = n / 5;

        // Select query indices: uniform in [0, n-2], matching C++ dist1(1, n-1) (1-indexed).
        let mut rng1 = SmallRng::seed_from_u64(42);
        let rands1: Vec<usize> = (0..times_test).map(|_| rng1.gen_range(0..n - 1)).collect();

        // Rank query values: uniform in [data[0], data[n-1]-1], matching C++.
        // gen_range(a..b) is [a, b), so gen_range(data_min..data_max) = [data_min, data_max-1].
        let data_min = data[0] as usize;
        let data_max = data[n - 1] as usize; // exclusive upper bound for gen_range
        let mut rng2 = SmallRng::seed_from_u64(42);
        let rands2: Vec<usize> = (0..times_test)
            .map(|_| rng2.gen_range(data_min..data_max))
            .collect();

        // benchmark a built EF structure (select + rank), return (select_ns, bpk, rank_ns).
        macro_rules! bench_ef {
            ($ef:expr) => {{
                let ef = $ef;
                let bpk = ef.mem_size(SizeFlags::default()) as f64 * 8.0 / n as f64;

                let start = Instant::now();
                let mut cnt: usize = 0;
                for &i in &rands1 {
                    cnt = cnt.wrapping_add(unsafe { ef.get_unchecked(i) });
                }
                core::hint::black_box(cnt);
                let select_ns = start.elapsed().as_nanos() as f64 / times_test as f64;

                let start = Instant::now();
                let mut cnt: usize = 0;
                for &x in &rands2 {
                    cnt = cnt.wrapping_add(ef.rank(x));
                }
                core::hint::black_box(cnt);
                let rank_ns = start.elapsed().as_nanos() as f64 / times_test as f64;
                (select_ns, bpk, rank_ns)
            }};
        }

        let build_start = Instant::now();
        let mut efb = EliasFanoBuilder::new(n, u);
        for &v in &data {
            efb.push(v as usize);
        }
        let ef = efb.build_with_seq_and_dict().try_into_unaligned()?;
        let build_ns = build_start.elapsed().as_nanos() as f64 / n as f64;

        let (sel, bpk, rank) = bench_ef!(ef);

        writeln!(
            out,
            ",{},{},{},{}",
            fmt_sci(build_ns),
            fmt_sci(sel),
            fmt_sci(bpk),
            fmt_sci(rank)
        )?;
    }

    Ok(())
}
