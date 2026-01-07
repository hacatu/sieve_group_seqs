/*
Group Number Sequence Sieve
Copyright 2026 Gabriel Eiseman
This program is released under the MPL 2 License

This program finds numbers k where gnu(k+1) (the number of non-isomorphic groups of size k+1) is 1,
gnu(k+2) = 2, gnu(k+3) = 3, gnu(k+4) = 4, gnu(k+5) = 5, gnu(k+6) = 6, and gnu(k+7) = 7.

It's impossible for gnu(k+8) to also equal 8, so 7 is the maximum.

conditions for when gnu(n) = 1 thru 3 are given by Olsson 2006 "Three-group numbers"
https://web.math.ku.dk/~olsson/manus/three-group-numbers.pdf

conditions for 4 and 5 are given by Miller 1932 "Orders for which There Exist Exactly Four or Five Groups"
https://www.pnas.org/doi/epdf/10.1073/pnas.18.7.511

conditions for 6 and 7 are given by Mahmoud 2024 "Orders for which there exist exactly six or seven groups"
https://arxiv.org/pdf/2405.04794

it is also worth studying holder's formula, which is presented in the 67 paper and Conway, Dietrich, and O'Brien 2008
"Counting groups: gnus, moas and other exotica"
https://www.math.auckland.ac.nz/~obrien/research/gnu.pdf

The main explanation of the background can be found here:
"Does the sequence 1,2,3,4,5,6 appear in the number of groups of order n up to isomorphism?"
https://math.stackexchange.com/questions/5062059
*/

#![feature(file_buffered)]
mod group_count_filters;
mod segment_sieves;

use std::{borrow::Borrow, cmp::min, fs::File, sync::atomic::{self, AtomicU64}, thread::scope};
use itertools::Itertools as _;
use primesieve_wrapper::{generate_primes, get_num_threads};
use crate::{group_count_filters::*, segment_sieves::*};

/// Get the L2 cache size in bytes, because this is the TOTAL number of bytes all buckets should take up.
/// I tried to have this function return the cache size whenever possible, it should work on "any" linux system,
/// any mac with osx 10+ (but I can't test this), and any x86_64 system.
/// On other systems, we just guess that 512 KB is a good bucket size; feel free to recompile with a different value if running on such a system.
fn get_cache_size() -> Option<u64> {
	let res = {
	#[cfg(any(target_os="linux", target_os="android"))]
	unsafe {
		libc::sysconf(libc::_SC_LEVEL2_CACHE_SIZE)
	}
	#[cfg(any(target_os="macos", target_os="ios"))]
	unsafe {
		let mut res = MaybeUninit::<i64>::zeroed();
		libc::sysctlbyname("hw.perflevel0.l2cachesize\0".as_ptr().cast(), res.as_mut_ptr().cast(), size_of::<i64>(), std::ptr::null_mut(), 0);
		res.assume_init()
	}
	#[cfg(not(any(target_os="linux", target_os="android", target_os="macos", target_os="ios")))]
	{
		#[cfg(target_arch="x86_64")]
		{
			cache_size::l2_cache_size().map_into().unwrap_or_default()
		}
		#[cfg(not(target_arch="x86_64"))]
		{
			0
		}
	}};
	if res <= 0 {
		None
	} else {
		Some(res as _)
	}
}

// this wrapper function forces k (the sequence length) to be a compile time constant, and forces the compiler to produce specialized code
// for every value of K from 2 to 7.
fn sieve_wrapper<const K: u64>(primes: &[u32], m_a: u64, m_b: u64, num_threads: u64, bucket_size: u64, progress_ticks: u64) -> Vec<u64> {
	const BATCH_SIZE: u64 = 200_000;
	let counter = AtomicU64::new(m_a);
	scope(|s|{
		let mut workers = Vec::new();
		for _ in 0..num_threads {
			workers.push(s.spawn(||{
				let mut bucket = unsafe { Box::new_zeroed_slice(bucket_size.div_ceil(8)as _).assume_init() };
				let mut out = Vec::new();
				let mut filtered = Vec::new();
				loop {
					let m_a1 = counter.fetch_add(bucket_size, atomic::Ordering::Relaxed);
					if out.len() as u64 >= BATCH_SIZE || m_a1 >= m_b {
						filtered.extend(out.drain(..).filter(|&n|match K {
							2 => is_gnu_1(n+1),
							3 | 4 => is_gnu_3(n+3)&&is_gnu_1(n+1),
							5 | 6 => is_gnu_3(n+3)&&is_gnu_5(n+5)&&is_gnu_1(n+1),
							7 => is_gnu_3(n+3)&&is_gnu_7(n+7)&&is_gnu_5(n+5)&&is_gnu_1(n+1),
							_ => unreachable!()
						}));
					}
					if m_a1 >= m_b {
						break;
					}
					if progress_ticks != 0 && (m_a1 - m_a + bucket_size)*progress_ticks/(m_b - m_a) != (m_a1 - m_a)*progress_ticks/(m_b - m_a) {
						println!("{}/{progress_ticks}", (m_a1 - m_a)*progress_ticks/(m_b - m_a));
					}
					match K {
						2 | 3 => sieve_segment_2(primes, m_a1, min(m_a1 + bucket_size, m_b), &mut bucket, &mut out),
						4 | 5 => sieve_segment_4(primes, m_a1, min(m_a1 + bucket_size, m_b), &mut bucket, &mut out),
						6 | 7 => sieve_segment_6(primes, m_a1, min(m_a1 + bucket_size, m_b), &mut bucket, &mut out),
						_ => unreachable!()
					}
				}
				filtered
			}));
		}
		workers.into_iter().map(|h|h.join().unwrap()).flatten().collect()
	})
}

/// Sieve to find candidates and then filter them in batches to find all sequences in range.
/// Returns a vec of all n in the range [n_min, n_max] where gnu(n+i) = i for i=1,...,k.
/// Currently, k must be >1.  8+ is impossible so it will just return an empty vector.
/// 1 is just looking for cyclic numbers.
/// This just calls `sieve_segment_6` (or 2 or 4) on all available threads until the whole range has been processed.
/// `progress_ticks` sets how often to print `x/y' to give an indication of progress.
fn find_seqs(n_min: u64, n_max: u64, k: u64, progress_ticks: u64) -> Vec<u64> {
	let num_threads = get_num_threads() as u64;
	let preffered_bucket_size = match get_cache_size() {
		Some(l2_bytes) => {
			println!("Detected L2 cache size {l2_bytes} B");
			l2_bytes * 4
		},
		None => {
			println!("WARNING: Could not detect L2 cache size, falling back to 512 KB/thread");
			1 << 22
		}
	};
	let bucket_size = min(preffered_bucket_size, n_max.isqrt());
	let mut res = match k {
		0 | 1 => panic!("Only k > 1 is supported"),
		2 | 3 => {
			let m_max = n_max/4;
			let m_min = (n_max + 3)/4;
			if m_min > m_max {
				return Vec::new();
			}
			let primes = generate_primes(3, (4*m_max + 1).isqrt()as _);
			if k == 3 {
				sieve_wrapper::<3>(primes.borrow(), m_min, m_max + 1, num_threads, bucket_size, progress_ticks)
			} else {
				sieve_wrapper::<2>(primes.borrow(), m_min, m_max + 1, num_threads, bucket_size, progress_ticks)
			}
		},
		4 | 5 => {
			let m_max = (n_max - 24)/48;
			let m_min = (n_min + 23)/48;
			if m_min > m_max {
				return Vec::new();
			}
			let primes = generate_primes(5, (48*m_max + 25).isqrt()as _);
			if k == 5 {
				sieve_wrapper::<5>(primes.borrow(), m_min, m_max + 1, num_threads, bucket_size, progress_ticks)
			} else {
				sieve_wrapper::<4>(primes.borrow(), m_min, m_max + 1, num_threads, bucket_size, progress_ticks)
			}
		},
		6 | 7 => {
			let m_max = (n_max - 72)/144;
			let m_min = (n_min + 71)/144;
			if m_min > m_max {
				return Vec::new();
			}
			let primes = generate_primes(5, (144*m_max + 73).isqrt()as _);
			if k == 7 {
				sieve_wrapper::<7>(primes.borrow(), m_min, m_max + 1, num_threads, bucket_size, progress_ticks)
			} else {
				sieve_wrapper::<6>(primes.borrow(), m_min, m_max + 1, num_threads, bucket_size, progress_ticks)
			}
		},
		_ => {
			return Vec::new();
		}
	};
	if progress_ticks != 0 {
		println!("{progress_ticks}/{progress_ticks}");
	}
	res.sort();
	res
}

fn main() {
	let [count, n_min, n_max, progress_ticks] = match std::env::args().skip(1).map(|s|s.parse()).try_collect::<_, Vec<_>, _>() {
		Ok(v) => {
			if v.len() != 4 {
				println!("Incorrect number of arguments!");
				println!("Call like `sieve_group_seqs <sequence length> <n_min> <n_max> <progress ticks>', where:");
				println!("  - <sequence length> should be 6 or 7 (longer sequences are impossible and shorter sequences are not currently supported");
				println!("  - <n_min> and <n_max> should be integers with 0 <= n_min <= n_max");
				println!("  - <progress ticks> specifies how many times to print `x/y' while running, to indicate how close to done the program is.  0 disables printing.");
				return;
			}
			v
		},
		Err(e) => {
			println!("Failed to parse argument: {e}.  Run with no arguments for usage information");
			return;
		}
	}[..] else {
		unreachable!()
	};
	if count < 2 {
		println!("Counts above 7 are impossible and counts below 2 are not currently supported");
		return;
	}
	if n_max < n_min {
		println!("Invalid range, n_min should be <= n_max");
		return;
	}
	if progress_ticks > 1000 && progress_ticks > (n_max - n_min)/(144 << 22) {
		println!("progress_ticks is unreasonably high and probably higher than the actual number of segments; limit it to max(1000, (n_min - n_max)/{})", 144 << 22);
		return;
	}
	let res = find_seqs(n_min, n_max, count, progress_ticks);
	println!("Found {} length {count} sequences from {n_min} to {n_max}", res.len());
	let filename = format!("result-{count}-{n_min}-{n_max}.json");
	serde_json::to_writer(File::create_buffered(&filename).unwrap(), &res).unwrap();
	println!("Saved to '{filename}'");
}

#[cfg(test)]
mod tests {
    use std::{fs::{self, File}, io::BufRead};
    use itertools::Itertools as _;

    use crate::{find_seqs, is_gnu_1, is_gnu_2, is_gnu_3, is_gnu_4, is_gnu_5, is_gnu_6, is_gnu_7};
	static FILTERS: [fn(u64) -> bool; 7] = [is_gnu_1 as _, is_gnu_2 as _, is_gnu_3 as _, is_gnu_4 as _, is_gnu_5 as _, is_gnu_6 as _, is_gnu_7 as _];

	fn get_oeis_seq_cached(a: u64) -> impl Iterator<Item=u64> {
		let filename = format!("resources/b{a:06}.txt");
		if !fs::exists(&filename).unwrap_or_default() {
			fs::write(&filename, ureq::get(format!("https://oeis.org/A{a:06}/b{a:06}.txt")).call().unwrap().body_mut().read_to_string().unwrap()).unwrap();
		}
		File::open_buffered(&filename).expect("file should exist").lines().filter_map(|line|{
			let line = line.expect("file should be readable");
			line.split(' ').skip(1).filter_map(|s|u64::from_str_radix(s, 10).ok()).next()
		})
	}

	#[test]
	fn check_gnu_filters() {
		let a_nums = [
			3277, // 1
			54395, // 2
			55561, // 3
			54396, // 4
			54397, // 5
			135850, // 6
			249550, // 7 (list is short)
		];
		for (i, a) in a_nums.into_iter().enumerate() {
			println!("Checking gnu{}", i+1);
			let mut last_good = 0;
			for (expected, observed) in get_oeis_seq_cached(a).zip((1..).filter(|&n|FILTERS[i](n))) {
				assert!(expected == observed, "gnu{} mismatch after {last_good}: expected {expected} but got {observed}", i+1);
				last_good = expected;
			}
		}
	}

	#[test]
	fn check_existing_oeis() {
		let a_nums = [373648, 373649, 373650]; // length 2-4 sequences (5 has a mistake and I already checked it by hand)
		for (i, a) in a_nums.into_iter().enumerate() {
			let i = i + 2;
			let expected_seq = get_oeis_seq_cached(a).collect_vec();
			let c = expected_seq.len();
			let observed_seq = find_seqs(0, expected_seq.iter().max().copied().unwrap(), i as _, 0);
			let mut last_good = 0;
			for (e, o) in expected_seq.into_iter().zip(observed_seq) {
				assert!(e == o, "after {last_good}: expected next {i}-sequence should be {e} but got {o}");
				last_good = e;
			}
			println!("all {c} terms matched for {i}-sequence");
		}
	}

	#[test]
	fn check_7_case_9() {
		// Cases I-VIII are checked by the oeis data, but the smallest example for case IX is larger than the max listed for that sequence,
		// so I constructed an example (not necessarily the smallest):
		// recall IX: n = p^2qrs where (p + 1) is a multiple of q, (p + 1) is a multiple of r, and (s - 1) is a multiple of p but not p^2
		// so pick q = 3, r = 5, p = 29, s = 29, check it meets these constraints, and compute n = 744285.
		// Then I double checked this in gap using the cubefree package https://github.com/gap-packages/cubefree
		// On arch linux, this is part of the `gap-packages` package, so install it, run gap, and run `LoadPackage("cubefree");NumberCFGroups(744285);`
		// Obviously this would not help for non-cubefree n
		assert!(is_gnu_7(744285));
	}

	#[test]
	fn spot_check() {
		for n in [5973822114120,6305771634120,28058687347320,48414128744520,74556478687320,84300170172120,142366076070120,153432090884520,207916418382120,339858901929960,375278016786120,415728077892120,426239950426920] {
			assert!(FILTERS.iter().enumerate().all(|(i, f)|f(n+i as u64+1)));
		}
	}
}

