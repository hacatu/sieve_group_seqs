# Group Number Sequence Sieve

Motivated by "[Does the sequence 1,2,3,4,5,6 appear in the number of groups of order n up to isomorphism?](https://math.stackexchange.com/questions/4931487)" (and "[Show that if there are 1,2,3,4 groups of order n+1,n+2,n+3,n+4, then 24∣n](https://math.stackexchange.com/questions/5062059)").

This program finds numbers `n` where `gnu(n+1)=1, ..., gnu(n+7)=7` (sequence length 7), or similar sequences with length 6.

`gnu(n)` is the number of groups of order (size) `n`, up to isomorphism.

The `n` values where a sequence of length 7 exist are a (quite sparse) subset of those where a sequence of length 6 appears,
which are a subset of those where a sequence of length 5 appears, etc.

Sequences of length 8 don't exist, because we can show `n = 144m + 72` for some `m` for 6+ length sequences,
but then `n + 8` is a multiple of 16 and has more than 8 groups.

## Building

Requirements: Rust and Git

```
git clone "https://github.com/hacatu/sieve_group_seq"
cd sieve_group_seq
cargo b --profile release-trapv
```

`release-trapv`: optimized build with runtime overflow checking (use this)
`release-lto`: optimized build without overflow checking (don't use this)

## Running

```
target/release-trapv/sieve_group_seq <sequence length> <n_min> <n_max> <number of times to print a progress update, eg 90/100>
```

Outputs to `result-<sequence length>-<n_min>-<n_max>.json`.  No attempt is made to check if the file already exists or some overlapping file exists.

The program will automatically use all CPU cores + hyperthreading.

To run the program on a supercomputer cluster, just run many instances with sequential `<n_min> <n_max>` bounds (the bounds are inclusive).

For example, to sieve up to `10^16` on 10 computers, run
```
target/release-trapv/sieve_group_seq 7 $((x*10**15 + 1)) $(((x+1)*10**15)) 0
```

on each of the 10 computers in the cluster, where `x` is replaced with numbers from `0` up to `9`.

## Additional Background

Conditions for when gnu(n) = 1 thru 3 are given by Olsson 2006 [Three-group numbers](https://web.math.ku.dk/~olsson/manus/three-group-numbers.pdf).
Conditions for 4 and 5 are given by Miller 1932 [Orders for which there exist exactly Four or Five Groups](https://www.pnas.org/doi/epdf/10.1073/pnas.18.7.511).
Conditions for 6 and 7 are given by Mahmoud 2024 [Orders for which there exist exactly Six or Seven Groups](https://arxiv.org/pdf/2405.04794).

The second paper is pretty ... bad, but the others are good.  The code has explanations of all the cases, which is helpful for understanding 4 and 5.

It is also worth studying holder's formula, which is presented in Mahmoud 2024 and Conway 2008 [Counting groups: gnus, moas and other exotica](https://www.math.auckland.ac.nz/~obrien/research/gnu.pdf).

## Testing

Some tests are included in `src/main.rs`.
To run all tests, do `cargo test`.
To run all tests and generate a coverage report, install `cargo-llvm-cov` and run it:
```
cargo install cargo-llvm-cov
cargo llvm-cov --all-features --html --open
```
DON'T do this unless you want to verify that the tests actually test all the code.
`cargo-llvm-cov` has over 500 dependencies.
But it doesn't actually take too long to install and run the tests.

It will prompt you to install more packages when you run it.

All the `is_gnu_k` functions match all terms in oeis, and for all but `k=7` this ensures full branch coverage.

For `k=7`, I found that case IX from Mahmoud's paper wasn't covered, so I constructed a small example for that case and checked it actually has 7 groups using gap.

See `tests::check_7_case_9` in `main.rs` for more info.

## Future Directions

Currently, only sequences of length 2+ are supported (but sequences of length 8+ don't exist, so the program just immediately returns an empty list).
Three different but very similar sieves are used for 2/3, 4/5, and 6/7.
There are 1, 2, and 3 bad residues respectively mod each prime > 3, and that's the only real difference.

We could also find "sequences" of length 1, which just correspond to `n` where `n+1` has `1` group.
This is just [A003277](https://oeis.org/A003277), so there isn't really a need to support this case.
If we wanted to, we could make a sieve to totally handle this case:

1. Initialize `slot[n] = {g = 1, r = n}` for every odd `n` in the range (the only even `n` with `gnu(n) = 1` is `2`)
2. For every sieving prime `p`, for every multiple `m`:
3. If `m` is divisible by `p^2`, set `slots[m].g = 0`, marking the slot as invalid.
4. Otherwise, multiply `slots[m].g` by `p - 1` and divide `slots[m].r` by `p`.
5. At the end, for every odd `n` in the range, if `slots[n].g != 0` and `slots[n].r != 1`, then `slots[n].r` is the unique large prime factor of `n`, so multiply `slots[n].g` by `slots[n].r - 1`.
6. Now `slots[n].g` holds `phi(n)` for all squarefree `n` and `0` for squareful `n`, so `gnu(n) = 1` iff `slots[n].g != 0` and `gcd(slots[n].g, n) == 1`.

Also, we can't use a standard wheel optimization like `m mod 6` because we've already incorporated the maximal restrictions mod powers of 2 and 3,
but we could do an `m mod 5`, `m mod 35`, or `m mod 385` wheel optimization.  Conveniently, `35` and `385` lead to `8` and `64` coprime residues
respectively.  In theory this could skip 60%-75% of the `m` values, but it would make the code significantly more complicated.

