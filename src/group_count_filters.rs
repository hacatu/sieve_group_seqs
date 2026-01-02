/*
Factorization based filters to check if n has k non-isomorphic groups, for k=1,...,7
Copyright 2026 Gabriel Eiseman
This program is released under the MPL 2 License

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



For squarefree n, gnu(n) can be computed using Holder's formula, which encodes a connection between the number of groups of size
n and the Holder graph for n.

The Holder graph is a directed graph with labelled nodes and unlabelled edges.  Each node is labelled with a prime, and there is an edge from
p to q iff p | q - 1.

The number of groups of size n is the sum over (all subsets of the Holder graph's nodes where every node in the subset has at least one out-neighbor
outside the subset) of (product over all nodes p in subset of (p^(number of p's out-neighbors outside of the subset) - 1)/(p-1)).

So you can see that we don't actually have to restrict the subsets we sum over, it's just that subsets with at least one node that DOESN'T have an
out-neighbor outside the subset will contribute 0.

This leads to all sorts of lower bounds and is where a ton of the analysis in the papers comes from.

We can reconstruct all the conditions for squarefree n for gnu(n) = k that the papers give.



For cubefree n, all is not yet lost.  Mahmoud's paper gives a good explanation of the extended Holder graph, where nodes are prime powers,
so if p^2 | n we have a node labelled p^2 instead of p.  This graph also has labelled edges: an edge from p to q is normal ("strong"),
but an edge from p^2 to q can be strong if p^2 | q-1 or weak if p | q-1 but p^2 does not.  And an edge from q to p^2 is strong if q | p-1 or
weak if q | p+1.  I'm not sure about edges from p^2 to q^2.



For n that are squarefree aside from a single prime cube, we need to understand at least roughly what relations or p | q + 1 when q^2 | n means,
and then understand the groups of order p^5.

Ordinary relations correspond to when we can form nontrivial semidirect products when a group of size p acts on a group of size q,
because q-1 is the size of the automorphism group of the cyclic group Cq.  Similarly, if a prime p divides q+1, then a nontrivial sdp exists
with Cp acting on C(q^2) but not when acting on Cq x Cq.  If p divides q^2 + q + 1, then a nontrivial sdp exists with Cp acting on C(q^3).
The p | q - 1 style relations are encoded by the Holder graph, which is why it can be used to compute gnu(n) for all squarefree n.
Similarly, the extended Holder graph whose nodes are p or p^2 and edges can be "weak" (ie p divides q-1 but p^2 does not, or q divides p+1 but not p-1)
can be used to compute gnu(n) for all cubefree n (see Mahmoud 2024).  I'm not sure if this generalizes more, but we only have to deal with very limited
cube cases in this analysis.
There are 5 groups of size p^3 (see https://kconrad.math.uconn.edu/blurbs/grouptheory/groupsp3.pdf):
- C(p^3), the cyclic group of order p^3
- Cp x C(p^2)
- Cp x Cp x Cp
- Heisenberg(Z_p), the Heisenberg group mod p, which we can represent using three 3x3 matrices x = [[1, 1, 0], [0, 1, 0], [0, 0, 1]], y = [[1, 0, 1], [0, 1, 0], [0, 0, 1]], and z = [[1, 0, 0], [0, 1, 1], [0, 0, 1]].
- Cp sdp C(p^2), (sdp means semi-direct product, since here we have a nontrivial sdp)
The first three groups are abelian, while the last two are not.  Also, the last two are slightly different for p=2.
We can count the number of automorphisms for these using a similar argument to https://en.wikipedia.org/wiki/Group_isomorphism#Automorphisms,
and we get that they have p^2(p-1), p^3(p-1)^2, (p-1)(p^2 + p + 1), p^3(p-1)^3(p+1), and p^3(p-1)^2 respectively.
*/

#![allow(dead_code)]
use std::mem::MaybeUninit;
use flint_sys::ulong_extras::{n_factor, n_factor_init, n_factor_t};
use num_integer::gcd;


/// Compute "reduced phi", ie if p^a is the power of p dividing n, it contributes p-1 instead of p^(a-1)(p-1).
fn red_phi(fxn: &n_factor_t) -> u64 {
	fxn.p[..fxn.num as usize].iter().map(|&p|p-1).product()
}

/// Get a bitmask of the indexes j where fxn.p[j] divides fxn.p[i]-1.
/// This gets the out edges of a node in the holder graph in a bit-packed format.
/// Note that fxn.p is guaranteed to be ascending.
/// The bitmask is 16 bits because a 64 bit number can have at most 15 distinct prime factors.
fn mult_mask(fxn: &n_factor_t, i: usize) -> u16 {
	if i >= fxn.num as usize {
		return 0;
	}
	let p = fxn.p[i];
	(i+1..fxn.num as usize).filter(|&j|(fxn.p[j]-1).is_multiple_of(p)).fold(0, |a,j|a|1<<j)
}

/// Get a vec of (i, mask) pairs for all primes p = fxn.p[i] that have relations, ie where `mult_mask` does not return 0.
/// This gets the whole holder graph in a bit-packed format.
/// * `g` - should be `red_phi` and is used to identify the primes that have out edges, since p has out edges iff it divides g.
/// * `p_max` - the max number of primes that should be allowed to have out edges.  If more than this many are found, return an empty vec.
/// * `r_max` - the max number of relations that should be found.  If more than this are found, return an empty vec.
/// This function is used in many of the `is_gnu_k` functions, when we know that a holder graph with too many non-sink nodes or edges would definitely not be valid.
fn gather_mults(fxn: &n_factor_t, mut g: u64, mut p_max: usize, mut r_max: usize) -> Vec<(u8, u16)> {
	let mut res = Vec::with_capacity(p_max);
	if p_max == 0 {
		return res;
	}
	for (i, &p) in fxn.p[..fxn.num as usize].iter().enumerate() {
		if g.is_multiple_of(p) {
			g /= p;
			if p_max == 0 {
				res.clear();
				return res;
			}
			p_max -= 1;
			let mask = mult_mask(fxn, i);
			r_max = match r_max.checked_sub(mask.count_ones() as _) {
				Some(v) => v,
				None => { res.clear(); return res; }
			};
			res.push((i as _, mask));
		}
	}
	res
}

/// Check if gnu(n) = 1.
/// Per Olsson 2006, this happens iff n is squarefree and has no relations
pub fn is_gnu_1(n: u64) -> bool {
	unsafe {
		let mut fxn = MaybeUninit::uninit();
		n_factor_init(fxn.as_mut_ptr());
		let mut fxn = fxn.assume_init();
		n_factor(&mut fxn, n, 0);
		let mut phi = 1;
		for i in 0..fxn.num as usize {
			if fxn.exp[i] != 1 {
				return false;
			}
			phi *= fxn.p[i] - 1;
		}
		gcd(phi, n) == 1
	}
}

/// Check if gnu(n) = 2.
/// Per Olsson 2006, this happens iff on of the following cases holds:
/// * I: n is squarefree and has exactly one relation p | q-1
/// * II: n is squarefree aside from one square prime factor p, has no relations, and has no prime factor dividing p+1
pub fn is_gnu_2(n: u64) -> bool {
	unsafe {
		let mut fxn = MaybeUninit::uninit();
		n_factor_init(fxn.as_mut_ptr());
		let mut fxn = fxn.assume_init();
		n_factor(&mut fxn, n, 0);
		let mut q = None;
		let mut phi_red = 1;
		for i in 0..fxn.num as usize {
			match (fxn.exp[i], q) {
				(1, _) => (),
				(2, None) => q = Some(fxn.p[i]),
				_ => return false
			}
			phi_red *= fxn.p[i]-1;
		}
		let g = gcd(n, phi_red);
		match q {
			None => {
				for (i, &p) in fxn.p[..fxn.num as usize].iter().enumerate() {
					if p == g {
						return mult_mask(&fxn, i).count_ones() == 1;
					} else if (g).is_multiple_of(p) {
						return false;
					}
				}
				false
			},
			Some(q) => g == 1 && fxn.p[..fxn.num as usize].iter().all(|&p|!(q+1).is_multiple_of(p))
		}
	}
}

/// Check if gnu(n) = 3.
/// Per Olsson 2006, this happens iff one of the following cases holds:
/// * I: n is squarefree and has exactly two relations p | q - 1 and q | r - 1
/// * II: n is squarefree aside from one square prime factor p, has no relations, and has exactly one prime divisor dividing p + 1.
pub fn is_gnu_3(n: u64) -> bool {
	unsafe {
		let mut fxn = MaybeUninit::uninit();
		n_factor_init(fxn.as_mut_ptr());
		let mut fxn = fxn.assume_init();
		n_factor(&mut fxn, n, 0);
		let mut q = None;
		let mut phi_red = 1;
		for i in 0..fxn.num as usize {
			match (fxn.exp[i], q) {
				(1, _) => (),
				(2, None) => q = Some(fxn.p[i]),
				_ => return false
			}
			phi_red *= fxn.p[i]-1;
		}
		let mut g = gcd(n, phi_red);
		match q {
			None => {
				// we need to have a chain of divisibility p1 | p2 - 1, p2 | p3 - 1, with p1 and p2 not dividing anything else
				// and p1 * p2 = gcd(n, phi_red)
				let mut remaining = 2;
				for (i, &p) in fxn.p[..fxn.num as usize].iter().enumerate() {
					if p == g {
						return remaining == 1 && mult_mask(&fxn, i).count_ones() == 1;
					} else if g.is_multiple_of(p) {
						g /= p;
						if remaining != 2 || !(g-1).is_multiple_of(p) || mult_mask(&fxn, i).count_ones() != 1 {
							return false;
						}
						remaining -= 1;
					}
				}
				false
			},
			Some(q) => g == 1 && fxn.p[..fxn.num as usize].iter().copied().filter(|&p|(q+1).is_multiple_of(p)).count() == 1
		}
	}
}

/// Check if gnu(n) = 4.
/// Per Miller 1932, this happens iff one of the following cases holds (differently named primes are all different):
/// * I: n is squarefree and there are exactly two relations p | q - 1 and r | q - 1
/// * II: n is squarefree and there are exactly two relations p | q - 1 and r | s - 1
/// * III: n = 2pq with p and q unrelated odd primes (2 is always related to any odd prime)
/// * IV: n is squarefree aside from two square factors p and q, there are no relations, and p+1 and q+1 are not divisible by any prime divisor of n
/// * V: n is squarefree aside from one square factor p, there is exactly one relation p | q - 1, p^2 does not divide q-1, and p+1 is not divisible by any prime divisor of n
/// * VI: n is squarefree aside from one square factor p, there is exactly one relation q | r - 1, and no prime divides p+1
pub fn is_gnu_4(n: u64) -> bool {
	unsafe {
		let mut fxn = MaybeUninit::uninit();
		n_factor_init(fxn.as_mut_ptr());
		let mut fxn = fxn.assume_init();
		n_factor(&mut fxn, n, 0);
		let mut sqmask = 0u32;
		for (i, &e) in fxn.exp[..fxn.num as usize].iter().enumerate() {
			if e > 2 {
				return false;
			}
			if e == 2 {
				if sqmask.count_ones() == 2 {
					return false;
				}
				sqmask |= 1 << i;
			}
		}
		let g = gcd(red_phi(&fxn), n);
		if sqmask == 0 {
/*
- case I: n is squarefree and there are exactly two relations p | q - 1 and r | q - 1
- case II: n is squarefree and there are exactly two relations p | q - 1 and r | s - 1
- case III: n = 2pq with p and q unrelated odd primes (2 is always related to any odd prime)
*/
			if fxn.num == 3 && n&1 == 0 && !(fxn.p[2]-1).is_multiple_of(fxn.p[1]) { // case III
				return true;
			}
			let vvv = gather_mults(&fxn, g, 2, 2);
			if vvv.len() != 2 {
				return false;
			}
			let m = vvv[0].1 | vvv[1].1;
			if m.is_power_of_two() { // case I, both the prime at index vvv[0].0 and vvv[1].0 divide only the prime corresponding to mask m
				return true;
			}
			// otherwise, we have 2 relations to 2 different rhs's, so we need to check neither lhs occurs as an rhs (p, q, r, s are all distinct)
			(m | (1 << vvv[0].0) | (1 << vvv[1].0)).count_ones() == 4
		} else if sqmask.count_ones() == 2 {
/*
- case IV: n is squarefree aside from two square factors p and q, there are no relations, and p+1 and q+1 are not divisible by any prime divisor of n
           (so compute gcd((p1+1)(p2+1), n) and check that it is 1)
*/
			g == 1 &&
			gcd((fxn.p[sqmask.trailing_zeros()as usize] + 1)*(fxn.p[31-sqmask.leading_zeros()as usize] + 1), n) == 1
		} else {
/*
- case V: n is squarefree aside from one square factor p, there is exactly one relation p | q - 1, p^2 does not divide q-1, and p+1 is not divisible by any prime divisor of n
- case VI: n is squarefree aside from one square factor p, there is exactly one relation q | r - 1, and no prime divides p+1
*/
			let i1 = sqmask.trailing_zeros()as usize;
			let p1 = fxn.p[i1];
			if g != 1 && (p1*p1).is_multiple_of(g) {
				// (only case V is possible)
				let mut found_pi = false;
				for i in 0..fxn.num as usize {
					if i == i1 {
						continue;
					}
					if (p1 + 1).is_multiple_of(fxn.p[i]) {
						return false;
					}
					if (fxn.p[i]-1).is_multiple_of(p1) {
						if found_pi || (fxn.p[i]-1).is_multiple_of(p1*p1) {
							return false;
						}
						found_pi = true;
					}
				}
				found_pi
			} else {
				// (only case VI is possible)
				if fxn.p[..fxn.num as usize].iter().any(|&p|(p1+1).is_multiple_of(p)) {
					return false;
				}
				let vvv = gather_mults(&fxn, g, 1, 1);
				vvv.len() == 1 && vvv[0].1.trailing_zeros() != i1 as u32
			}
		}
	}
}

/// Check if gnu(n) = 5
/// Per Miller 1932, this happens iff one of the following holds (differently named primes are all different):
/// * I: n is squarefree with p = 3, p | q - 1, p | r - 1, and no other relations
/// * II: n is squarefree with p | q - 1, r | q - 1, and q | s - 1, and no other relations
/// * III: n is squarefree with p | q - 1, q | r - 1, and r | s - 1, and no other relations
/// * IV: n is squarefree aside from one cube factor p, there are no relations, and no prime divisor of n divides p+1 or p^2 + p + 1
/// * V: n = 2*p^2
/// * VI: n = 12
/// * VII: n is squarefree aside from one square factor p, p divides q - 1, no other relations exist, and EITHER
///        p^2 | q - 1 and no prime divisor of n divides p + 1, OR p^2 does not divide q - 1 and exactly one prime divisor of n divides p + 1
/// * VIII: n is squarefree aside from one square factor p, no relations exist, and p + 1 is divisible by exactly two prime divisors of n
/// * IX: n is squarefree aside from one square factor p, q | r - 1, no other relations, and r | p + 1 but no other prime divisor of n divides p + 1
pub fn is_gnu_5(n: u64) -> bool {
	unsafe {
		let mut fxn = MaybeUninit::uninit();
		n_factor_init(fxn.as_mut_ptr());
		let mut fxn = fxn.assume_init();
		n_factor(&mut fxn, n, 0);
		let mut i1 = None;
		for (i, &e) in fxn.exp[..fxn.num as usize].iter().enumerate() {
			if e > 3 {
				return false;
			}
			if e > 1 {
				if i1.is_some() {
					return false;
				}
				i1 = Some(i);
			}
		}
		let g = gcd(red_phi(&fxn), n);
		let Some(i1) = i1 else {
/*
here, n is squarefree, so there are three cases:
- case I: n is squarefree with p = 3, p | q - 1, p | r - 1, and no other relations
- case II: n is squarefree with p | q - 1, r | q - 1, and q | s - 1, and no other relations
- case III: n is squarefree with p | q - 1, q | r - 1, and r | s - 1, and no other relations
*/
			if g == 3 { // only case I is possible, since case II and III need g to have 3 prime divisors
				return fxn.p[..fxn.num as usize].iter().copied().filter(|&p|(p-1).is_multiple_of(3)).count() == 2;
			}
			let vvv = gather_mults(&fxn, g, 3, 3);
			let [(_, m1), (i2, m2), (i3, _)] = vvv[..] else {
				return false;
			};
			// EITHER m1 = 1 << i3 AND m2 = 1 << i3 (case II)
			// OR m1 = 1 << i2 AND m2 = 1 << i3 (case III)
			return m2 == 1 << i3 && (m1 & ((1 << i2) | (1 << i3)) != 0);
		};
		let p1 = fxn.p[i1];
		if fxn.exp[i1] == 3 {
/*
here, n is squarefree aside from one cube factor
- case IV: n is squarefree aside from one cube factor p, there are no relations, and no prime divisor of n divides p+1 or p^2 + p + 1
*/
			g == 1 && fxn.p[..fxn.num as usize].iter().all(|&p|!((p1+1)*(p1*p1+p1+1)).is_multiple_of(p))
		} else {
/*
otherwise, in the last 5 cases, n is squarefree aside from one square factor p
- case V: n = 2*p^2
- case VI: n = 12
- case VII: p | q - 1, there are no other relations, and EITHER
            (p^2 | q - 1 and no prime divisor of n divides p + 1),
			OR (p^2 does not divide q - 1, and exactly one prime divisor of n divides p + 1)
- case VIII: there are no relations and p + 1 is divisible by exactly 2 prime divisors of n
- type IX: p + 1 divisible by exactly 1 prime divisor q of n, and r | q - 1, with no other relations
*/
			if n == 2*p1*p1 || n == 12 {
				return true; // case V or VI
			}
			if g != 1 && (p1*p1).is_multiple_of(g) { // all relations are p1 | pi - 1, so check case VII
				// need to check that p1 | pi - 1 for a UNIQUE pi,
				// and then EITHER
				// p1^2 | pi - 1, AND no OTHER p | p1 + 1 (but pi would be allowed to), OR
				// p1^2 !| pi - 1, AND "ANY" p | p1 + 1
				let mut pi_mask = 0u32;
				let mut other_div = false;
				for (i, &p) in fxn.p[..fxn.num as usize].iter().enumerate() {
					if (p-1).is_multiple_of(p1) {
						pi_mask |= 1 << i;
						if !pi_mask.is_power_of_two() {
							break;
						}
					} else if (p1+1).is_multiple_of(p) {
						other_div = true;
					}
				}
				if pi_mask.is_power_of_two() {
					let pi = fxn.p[pi_mask.trailing_zeros() as usize];
					if (pi - 1).is_multiple_of(p1*p1) {
						if !other_div {
							return true;
						}
					} else {
						if other_div || (p1+1).is_multiple_of(pi) {
							return true;
						}
					}
				}
			}
			if g == 1 { // there are no relations, so check type 8
				if fxn.p[..fxn.num as usize].iter().copied().filter(|&p|(p1+1).is_multiple_of(p)).count() == 2 {
					return true; // type 8
				}
			}
			/*
			otherwise, check case IX, so we need to test
			- is p1 + 1 divisible by exactly one pi
			- is g equal to one of the primes pj
			- does pj divide pi - 1 but no other prime
			*/
			let mut pi_mask = 0u32; // indices of primes pi that divide p1 + 1
			let mut pj_mask = 0u32; // indices of primes pj where pj-1 is divisible by g. if g = pj and this mask = pi_mask, then the conditions are met
			let mut g_found = false;
			for (i, &p) in fxn.p[..fxn.num as usize].iter().enumerate() {
				if (p1+1).is_multiple_of(p) {
					pi_mask |= 1 << i;
				}
				if (p-1).is_multiple_of(g) {
					pj_mask |= 1 << i;
				}
				if g == p {
					g_found = true;
				}
			}
			g_found && pi_mask.is_power_of_two() && pj_mask == pi_mask
		}
	}
}

/// Find the path in the Holder graph from the ith prime to a sink.
/// For example, if pi | pj - 1, pj | pk - 1, and px | py - 1, with no other relations,
/// then the path starting at the ith prime goes i, j, k.
/// This function returns (last index, bitmask of all indices besides the last).
/// So if pi has no relations, it would return (i, 0).
/// SHOULD NOT be called if pi or some other prime in the chain is related to more than one other prime.
/// Generally we know the graph has this property if we call `gather_mults(p_max, r_max)` and find the actual returned length equals `r_max`,
/// since then each node has at most one relation by the pidgeonhole principle.
/// `masks` is a slice of bitmasks of out-edges for every prime, ie if pi | pj - 1, then the jth bit of `masks[i]` should be 1.
/// This slice is generally build by walking the result of `gather_mults` and populating a length 15 array (since 64 bit n has at most 15 distinct prime divisors).
/// This is only used in a couple places so there is no helper function to do that.
fn trace_chain(mut i: u8, masks: &[u16]) -> (u8, u16) {
	let mut visited_mask = 0;
	loop {
		let m = masks[i as usize];
		if m == 0 {
			return (i, visited_mask);
		}
		visited_mask |= 1 << i;
		debug_assert!(m.is_power_of_two());
		i = m.trailing_zeros() as _;
	}
}

/// Find all prime divisors of n that divide p+1 and return a bitmask of their indices.
/// This is similar to `mult_mask`, except that returns a mask of normal out-edges, whereas this returns a mask of weak in-edges.
/// A normal edge is from p to q when p | q - 1.  A weak edge in this case is from p to q when p | q + 1 and q^2 divides n.
/// The concept of weak edges is from Mahmoud 2024's extended Holder graphs, and p || q - 1, that is, p | q - 1 but p^2 does not (when p^2 divides n)
/// is also called a weak edge.
/// But this function is mainly a helper to check that the prime divisors of n which divide p + 1 look the way we want in whatever case we are testing.
fn weak_inedge_mask(fxn: &n_factor_t, p: u64) -> u16 {
	let mut res = 0;
	for (i, &r) in fxn.p[..fxn.num as usize].into_iter().enumerate() {
		if (p + 1).is_multiple_of(r) {
			res |= 1 << i;
		}
	}
	res
}

/// Check if gnu(n) = 6
/// Per Mahmoud 2024, there are "7" cases, but really there are 12 if we count all subcases.
/// 
/// This paper explains the conditions for "cyclic-free" n, which means in the following table, when we say "n=2pq" for example,
/// we mean "n=2pqm where m is squarefree and has no relations (ie, if a is a prime dividing the mentioned part, and b is a prime dividing m,
/// or both a and b divide m, then neither a∣b−1 nor b∣a−1)
/// (this includes half relations, that is, we should not have b∣a+1 if a2 divides the mentioned part,
/// and we should not have b∣a2+a+1 if a3 divides the mentioned part)".
/// 
/// Also within this table, s, p, q, and r are distinct odd primes, and there are NO relations (as described in my comment above) other than those mentioned.
/// gnu(n) = 6 iff one of the following:
/// * I: n=pqrs where q|r−1, s|r−1, and p|q−1
/// * II: n=2pq where p|q−1
/// * III: n=3p^2 where p is not 3 and 3|p−1
/// * IV: n=p^3q where q|p^2+p+1
/// * V.a: n=p^2qr where p||q−1 (that is, p divides but p^2 does not divide) and q|r−1
/// * V.b: n=p^2qr where q=3, 3|p+1, and 3|r−1 (the paper omits this last condition due to a typo)
/// * V.c: n=p^2qr where p^2|r−1 and q|p+1
/// * VI: n=p^2q^2 where p||q+1
/// * VII: n=n1n2 where n1 and n2 are arithmetically independent (have no relations as described above) and each has two independent choices:
///     * a: n1=q1r1 where q1|r1−1; OR n1=q1^2
///     * b: n2=p2q2r2 where q2|r2−1 and p1|q1−1; OR n2=p2^2q2 and q2∣p2+1
/// (I refer to the choices left, left as "VII.1"; left, right as "VII.2"; right, left as "VII.3"; and right, right as "VII.4")
pub fn is_gnu_6(n: u64) -> bool {
	unsafe {
		let mut fxn = MaybeUninit::uninit();
		n_factor_init(fxn.as_mut_ptr());
		let mut fxn = fxn.assume_init();
		n_factor(&mut fxn, n, 0);
		let mut sqmask = 0u16;
		let mut sqmask_pow = 0;
		for (i, &e) in fxn.exp[..fxn.num as usize].iter().enumerate() {
			if e > 3 {
				return false;
			}
			if e > 1 {
				if sqmask_pow == 0 {
					sqmask_pow = e;
				} else if sqmask_pow != e {
					return false;
				}
				if sqmask.count_ones() == 2 {
					return false;
				}
				sqmask |= 1 << i;
			}
		}
		let g = gcd(red_phi(&fxn), n);
		if sqmask_pow == 0 { // squarefree cases
/*
I: n = pqrs where (r - 1) is divisible by q and s, and (q - 1) is divisible by p
II: n = 2pq where (q - 1) is divisible by p
VII.1: n = n1 n2 with n1, n2 arithmetically independent
  and n1 = q1r1 and n2 = p2q2r2 where (r1 - 1) is divisible by q1 and (r2 - 1) is divisible by q2 and (q2 - 1) is divisible by p2

other prime divisors can exist, BUT no other relations (p | q-1)
*/
			let vvv = gather_mults(&fxn, g, 3, 3);
			if fxn.num == 3 && fxn.p[0] == 2 && (fxn.p[2] - 1).is_multiple_of(fxn.p[1]) {
				return true; // case II
			}
			if vvv.len() != 3 {
				return false;
			}
			let mut masks = [0u16; 15];
			for &(ii, mi) in &vvv {
				masks[ii as usize] = mi;
			}
			let (i1, visited_mask) = trace_chain(vvv[0].0, &masks[..]);
			// trace_chain follows out edges until it reaches a sink, with i1 being the sink and visited_mask being the pre-sinks
			if visited_mask.count_ones() == 3 {
				return false; // everything is in one chain
			} // otherwise, at least one node was not in the chain
			for &(i, _) in &vvv[1..] {
				if (visited_mask >> i) & 1 == 0 {
					// we are guaranteed to eventually hit this if statement
					let (i2, visited_mask2) = trace_chain(i, &masks[..]);
					// the masks are disjoint, so we don't need to check that, but we do need to check that their lengths add to 3
					return visited_mask.count_ones() + visited_mask2.count_ones() == 3 && visited_mask & (1 << i2) == 0 && visited_mask2 & (1 << i1) == 0;
					// case I (when i1 == i2) and case VII.1 (when i1 != i2)
				}
			}
			unreachable!("If we found 3 primes with 3 relations, each prime must have out degree 1, and if we traced a chain of length < 3 there must be an unvisted prime");
		} else if sqmask_pow == 2 { // cubefree cases
/*
III: n = 3p^2 where (p - 1) is divisible by 3
V.a: n = p^2qr where (q - 1) is divisible by p but not p^2, and (r - 1) is divisible by q
V.b: ^^^       where q = 3 and (p + 1) is divisible by 3
V.c: ^^^       where (r - 1) is divisible by p^2 and (p + 1) is divisible by q
VI: n = p^2q^2 where (q + 1) is divisible by p but not p^2
VII.2: n = q1r1p2^2q2 where (r1 - 1) is divisible by q1 and (p2 + 1) is divisible by q2
VII.3: n = q1^2p2q2r2 where (r2 - 1) is divisible by q2 and (q2 - 1) is divisible by p2
VII.4: n = q1^2p2^2q2 where (p2 + 1) is divisible by q2
*/
			if sqmask.count_ones() > 2 {
				return false;
			} else if sqmask.count_ones() == 2 { // we have two square divisors, so only VI and VII.4 are possible
				if g != 1 { // both of these cases only have weak edges p | (q + 1)
					return false;
				}
				let p = fxn.p[sqmask.trailing_zeros() as usize];
				let q = fxn.p[(15-sqmask.leading_zeros())as usize];
				// these masks are the converse of the normal masks: they have bits set for r when (p + 1) is divisible by r (rather than for p)
				let p_mask = weak_inedge_mask(&fxn, p);
				let q_mask = weak_inedge_mask(&fxn, q);
				if p_mask.count_ones() + q_mask.count_ones() > 1 {
					return false;
				}
				if q_mask != 0 {
					// case VII.4 || case VI
					return q_mask & sqmask == 0 || !(q + 1).is_multiple_of(p*p);
				} // otherwise, (q + 1) cannot be a multiple of p, so we just need to check that p_mask is not 0
				return p_mask != 0;
			} // otherwise, we have one square divisor, so we could be in case III, V.a, V.b, V.c, VII.2, or VII.3
			let vvv = gather_mults(&fxn, g, 2, 2);
			let p = fxn.p[sqmask.trailing_zeros() as usize];
			let p_mask = weak_inedge_mask(&fxn, p);
			if vvv.len() == 1 && vvv[0].1.is_power_of_two() && n.is_multiple_of(3) && (p-1).is_multiple_of(3) && p_mask == 0 {
				return true; // case III
			}
			if vvv.len() == 2 && vvv[0].0 == sqmask.trailing_zeros() as u8 && vvv[0].1 == 1 << vvv[1].0 && !(fxn.p[vvv[1].0 as usize] - 1).is_multiple_of(p*p) && p_mask == 0 {
				return true; // case V.a
			}
			if vvv.len() == 1 && fxn.p[vvv[0].0 as usize] == 3 && vvv[0].1.is_power_of_two() && vvv[0].1 != sqmask && p_mask.is_power_of_two() && (p + 1).is_multiple_of(3) {
				return true; // case V.b !!! TODO: this case appears to have a mistake since r is unused in the paper
			}
			if vvv.len() == 1 && vvv[0].0 == sqmask.trailing_zeros() as u8 && vvv[0].1.is_power_of_two() && (fxn.p[vvv[0].1.trailing_zeros()as usize]-1).is_multiple_of(p*p) && p_mask.is_power_of_two() && p_mask != vvv[0].1 {
				return true; // case V.c
			}
			if vvv.len() == 1 && vvv[0].1.is_power_of_two() && vvv[0].1 != sqmask && 1 << vvv[0].0 != sqmask && p_mask.is_power_of_two() && p_mask != vvv[0].1 && p_mask != 1 << vvv[0].0 {
				return true; // case VII.2
			}
			if vvv.len() == 2 && vvv[0].1 == 1 << vvv[1].0 && p_mask == 0 && vvv[0].1 != sqmask && vvv[1].1 != sqmask && 1 << vvv[0].0 != sqmask {
				return true; // case VII.3
			}
			return false;
		} else { // cubeful cases
/*
IV: n = p^3q where (p^3-1) is divisible by q
recall that (p^2-1) should not be divisible by q though, so we could "simplify" to q divides p^2 + p + 1 but not p + 1 or p - 1.
*/
			let p = fxn.p[sqmask.trailing_zeros() as usize];
			if sqmask.count_ones() != 1 || g != 1 || weak_inedge_mask(&fxn, p) != 0 {
				return false;
			}
			return fxn.p[..fxn.num as usize].iter().copied().filter(|&q|(p*p+p+1).is_multiple_of(q)).count() == 1;
		}
	}
}

/// Check if gnu(n) = 7
/// Per Mahmoud 2024, there are 9 cases where this holds.
/// In the following table, p, q, r, s are distinct odd primes, and when we write "n=blah", we mean "n=blah * m", where m is squarefree and has no
/// relations (including extended relations or cube relations) with "blah" or itself.  Additionally, there are NO relations (including extended relations
/// or cube relations) other than those described.
/// And if explicit primes like 3 or 5 appear, then any named primes like p, q must be distinct from those.
/// * I: n = 5pq where 5|p-1, 5|q-1
/// * II: n = 3qr where 3|q-1, 3|r-1, and q|r-1
/// * III: n = 3pqr where 3|p-1, 3|q-1, and q|r-1
/// * IV: n = 5p^2 where 5|p-1
/// * V: n = p^3q where q|p+1
/// * VI: n = 5p^2q where 5|q-1 and 5|p+1
/// * VII: n = p^2qr where p^2|q-1 and q|r-1
/// * VIII: n = p^2q^2 where p^2|q+1
/// * IX: n = p^2qrs where q|p+1, r|p+1, and p||s-1
pub fn is_gnu_7(n: u64) -> bool {
	unsafe {
		let mut fxn = MaybeUninit::uninit();
		n_factor_init(fxn.as_mut_ptr());
		let mut fxn = fxn.assume_init();
		n_factor(&mut fxn, n, 0);
		let mut sqmask = 0u16;
		let mut sqmask_pow = 0;
		for (i, &e) in fxn.exp[..fxn.num as usize].iter().enumerate() {
			if e > 3 {
				return false;
			}
			if e > 1 {
				if sqmask_pow == 0 {
					sqmask_pow = e;
				} else if sqmask_pow != e {
					return false;
				}
				if sqmask.count_ones() == 2 {
					return false;
				}
				sqmask |= 1 << i;
			}
		}
		let g = gcd(red_phi(&fxn), n);
		if sqmask_pow == 0 { // sqfree cases
/*
I: n = 5pq where (p - 1) is divisible by 5 and (q - 1) is divisible by 5 (NB: THIS IS THE FIRST INSTANCE WHERE A NODE HAS OUT DEGREE > 1)
II: n = 3qr where (q - 1) is divisible by 3, (r - 1) is divisible by 3, and (r - 1) is divisible by q
III: n = 3pqr where (p - 1) is divisible by 3, (q - 1) is divisible by 3, and (r - 1) is divisible by q
*/
			let vvv = gather_mults(&fxn, g, 2, 3);
			if vvv.len() == 0 {
				return false;
			} else if vvv.len() == 1 {
				return fxn.p[vvv[0].0 as usize] == 5 && vvv[0].1.count_ones() == 2;
				// case I
			} // otherwise there are 2 primes that have 2-3 total relations, so we can only have case II or III
			if fxn.p[0] != 3 {
				return false; // these cases can't happen if n is a multiple of 2, and certainly can't happen if n is NOT a multiple of 3
			}
			if vvv[0].0 != 0 || vvv[0].1.count_ones() != 2 {
				return false; // again, both cases require 3 to have out-edges to two other primes
			}
			// now, the last edge is from "q" to "r", where "q" MUST be an out-neighbor of 3, but "r" is unrestricted
			return (1 << vvv[1].0) & vvv[0].1 != 0; // is the other prime with an out edge an out neighbor of 3?
		} else if sqmask_pow == 2 { // cubefree cases
/*
IV: n = 5p^2 where (p - 1) is divisible by 5
VI: n = 5p^2q where (q - 1) is divisible by 5 and (p + 1) is divisble by 5
VII: n = p^2qr where (q - 1) is a multiple of p^2 and (r - 1) is a multiple of q
VIII: n = p^2q^2 where (q + 1) is a multiple of p^2
IX: n = p^2qrs where (p + 1) is a multiple of q, (p + 1) is a multiple of r, and (s - 1) is a multiple of p but not p^2
*/
			let p = fxn.p[sqmask.trailing_zeros() as usize];
			if sqmask.count_ones() > 2 {
				return false;
			} else if sqmask.count_ones() == 2 { // we can only be in case VIII
				let q = fxn.p[(15 - sqmask.leading_zeros()) as usize];
				return g == 1 && weak_inedge_mask(&fxn, p) == 0 && weak_inedge_mask(&fxn, q).is_power_of_two() && (q + 1).is_multiple_of(p*p);
			} // otherwise p is the only square factor and we can only be in case IV, VI, VII, or IX
			let vvv = gather_mults(&fxn, g, 2, 2);
			if vvv.len() == 0 {
				return false;
			} else if vvv.len() == 2 { // can only be case VII
				return vvv[0].1 == 1 << vvv[1].0 && (fxn.p[vvv[1].0 as usize] - 1).is_multiple_of(p*p) && weak_inedge_mask(&fxn, p) == 0;
			} // otherwise, we have one prime with out edges and we can only be in case IV, VI, or IX
			if !vvv[0].1.is_power_of_two() {
				return false;
			}
			if 1 << vvv[0].0 != sqmask && fxn.p[vvv[0].0 as usize] == 5 {
				if vvv[0].1 == sqmask {
					if weak_inedge_mask(&fxn, p) == 0 {
						return true; // case IV
					}
				} else {
					if weak_inedge_mask(&fxn, p) == 1 << vvv[0].0 {
						return true; // case VI
					}
				}
			}
			// case IX does not preclude p from being 5, which is why we don't return false within the above if statement if cases IV/VI fail
			if 1 << vvv[0].0 == sqmask && !(fxn.p[vvv[0].1.trailing_zeros() as usize] - 1).is_multiple_of(p*p) {
				let p_mask = weak_inedge_mask(&fxn, p);
				if p_mask.count_ones() == 2 && vvv[0].1 & p_mask == 0 {
					return true;
				}
			}
			return false;
		} // cubeful cases
/*
V: n = p^3q where (p + 1) is divisible by q
there should be no relations (g == 1), (p + 1) should only be divisible by q, and I think we need to check that no prime
divides (p^2 + p + 1) (the paper does not mention this directly since q can't divide this since it is coprime to p + 1, and the rest
of the prime divisors are not mentioned since they are "arithmetically independent")
*/
		if !sqmask.is_power_of_two() || g != 1 {
			return false;
		}
		let p = fxn.p[sqmask.trailing_zeros() as usize];
		return weak_inedge_mask(&fxn, p).is_power_of_two() && fxn.p[..fxn.num as usize].iter().all(|&r|!(p*p + p + 1).is_multiple_of(r));
	}
}

