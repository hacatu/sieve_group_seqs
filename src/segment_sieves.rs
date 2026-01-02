/*
Core sieve implementations to find candidates for length 2+, 4+, and 6+ sequences
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



We can prove that k must be 8 mod 16 and 0 mod 9, so k = 144m + 72.
This is explained in detail at https://math.stackexchange.com/questions/5062059, but the short version is
- Since gnu(k + 1) = 1, k + 1 = 2 or k + 1 is odd.  But if k + 1 = 2, then k + 2 = 3, and so gnu(k + 2) would be 1, not 2.  So 2 | k.
- Suppose k = 4m + 2.  Then gnu(k + 2) = gnu(4(m+1)).  But if 8 | 4(m+1), then gnu(4(m+1)) >= gnu(8) = 5, and if 8 does not divide 4(m+1),
  then m+1 is odd and gnu(4(m+1)) >= 4*gnu(m+1) >= 4.  But we need gnu(k+2) to be 2, so we must have 4 | k, not k = 4m + 2.
- Suppose k = 8m + 4.  Then 8 | k + 4, so gnu(k + 4) >= 5, which is too big, so we must have 8 | k.
- Now, we want to show 3 | k, and we know 8 | k, so we will prove that k = 24m + 8 and k = 24m + 16 do not work.  If k = 24m + 8, then k + 4 = 12(2m + 1),
  so gnu(k + 4) >= gnu(12) = 5, but this is too big.  If k = 24m + 16, then k+2 = 6(4m + 3).  There must be some number of the form 6p that divides 6(4m + 3),
  since 4m + 3 cannot be 1.  (Also, 4m + 3 is never even, and we can check it can't be 9, so even if it is a power of 3 that would not work).  But
  gnu(6(4m + 3)) >= gnu(6p), and we can list at least 4 non-isomorphic groups of size 6p (C2 x C3 x Cp, C3 x D2p, Cp x D6, D6p),
  so gnu(k+2) >= 4 which is too big.
- If k is a multiple of 16, then k = 16m, and k + 4 = 4(4m + 1).  Refer to the doc comment on `is_gnu_4`: there are 6 cases, but only 3 of them can apply to
  non-squarefree numbers.  Of these, one requires no relations, but 2 is related to any odd prime.  Another requires one relation which must be between two
  non-squared primes, but if we have two non-squared primes then 2 would be related to both of them (more than one relation).  And the last one requires
  that there is one relation (only possible if 4m + 1 is prime) and that the squared prime p^2 not divide any q-1, but if 4m+1 is prime then 2^2 divides 4m+1-1.
  So we have proven 16 does not divide k.
- The last condition, that 9 divides k, requires gnu(k+6)=6.  If we assume that k = 144m+24 or 144m+120, then k+6 is 6(24m+5) or 18(8m+7).
  But per Mahmoud 2024, if k+6 is even and has gnu 6, only case II is possible, which requires k+6 = 2pq.  For 18(8m+7) this is impossible because 18 isn't
  squarefree.  For 6(24m+5), this is possible iff 24m+5 is prime, but then there is NOT a relation between 3 and 24m+5, and we need there to be one.
  So k cannot be 3 or 6 mod 9, it must be 0.

And altogether we have k = 144m + 72.

Then if we look at k+2, k+4, and k+6, and analyze what cases in the above papers are actually possible,
we find that
- 144m + 73 must be squarefree
- 72m + 37 must be prime
- 36m + 19 must be prime
- 24m + 13 must be prime
and for any m where this holds, the gnu(k+2) = 2, gnu(k+4) = 4, and gnu(k+6) = 6 condtions must hold,
so we just have to check the gnu(k+odd) condtions.

For length 2-3 and 4-5 sequences, we can find similar conditions, but the modulus is less and less constrained.
*/

use std::cmp::max;

/// Compute the inverse of 144 mod p, by combining the inverses of 2 and 3 mod p which are easy to compute.
/// 2^-1 = (p + 1)/2, and 3^-1 is (2p + 1)/3 if p is 1 mod 3, or (p + 1)/3 if p is 2 mod 3.
/// This is only used for sieving primes p>3, so 144 is guaranteed to be invertible.
fn inv144(p: u32) -> u32 {
	let p = p as u64;
	let inv2 = (p + 1)/2;
	let inv3 = if p%3 == 1 { 2*p + 1 } else { p + 1 } / 3;
	((inv2.pow(2)%p*inv3%p).pow(2)%p) as _
}

/// Lift the inverse of 144 mod p to the inverse mod p^2, using hensel's lifting lemma and avoiding 128 bit math.
/// This takes 3 32-bit widening multiplications and 2 64 bit divisions.
/// The first division could be done with 32 bits for n_max up to 889600953739896, but it would only take about 81 minutes
/// to exceed this bound on my computer, so we definitely want to support larger numbers.
fn lift_inv144(p: u32, x: u32) -> u64 {
	let p = p as u64;
	let x = x as u64;
	let q = (144*x - 1)/p; // must be divisible because x is the inverse of 144 mod p
	// hensel lifting picks which value of k makes 144*(x + k*p) - 1 = 0 mod p^2:
	// by definition, 144*x - 1 = p*q, so we have p*q + 144*k*p = 0 mod p^2
	// and we can divide through by p now to get q + 144*k = 0 mod p, and solve for k
	let k = (p - x)*q%p;
	x + k*p
}

/// Core function for the segmented sieve: find all "candidate" n corresponding to m in the interval [m_a, m_b).
/// This is just the inner loop for `find_67_seqs`, which manages efficiently splitting up the whole interval and running this function in multiple threads.
/// n = 144*m + 72, and it is a candidate if
/// * 144*m + 73 is squarefree
/// * 72*m + 37 is prime
/// * 36*m + 19 is prime
/// * 24*m + 13 is prime
/// and for such an n, we just have to check that the ODD gnu(n+i)=i are correct; the even constraints are GUARANTEED.
/// It is the caller's responsibility (ie, the function `find_67_seqs`) to
/// * compute the sieving primes
/// * split the target range into segments and call this function on them
/// * allocate a bitarray `bucket` of a reasonable size
/// * consume candidate n values from `out` and filter them to find the ones that also satisfy the odd constraints
pub fn sieve_segment_6(sieving_primes: &[u32], m_a: u64, m_b: u64, bucket: &mut [u8], out: &mut Vec<u64>) {
	let nbits = (m_b - m_a)as usize;
	bucket[..nbits.div_ceil(8)].fill(0);
	let mut fallthrough = 4;
	for &p in sieving_primes {
		let p = p as u64;
		/*
		first, mark off all m where 144m + 73 is a multiple of p^2
		*/
		if p*p > 144*(m_b-1) + 73 {
			break
		}
		let inv144_p = inv144(p as _) as u64;
		let inv144_p2 = lift_inv144(p as _, inv144_p as _);
		let r = p*p - inv144_p2 * 73 % (p*p);
		/*
		m_a <= p**2 * q + r
		(m_a - r + p**2 - 1)//p**2 <= q
		p**2 * q + r <= m_b - 1
		q <= (m_b - r - 1)//p**2
		*/
		let q_a = (m_a + p*p - r - 1)/(p*p);
		let q_b = (m_b + p*p - r - 1)/(p*p);
		for q in q_a..q_b {
			let idx = (p*p*q + r - m_a)as usize;
			bucket[idx/8] |= 1 << (idx%8);
		}
		if fallthrough == 1 {
			continue;
		}
		/*
		now, we have to mark off all m where 72m + 37 is a multiple of p
		in addition to bounds like those for p**2 vs 144m + 73, we need to make sure
		that 72m + 37 >= p**2
		m >= (p**2 - 37 + 72 - 1)//72
		*/
		if p*p > 72*(m_b-1) + 37 {
			fallthrough = 1;
			continue;
		}
		let inv72_p = inv144_p*2%p;
		let mut r = inv72_p * (p-37%p) % p;
		let tmp = (p*p + 72 - 37 - 1)/72;
		let q_a = (max(m_a, tmp) + p - r - 1)/p;
		let q_b = (m_b + p - r - 1)/p;
		for q in q_a..q_b {
			let idx = (p*q + r - m_a)as usize;
			bucket[idx/8] |= 1 << (idx%8);
		}
		if fallthrough == 2 {
			continue;
		}
		/*
		now, we have to mark off all m where 36m + 19 is a multiple of p
		*/
		if p*p > 36*(m_b-1) + 19 {
			fallthrough = 2;
			continue;
		}
		if r < inv72_p {
			r = r + p - inv72_p;
		} else {
			r -= inv72_p;
		}
		// now r is inv72_p * -38 % p
		let tmp = (p*p + 36 - 19 - 1)/36;
		let q_a = (max(m_a, tmp) + p - r - 1)/p;
		let q_b = (m_b + p - r - 1)/p;
		for q in q_a..q_b {
			let idx = (p*q + r - m_a)as usize;
			bucket[idx/8] |= 1 << (idx%8);
		}
		if fallthrough == 3 {
			continue;
		}
		/*
		finally, we have to mark off all m where 24m + 13 is a multiple of p
		*/
		if p*p > 24*(m_b-1) + 13 {
			fallthrough = 3;
			continue;
		}
		if r < inv72_p {
			r = r + p - inv72_p;
		} else {
			r -= inv72_p;
		}
		// now r is inv72_p * -39 % p
		let tmp = (p*p + 24 - 13 - 1)/24;
		let q_a = (max(m_a, tmp) + p - r - 1)/p;
		let q_b = (m_b - r + p - 1)/p;
		for q in q_a..q_b {
			let idx = (p*q + r - m_a)as usize;
			bucket[idx/8] |= 1 << (idx%8);
		}
	}
	for m in m_a..m_b {
		let idx = (m - m_a)as usize;
		if bucket[idx/8] & (1 << (idx%8)) == 0 {
			out.push(144*m + 72);
		}
	}
}

/// Compute the inverse of 48 mod p, by combining the inverses of 2 and 3 mod p which are easy to compute.
/// 2^-1 = (p + 1)/2, and 3^-1 is (2p + 1)/3 if p is 1 mod 3, or (p + 1)/3 if p is 2 mod 3.
/// This is only used for sieving primes p>3, so 48 is guaranteed to be invertible.
fn inv48(p: u32) -> u32 {
	let p = p as u64;
	let inv2 = (p + 1)/2;
	let inv3 = if p%3 == 1 { 2*p + 1 } else { p + 1 } / 3;
	((inv2.pow(2)%p).pow(2)%p*inv3%p) as _
}

/// Lift the inverse of 48 mod p to the inverse mod p^2, using hensel's lifting lemma and avoiding 128 bit math.
/// This takes 3 32-bit widening multiplications and 2 64 bit divisions.
/// The first division could be done with 32 bits for n_max up to about 2.7 x 10^15, so we need 64 bits to support large n.
fn lift_inv48(p: u32, x: u32) -> u64 {
	let p = p as u64;
	let x = x as u64;
	let q = (48*x - 1)/p; // must be divisible because x is the inverse of 48 mod p
	// hensel lifting picks which value of k makes 48*(x + k*p) - 1 = 0 mod p^2:
	// by definition, 48*x - 1 = p*q, so we have p*q + 48*k*p = 0 mod p^2
	// and we can divide through by p now to get q + 48*k = 0 mod p, and solve for k
	let k = (p - x)*q%p;
	x + k*p
}

/// Core function for 4-5 length sequence sieve
/// n = 48*m + 24 is a candidate if
/// * 48*m + 25 is squarefree
/// * 24*m + 13 is prime
/// * 12*m + 7 is prime
pub fn sieve_segment_4(sieving_primes: &[u32], m_a: u64, m_b: u64, bucket: &mut [u8], out: &mut Vec<u64>) {
	let nbits = (m_b - m_a)as usize;
	bucket[..nbits.div_ceil(8)].fill(0);
	let mut fallthrough = 3;
	for &p in sieving_primes {
		let p = p as u64;
		/*
		first, mark off all m where 48m + 25 is a multiple of p^2
		*/
		if p*p > 48*(m_b-1) + 25 {
			break
		}
		let inv48_p = inv48(p as _) as u64;
		let inv48_p2 = lift_inv48(p as _, inv48_p as _);
		let r = if p == 5 { 0 } else { p*p - inv48_p2 * 25 % (p*p) };
		/*
		m_a <= p**2 * q + r
		(m_a - r + p**2 - 1)//p**2 <= q
		p**2 * q + r <= m_b - 1
		q <= (m_b - r - 1)//p**2
		*/
		let q_a = (m_a + p*p - r - 1)/(p*p);
		let q_b = (m_b + p*p - r - 1)/(p*p);
		for q in q_a..q_b {
			let idx = (p*p*q + r - m_a)as usize;
			bucket[idx/8] |= 1 << (idx%8);
		}
		if fallthrough == 1 {
			continue;
		}
		/*
		now, we have to mark off all m where 24m + 13 is a multiple of p
		in addition to bounds like those for p**2 vs 48m + 25, we need to make sure
		that 24m + 13 >= p**2
		m >= (p**2 - 13 + 24 - 1)//24
		*/
		if p*p > 24*(m_b-1) + 13 {
			fallthrough = 1;
			continue;
		}
		let inv24_p = inv48_p*2%p;
		let mut r = inv24_p * (p-13%p) % p;
		let tmp = (p*p + 24 - 13 - 1)/24;
		let q_a = (max(m_a, tmp) + p - r - 1)/p;
		let q_b = (m_b + p - r - 1)/p;
		for q in q_a..q_b {
			let idx = (p*q + r - m_a)as usize;
			bucket[idx/8] |= 1 << (idx%8);
		}
		if fallthrough == 2 {
			continue;
		}
		/*
		now, we have to mark off all m where 12m + 7 is a multiple of p
		*/
		if p*p > 12*(m_b-1) + 7 {
			fallthrough = 2;
			continue;
		}
		if r < inv24_p {
			r = r + p - inv24_p;
		} else {
			r -= inv24_p;
		}
		// now r is inv24_p * -14 % p
		let tmp = (p*p + 12 - 7 - 1)/12;
		let q_a = (max(m_a, tmp) + p - r - 1)/p;
		let q_b = (m_b + p - r - 1)/p;
		for q in q_a..q_b {
			let idx = (p*q + r - m_a)as usize;
			bucket[idx/8] |= 1 << (idx%8);
		}
	}
	for m in m_a..m_b {
		let idx = (m - m_a)as usize;
		if bucket[idx/8] & (1 << (idx%8)) == 0 {
			out.push(48*m + 24);
		}
	}
}

/// Core function for 2-3 length sequence sieve
/// n = 4*m is a candidate if
/// * 4*m + 1 is squarefree
/// * 2*m + 1 is prime
pub fn sieve_segment_2(sieving_primes: &[u32], m_a: u64, m_b: u64, bucket: &mut [u8], out: &mut Vec<u64>) {
	let nbits = (m_b - m_a)as usize;
	bucket[..nbits.div_ceil(8)].fill(0);
	let mut fallthrough = 2;
	for &p in sieving_primes {
		let p = p as u64;
		/*
		first, mark off all m where 4m + 1 is a multiple of p^2
		*/
		if p*p > 4*(m_b-1) + 1 {
			break
		}
		let inv4_p2 = (3*p*p + 1)/4;
		let r = p*p - inv4_p2;
		/*
		m_a <= p**2 * q + r
		(m_a - r + p**2 - 1)//p**2 <= q
		p**2 * q + r <= m_b - 1
		q <= (m_b - r - 1)//p**2
		*/
		let q_a = (m_a + p*p - r - 1)/(p*p);
		let q_b = (m_b + p*p - r - 1)/(p*p);
		for q in q_a..q_b {
			let idx = (p*p*q + r - m_a)as usize;
			bucket[idx/8] |= 1 << (idx%8);
		}
		if fallthrough == 1 {
			continue;
		}
		/*
		now, we have to mark off all m where 2m + 1 is a multiple of p
		in addition to bounds like those for p**2 vs 4m + 1, we need to make sure
		that 2m + 1 >= p**2
		m >= (p**2 - 1 + 2 - 1)//2
		*/
		if p*p > 2*(m_b-1) + 1 {
			fallthrough = 1;
			continue;
		}
		let inv2_p = (p + 1)/2;
		let r = inv2_p * (p-1) % p;
		let tmp = (p*p + 2 - 1 - 1)/2;
		let q_a = (max(m_a, tmp) + p - r - 1)/p;
		let q_b = (m_b + p - r - 1)/p;
		for q in q_a..q_b {
			let idx = (p*q + r - m_a)as usize;
			bucket[idx/8] |= 1 << (idx%8);
		}
	}
	for m in m_a..m_b {
		let idx = (m - m_a)as usize;
		if bucket[idx/8] & (1 << (idx%8)) == 0 {
			out.push(4*m + 1);
		}
	}
}

