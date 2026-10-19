---
title: Performance of dfloat
summary: Timings of the dfloat arithmetic, vector operations, elementary functions and complex functions, compared with nfloat and arb, and what they imply for choosing a format.
order: 3
---

The timings below measure the speed of the `dfloat` operations, to help
choose between the formats and between `dfloat`, `nfloat` and `arb_t`. All
numbers are nanoseconds per operation on one core of an Intel Xeon
(Skylake-AVX512 class, 2.8 GHz, a virtual machine), measured with the
programs in `src/dfloat/profile`. The absolute numbers vary by up to
1.5 times between runs on this machine, so only ratios within a table
are meaningful.

## Arithmetic {#arithmetic}

Scalar columns: calls on independent operands through the generic
ring. Vector columns: per element of a vector of length 1024
(`_gr_vec_add`, `_gr_vec_mul` and the dot product `_gr_vec_dot`), which
use SIMD four elements at a time. Suffixes: `s` a strong
(canonicalising) context, `f` a fast context (`DFLOAT_FAST`).
`nfloat` is at 64, 128, 192 and 256 bits and `arb_t` at $53N$ bits.

| Ring | add | mul | sqr | div | sqrt | vec_add | vec_mul | vec_dot |
|---|---|---|---|---|---|---|---|---|
| `d1` | 1.4 | 1.3 | 1.4 | 1.8 | 1.9 | 0.27 | 0.25 | 0.34 |
| `d1s` | 1.5 | 1.4 | 1.8 | 1.8 | 1.9 | 0.25 | 0.26 | 0.31 |
| `d1b` | 2.8 | 4.8 | 3.2 | 8.8 | 6.7 | 0.75 | 1.6 | 2.1 |
| `d1bf` | 2.4 | 4.3 | 3.5 | 8.6 | 6.9 | 0.77 | 1.4 | 1.9 |
| `nfloat` 64 | 6.4 | 4.5 | 3.5 | 23 | 21 | 3.2 | 1.6 | 2.0 |
| `arb` 53 | 24 | 20 | 17 | 55 | 76 | 21 | 19 | 8.6 |
| `d2` | 3.3 | 3.0 | 3.5 | 8.3 | 7.4 | 1.1 | 0.94 | 1.5 |
| `d2s` | 4.5 | 4.7 | 4.4 | 9.1 | 8.8 | 2.0 | 1.9 | 1.3 |
| `d2b` | 5.5 | 8.8 | 6.7 | 30 | 26 | 3.1 | 4.7 | 6.2 |
| `d2bf` | 6.4 | 9.3 | 7.2 | 32 | 29 | 3.1 | 4.1 | 5.5 |
| `nfloat` 128 | 6.7 | 5.5 | 3.9 | 44 | 88 | 3.8 | 2.7 | 4.8 |
| `arb` 106 | 36 | 22 | 20 | 80 | 114 | 34 | 21 | 11 |
| `d3` | 7.4 | 8.8 | 8.1 | 25 | 39 | 3.8 | 3.5 | 4.0 |
| `d3s` | 11 | 13 | 11 | 29 | 43 | 6.7 | 6.4 | 4.7 |
| `d3b` | 10 | 16 | 11 | 64 | 80 | 4.8 | 7.4 | 8.8 |
| `d3bf` | 12 | 20 | 14 | 76 | 84 | 4.8 | 7.6 | 8.3 |
| `nfloat` 192 | 7.2 | 7.2 | 5.0 | 84 | 118 | 4.4 | 6.7 | 8.1 |
| `arb` 159 | 41 | 33 | 32 | 80 | 130 | 40 | 32 | 25 |
| `d4` | 15 | 20 | 14 | 59 | 69 | 5.2 | 8.3 | 8.6 |
| `d4s` | 19 | 27 | 23 | 76 | 76 | 9.5 | 13 | 8.8 |
| `d4b` | 17 | 28 | 20 | 114 | 118 | 8.6 | 11 | 13 |
| `d4bf` | 17 | 30 | 20 | 114 | 118 | 8.8 | 11 | 13 |
| `nfloat` 256 | 7.2 | 8.3 | 6.4 | 99 | 118 | 5.2 | 8.3 | 11 |
| `arb` 212 | 41 | 37 | 35 | 99 | 153 | 41 | 35 | 26 |

Against `nfloat`, the plain expansions win clearly for $N \le 2$ (the
`d1` operations are single instructions), and the balls are comparable
for $N \le 2$ while carrying a rigorous error bound; the vector
operations and dot products are 1.5 to 3 times faster than `nfloat` for
$N \le 3$. For $N = 3, 4$, `nfloat` addition and multiplication (limb
arithmetic) are faster than the expansion operations, by a factor 2 for
products at $N = 4$. Division and square root of expansions are 1.5 to
4 times faster than `nfloat`. The balls cost 1.2 to 1.7 times the plain
numbers for products and 2 to 4 times for division and square root,
whose exact residuals cost a product and a difference.

Canonicalisation (`DFLOAT_STRONG`) costs nothing for $N = 1$, a check
(1 to 2 ns) for $N = 2$, and a dependent TwoSum chain for $N = 3, 4$:
about 5 and 8 ns per scalar operation (60% on a plain addition, 30% on
a ball addition, 25% on a ball product), and up to 100% on the cheapest
vector operations, whose results are canonicalised one at a time.

## Division and square roots {#division}

A Newton–Karp–Markstein division [@KM1997] (a reciprocal at half
precision, one product and one correction) is slower than the long
division at every $N$ (75 against 49 ns at $N = 4$ in one run): the
steps of the long division are cheap products of a single double with
the divisor, and it needs no full-precision square. Precomputing $1/y_0$
does not help either: the divisions are not on the critical path of the
products. Square roots, per element, in a run on a machine about twice
slower than in the table above (`b` denotes the ball versions):

| Format | sqrt | vec_sqrt | rsqrt | vec_rsqrt | b sqrt | b vec_sqrt | b rsqrt | b vec_rsqrt | vec_div |
|---|---|---|---|---|---|---|---|---|---|
| `d1` | 2.1 | 0.9 | 3.4 | 1.5 | 10 | 3.1 | 18 | 5.2 | 3.1 |
| `d2` | 12 | 3.7 | 10 | 3.1 | 48 | 16 | 57 | 18 | 6.4 |
| `d3` | 79 | 19 | 88 | 30 | 135 | 38 | 196 | 56 | 13 |
| `d4` | 122 | 34 | 150 | 45 | 227 | 53 | 312 | 89 | 28 |

## Elementary functions {#elementary}

Scalar calls through the generic ring and the vector functions
`_dN_vec_exp`, etc., per element of a vector of length 1024; the
`arb_t` vector versions are loops over the scalar functions.

| Ring | exp | vec_exp | sin | vec_sin | log | vec_log | atan | vec_atan |
|---|---|---|---|---|---|---|---|---|
| `d1` | 8.3 | 1.8 | 6.2 | 2.7 | 4.5 | 3.3 | 6.9 | 6.0 |
| `d1b` | 15 | 3.5 | 7.2 | 3.4 | 6.9 | 6.0 | 8.8 | 8.3 |
| `arb` 53 | 176 | 179 | 198 | 202 | 153 | 149 | 179 | 183 |
| `d2` | 80 | 14 | 76 | 25 | 59 | 23 | 74 | 22 |
| `d2b` | 88 | 15 | 80 | 26 | 58 | 22 | 76 | 24 |
| `arb` 106 | 214 | 210 | 244 | 248 | 210 | 210 | 248 | 236 |
| `d3` | 153 | 32 | 187 | 56 | 160 | 48 | 183 | 53 |
| `d3b` | 156 | 34 | 168 | 58 | 172 | 49 | 176 | 50 |
| `arb` 159 | 278 | 305 | 305 | 305 | 252 | 248 | 320 | 305 |
| `d4` | 305 | 62 | 320 | 111 | 305 | 92 | 366 | 99 |
| `d4b` | 305 | 58 | 320 | 111 | 305 | 92 | 382 | 103 |
| `arb` 212 | 336 | 351 | 412 | 382 | 305 | 305 | 382 | 382 |

The ball versions cost about the same as the plain ones for $N \ge 2$
(the radius is a static relative bound plus the propagated input
radius) and up to twice as much for $N = 1$, where the plain kernels are
a few instructions. For $N = 1$ the kernels are about as fast as the C
library on the same machine (`exp` 4.6, `sin` 6.7, `log` 3.3, `atan`
4.8 ns), whose near-correct rounding they do not attempt. Against
`arb_t`, the scalar functions are 20 to 40 times faster for $N = 1$,
3 times for $N = 2$ and about as fast for $N = 4$; the vector functions
50 to 100, 10 and 4 to 6 times. `sin_cos` costs about 10% more than
`sin`; `expm1` and `log1p` cost the same as `exp` and `log`, and the
compositions (`tan`, `pow`, ...) the sum of their parts.

## Polynomial multiplication {#polynomials}

Classical multiplication of two polynomials of length 256, which runs
on `_dN_vec_dot_rev`, costs per coefficient product 0.54, 3.4, 9.0 and
18 ns for `d1` to `d4`, and 3.5, 9.2, 16 and 25 ns for `d1b` to `d4b`
(against 1.0, 12, 25, 61 and 8.4, 9.7, 16, 27 ns with a scalar loop);
the gains start at length 16.

## Complex numbers {#complex}

Calls through the generic ring, against `acb_t` at $53N$ bits:

| Ring | mul | div | sqrt | exp | log | sin | tan | pow |
|---|---|---|---|---|---|---|---|---|
| `d1cb` | 35 | 118 | 113 | 41 | 99 | 94 | 162 | 182 |
| `acb` 53 | 86 | 360 | 512 | 652 | 766 | 1042 | 1362 | 1055 |
| `d2cb` | 43 | 148 | 170 | 202 | 270 | 313 | 457 | 544 |
| `acb` 106 | 103 | 386 | 447 | 492 | 612 | 771 | 1329 | 1287 |
| `d3cb` | 90 | 303 | 448 | 413 | 620 | 650 | 920 | 1174 |
| `acb` 159 | 120 | 576 | 670 | 687 | 724 | 1016 | 1433 | 1619 |
| `d4cb` | 119 | 412 | 478 | 695 | 971 | 960 | 1232 | 1617 |
| `acb` 212 | 104 | 418 | 529 | 706 | 694 | 952 | 1546 | 1581 |

The complex balls are several times faster than `acb_t` for $N \le 2$,
faster for $N = 3$ and on par for $N = 4$. Dot products of length 64,
per term (`loop`: a loop of complex `gr_mul` and `gr_add`):

| Format | real | complex | complex loop | real ball | complex ball | complex ball loop |
|---|---|---|---|---|---|---|
| `d1` | 0.36 | 1.06 | 3.99 | 1.69 | 7.8 | 23.5 |
| `d2` | 1.83 | 8.30 | 32.1 | 4.86 | 21.6 | 65.0 |
| `d3` | 5.06 | 22.1 | 69.7 | 9.88 | 36.7 | 95.7 |
| `d4` | 10.0 | 42.8 | 144 | 16.5 | 61.7 | 197 |

## Benchmarking notes {#benchmarking}

On Skylake-family processors, whose JCC erratum microcode keeps jumps
that cross or end on a 32-byte boundary out of the decoded-instruction
cache, the speed of short hot loops (dot products, short vector
operations) varies by up to 25% with the placement of the code, between
builds of identical instructions. Comparisons of such loops are best
made with the module assembled with
`-Wa,-mbranches-within-32B-boundaries`.

Everything that depends on $N$ is compiled four times and the kernels
are fully unrolled, so the object code is large (about 0.9 MB of text,
most of it the scalar and vector cores of the elementary functions and
of division and square root). Small operations are inlined only where
that pays off, and calls within the module never go through the
procedure linkage table; relative to inlining everything, this saves
36% of the text without measurable loss of speed.
