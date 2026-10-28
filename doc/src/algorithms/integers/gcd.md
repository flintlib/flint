---
title: Integer GCD
summary: Greatest common divisors of multiprecision integers by the binary algorithm, Lehmer's algorithm and the subquadratic half-GCD, in the matrix formulation shared with continued fractions; contents and least common multiples of vectors by divisibility tests; the port of GMP's implementation with transform-space matrix products behind fmpz_gcd and fmpz_xgcd, and FLINT's own half-GCD for continued fractions and rational reconstruction.
complexity: Euclid and Lehmer $O(n^2)$; half-GCD $O(M(n) \log n)$ for $n$-limb operands.
order: 30
---

The greatest common divisor of two integers is computed by the
Euclidean algorithm or one of its accelerations. For $n$-limb operands
the Euclidean algorithm takes $O(n \bits)$ division steps of decreasing
size and $O(n^2)$ limb operations in total; Lehmer's algorithm keeps this
complexity but replaces almost all multiprecision steps by word
operations; and the half-GCD algorithm of Knuth and Schönhage reduces the
cost to $O(M(n) \log n)$. The same algorithms compute continued
fractions, Bézout coefficients and rational reconstructions, which are
all byproducts of the sequence of quotients. The polynomial analogue is
described in
[](algorithms/polynomials/half-gcd).

## The Euclidean algorithm as a product of matrices {#matrices}

Starting from $r_0 = a > r_1 = b > 0$, the Euclidean algorithm computes
$r_{i+1} = r_{i-1} - q_i r_i$ with $q_i = \lfloor r_{i-1}/r_i \rfloor$
until $r_{k+1} = 0$, and then $r_k = \gcd(a, b)$. In matrix form,
$$
\begin{pmatrix} r_{i-1} \\ r_i \end{pmatrix}
= \begin{pmatrix} q_i & 1 \\ 1 & 0 \end{pmatrix}
\begin{pmatrix} r_i \\ r_{i+1} \end{pmatrix},
\qquad
\begin{pmatrix} a \\ b \end{pmatrix}
= M_i \begin{pmatrix} r_i \\ r_{i+1} \end{pmatrix},
\quad M_i = \prod_{j \le i} \begin{pmatrix} q_j & 1 \\ 1 & 0 \end{pmatrix}.
$$
The matrix $M_i$ has nonnegative entries and determinant $(-1)^i$, so
its inverse is again an integer matrix, and its entries are the
numerators and denominators of the convergents of the continued fraction
$a/b = [q_1; q_2, \ldots, q_k]$. The second row of $M_i^{-1}$ gives the
cofactors $s_i, t_i$ with $s_i a + t_i b = r_i$; at $i = k$ these are
Bézout coefficients. Each accelerated algorithm computes a product of
many quotient matrices at once from leading parts of the operands, and
applies it to the full operands.

The quotients that leading parts determine are characterised by a
statement about continued fractions of real numbers. For a matrix $M$
as above, write $M \cdot y = (m_{11} y + m_{12})/(m_{21} y + m_{22})$.

::: lemma {#lem-cfrac}
Let $M = \prod_{j=1}^{i} \left(\begin{smallmatrix} q_j & 1 \\ 1 & 0
\end{smallmatrix}\right)$ with integers $q_j \ge 1$. For every real
$y > 1$, the continued fraction expansion of $x = M \cdot y$ begins with
$q_1, \ldots, q_i$. In particular, all real numbers in $M \cdot Y$,
for a set $Y$ of reals greater than 1, have these first $i$ partial
quotients.
:::

::: proof
By induction on $i$. For $i = 1$, $x = q_1 + 1/y$ with $0 < 1/y < 1$, so
$\lfloor x \rfloor = q_1$. For $i > 1$, $x = q_1 + 1/x'$ with $x' = M'
\cdot y$, where $M'$ is the product of the last $i - 1$ matrices. By
induction $x'$ begins with $q_2 \ge 1$, and $x' > q_2 \ge 1$ since the
rest of its expansion is positive; hence $\lfloor x \rfloor = q_1$ and the
expansion of $x$ continues with that of $x'$.
:::

::: corollary {#cor-truncation}
Let $a > b > 0$ and $s \ge 0$, and write $a = a_h 2^s + a_l$ and
$b = b_h 2^s + b_l$ with $0 \le a_l, b_l < 2^s$. Every partial quotient
that is common to all real numbers in the open interval
$(a_h/(b_h + 1), (a_h + 1)/b_h)$ is a partial quotient of $a/b$.
:::

::: proof
$a/b > a_h 2^s / ((b_h + 1) 2^s)$ and $a/b < (a_h + 1) 2^s / (b_h 2^s)$,
so $a/b$ lies in the interval.
:::

Thus quotients can be computed from the top bits of $a$ and $b$, as long
as the computation checks that the quotients found are valid for the whole
interval; by [](#lem-cfrac), it suffices to check that the transformed
endpoints stay above 1.

## Binary GCD {#binary}

Stein's binary algorithm [@Ste1967] avoids divisions altogether: it
removes common factors of 2, then repeatedly replaces the larger of two
odd numbers by their difference with its factors of 2 removed. For
single words this is the fastest method on current processors, since a
subtraction, a count of trailing zeros and a shift are cheap. Stehlé and
Zimmermann [@SZ2004] gave a subquadratic "binary recursive" variant that
works from the least significant bits; it is not used in FLINT.

## Lehmer's algorithm {#lehmer}

Lehmer [@Leh1938] observed that the first quotients of $a/b$ are almost
always determined by the leading word of each operand. His algorithm
runs the Euclidean algorithm on the leading words $\hat a, \hat b$
(taken from the same bit position of $a$ and $b$), accumulating the
quotient matrix in single words, stops as soon as a quotient might
differ from that of the full numbers, and then applies the matrix to $a$
and $b$ with a few multiply-and-add passes. Each pass reduces the
operands by about one word at a cost of $O(n)$, so the total is $O(n^2)$
with a small constant.

Collins and Jebelean
[@Jeb1993b] showed that, with the remainders $\hat r_i$ and cofactors
$u_i, v_i$ of the single-word computation, the quotients $q_1, \ldots,
q_{i}$ are those of $a/b$ if
$$
\hat r_{i+1} \ge |v_{i+1}| \quad\text{and}\quad
\hat r_i - \hat r_{i+1} \ge |v_{i+1} - v_i|
$$
(with $u$ in place of $v$ for the other parity of $i$). This is the
interval condition of [](#cor-truncation) written in terms of the
cofactors. Using two words of each operand ("double-digit" Lehmer), as
GMP does in `mpn_hgcd2`, roughly doubles the progress per pass.

## The half-GCD algorithm {#hgcd}

The half-GCD of $a > b$ of $n$ limbs is a matrix $M$ such that $(a, b)^T
= M (a', b')^T$, where $a', b'$ are consecutive remainders with about
$n/2$ limbs; the entries of $M$ then have about $n/2$ limbs too. Knuth
[@Knu1970] and Schönhage [@Sch1971] observed that it can be computed
recursively from the top halves of $a$ and $b$.

::: algorithm Half-GCD {#alg-hgcd}
**Input:** $a > b > 0$ of $n$ limbs.

**Output:** a product $M$ of quotient matrices and $(a', b')$ with $(a,
b)^T = M (a', b')^T$ and $a'$ of about $n/2$ limbs.

1. If $n$ is below a threshold, run Lehmer's algorithm until the
   remainders have $n/2$ limbs and return.
2. Let $s$ be about $n/2$ limbs. Compute recursively the half-GCD $M_1$
   of the top parts $a_h, b_h$ (of about $n/2$ limbs), and set $(a_1,
   b_1)^T = M_1^{-1} (a, b)^T$; by [](#cor-truncation), after removing
   the last one or two quotients if necessary, $a_1, b_1$ are
   consecutive remainders of $a, b$ of about $3n/4$ limbs.
3. Do one Euclidean step on $(a_1, b_1)$.
4. Compute recursively the half-GCD $M_2$ of the top parts of the
   result, of about $n/2$ limbs, and apply it as in step 2.
5. Return $M = M_1 Q M_2$, with $Q$ the matrix of step 3, and the
   remainders.
:::

Each level does two recursive calls of half size and a constant number
of products of $n/2$-limb numbers and matrices, so the cost is $H(n) =
2H(n/2) + O(M(n))$, which is $O(M(n) \log n)$. The GCD is obtained by
applying the half-GCD repeatedly, each time halving the size, at the
same total cost. The correctness of the truncations in steps 2 and 4,
which is delicate in Schönhage's original formulation, was simplified by
Möller [@Mol2008], whose version is the one in GMP. FLINT's own version
formulates it with intervals as in [](#cor-truncation), which makes the
check of step 2 a comparison of the transformed endpoints with 1.

## Contents and least common multiples of vectors {#content}

The content of a vector $(x_1, \ldots, x_n)$, the GCD of its entries, can
be computed as a chain of $n - 1$ GCDs, but the chain is dominated by the
latency of the GCDs when the entries are small. The running GCD $g$ can
only decrease, and after a few entries it is usually the content.
Testing each entry for divisibility by $g$ first, and computing a GCD
only when the test fails, replaces most GCDs by divisibility tests: a
test fails only when $g$ loses a prime factor, so at most once per prime
factor of $g$, counted with multiplicity. The tests are independent of each other, so they run
at the throughput of the multiplier.

For word-size $x$ and $g = 2^k h$ with $h$ odd, $g \mid x$ if and only if
the low $k$ bits of $x$ are zero and $x h^{-1} \bmod B \le \lfloor (B -
1)/h \rfloor$, where $B$ is the word base: one multiplication and one
comparison. If the test fails, the high word $t$ of $(x h^{-1} \bmod B)
\cdot h$ satisfies $x \equiv -B t \pmod h$, so $\gcd(x, h) = \gcd(t, h)$
starts from an operand already reduced modulo $h$. An entry of several
words is reduced modulo a word-size $h$ with independent multiplications,
by Hensel remainders or by folding its words with precomputed powers $B^i
\bmod h$. For a multi-word $g$, which is kept odd, the entries are tested
by Hensel divisibility tests, or for long $g$ by division with
remainder.

The order of the entries matters when the content is 1 or the entries
are large. For a polynomial, the leading and trailing coefficients are
the most likely to differ from the others in their factors, so the
entries at both ends of the vector come first. Single-word entries come
before larger ones, so that $g$ usually fits in a word before the large
entries are reduced; among the large entries, the two shortest come
first. Once $g$ has settled, the remaining tests can be split between
threads.

The least common multiple uses the same idea. The running LCM $L$ can
only grow, and when the entries are denominators that are all 1, all
equal or all divisors of a common denominator, every entry divides $L$.
Each entry $d$ is therefore first tested for $d \mid L$, and $L$ is
updated only when the test fails. While $L$ fits in a word the tests need
no division, and afterwards the small entries are tested in batches with
one remainder computation of $L$ per batch.

Compared with chained GCDs and LCMs (`fmpz_gcd`, `fmpz_lcm`), for vectors
of 1000 entries on one core of a Cascade Lake Xeon:

| Entries | Content | Speedup | Entries | Speedup |
|---|---|---:|---|---:|
| 62 bits | random (content 1) | 1.0 | all 1 | 5.0 |
| 62 bits | 16-bit content | 57 | equal, 62 bits | 18 |
| 62 bits | content of half the bits | 38 | divisors of a 62-bit $D$ | 6.1 |
| 62 bits | 6-bit cofactors | 17 | $d_i = i$ | 6.6 |
| 128 bits | 16-bit content | 2.2 | $d_i = i!$ | 3.4 |
| 640 bits | content of half the bits | 2.7 | random, 62 bits | 3.2 |

The gains are small where every entry needs a multiprecision GCD or
division anyway: random entries with content 1, large entries with a
word-size content, and random multi-word denominators.

## Use in FLINT {#use}

**`fmpz_gcd`, `fmpz_gcd_ui`, `fmpz_gcd3`.** An inline operand reduces the
problem to GMP's `mpn_gcd_1`, which divides the large operand by the
small one once and finishes with a binary GCD on words. Two large
operands go to `flint_mpz_gcd`, that is `flint_mpn_gcd`, which sets aside
the common power of two and calls `_flint_mpn_gcd_inplace`, a port of
GMP's `mpn_gcd` [@GMP] that uses FLINT's multiplication and division:
Lehmer's algorithm with double-word steps (`_flint_mpn_hgcd2`) for small
operands and Möller's half-GCD for large ones. Unlike GMP, each Lehmer
step applies its word matrix to both operands in a single pass, and the
$2 \times 2$ matrix products of the half-GCD, whose entries share
operands, are computed for large operands in transform space with
`fft_small` ([](fft_small/transform-space)), so that each operand is
transformed once and each result is transformed back once.

**`fmpz_xgcd`, `fmpz_xgcd_canonical_bezout`, `fmpz_gcdinv`,
`fmpz_invmod`.** `flint_mpz_gcdext` and `flint_mpz_invert`, built on
`_flint_mpn_gcdext_inplace`, the port of GMP's `mpn_gcdext` on the same half-GCD
with cofactor tracking; unlike GMP, it accumulates the matrices of
consecutive half-GCD steps from right to left before applying them to the
cofactors, instead of updating the full-length cofactors after every
step. For word-size operands, GMP's
`mpn_gcdext_1` (`fmpz_xgcd_canonical_bezout`) or FLINT's `n_gcdinv`
(`fmpz_gcdinv`, `fmpz_invmod`) are used ([](fmpz/gcd)). The modular
inverses of `mpn_mod` call `_flint_mpn_gcdext_inplace` directly.

**`_fmpz_vec_content`, `_fmpz_vec_lcm` and related functions.**
`flint_mpn_gcd_vec_1` and `flint_mpn_gcd_vec` compute contents as in
[](#content), of single limbs and of multi-limb values given by pointer.
`_fmpz_vec_content`, `_fmpz_vec_content_chained`, `_fmpz_vec_lcm` and
`_fmpz_vec_lcm_strided` read the entries of an `fmpz` vector in place,
small entries in blocks and large entries by pointer to their limbs, and
use threads only for long vectors of large entries. They give the
contents and canonical forms of `fmpz_poly`, `fmpq_poly`, `fmpz_mat`,
`fmpz_mpoly`, `fmpq_mpoly` and `fmpz_mpoly_q`, and the common
denominators of `_fmpq_vec_get_fmpz_vec_fmpz` and the
`fmpq_mat_get_fmpz_mat` functions, whose numerators
`_fmpq_vec_get_fmpz_vec_scaled` forms without divisions for denominators
equal to the LCM or to 1.

**`n_gcd`, `n_xgcd`, `n_gcdinv`.** Word GCDs: `n_gcd` calls GMP's
binary `mpn_gcd_11`; `n_xgcd` and `n_gcdinv` are FLINT's own extended
Euclidean algorithms with special handling of small quotients.

**`fmpz_xgcd_partial`.** FLINT's own Lehmer algorithm with the
Jebelean condition of [](#lehmer), which stops when the remainder falls
below a bound $L$ and returns the last two remainders and cofactors.
It is the partial Euclidean algorithm of Shanks's NUCOMP and NUDUPL for
binary quadratic forms in `qfb`, and is also used by Cornacchia's
algorithm in `ecpp`.

**`_fmpq_hgcd` and `_fmpq_ball_get_cfrac`.** FLINT's own integer
half-GCD in the interval formulation of [](#lem-cfrac): it maintains a
matrix $M$ and a pair $(x_a, x_b)$ representing the interval $M \cdot
(x_a/(x_b + 1), (x_a + 1)/x_b)$, and generates quotients valid for the
whole interval. It uses single Euclidean steps when the remaining size is
below 4 limbs, a Lehmer algorithm with a two-limb half-GCD
(`_uiui_hgcd`) up to 500 limbs, and above that recursion on the top
halves as in [](#alg-hgcd). It computes continued fractions
(`fmpq_get_cfrac`), the simplest rational between two numbers
(`fmpq_simplest_between`), Dedekind sums (`fmpq_dedekind_sum`) and
rational reconstruction (`fmpq_reconstruct_fmpz_2`, which also has its
own Lehmer loop for the final steps).
[](algorithms/number-theory/continued-fractions) and
[](algorithms/number-theory/rational-reconstruction) describe these
applications.

::: functions fmpz_gcd fmpz_xgcd_canonical_bezout fmpz_xgcd_partial _fmpz_vec_content _fmpz_vec_lcm fmpq_get_cfrac fmpq_reconstruct_fmpz_2 :::
