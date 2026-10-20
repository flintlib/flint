---
title: Integer GCD
summary: Greatest common divisors of multiprecision integers by the binary algorithm, Lehmer's algorithm and the subquadratic half-GCD, in the matrix formulation shared with continued fractions; GMP's implementations behind fmpz_gcd and fmpz_xgcd, and FLINT's own half-GCD for continued fractions and rational reconstruction.
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

## Use in FLINT {#use}

**`fmpz_gcd`, `fmpz_gcd_ui`, `fmpz_gcd3`.** An inline operand reduces the
problem to GMP's `mpn_gcd_1`, which divides the large operand by the
small one once and finishes with a binary GCD on words. Two large
operands go to GMP's `mpz_gcd`, that is `mpn_gcd`: Lehmer's algorithm
with double-word steps (`mpn_hgcd2`) below GMP's threshold
`GCD_DC_THRESHOLD` (several hundred limbs) and Möller's half-GCD above
[@GMP].

**`fmpz_xgcd`, `fmpz_xgcd_canonical_bezout`, `fmpz_gcdinv`,
`fmpz_invmod`.** GMP's `mpz_gcdext` and `mpz_invert`, built on the same
`hgcd` machinery with cofactor tracking; for word-size operands, GMP's
`mpn_gcdext_1` (`fmpz_xgcd_canonical_bezout`) or FLINT's `n_gcdinv`
(`fmpz_gcdinv`, `fmpz_invmod`) ([](fmpz/gcd)).

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

::: functions fmpz_gcd fmpz_xgcd_canonical_bezout fmpz_xgcd_partial fmpq_get_cfrac fmpq_reconstruct_fmpz_2 :::
