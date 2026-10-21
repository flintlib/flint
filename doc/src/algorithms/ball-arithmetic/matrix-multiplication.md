---
title: Ball matrix multiplication
summary: How arb_mat_mul and acb_mat_mul multiply ball matrices by splitting the inner dimension into blocks of entries with similar exponents, multiplying the midpoints exactly as integer matrices, and bounding the radius product in double precision, with the resulting error bounds.
complexity: For $n \times n$ matrices at $p$ bits: one exact integer matrix product of $O(p)$-bit entries per block (typically one block), using the fast algorithms of fmpz_mat, plus $O(n^3)$ double operations for the radii, against $O(n^3)$ multiplications at $p$ bits for the classical algorithm.
order: 80
---

The product $C = AB$ of an $M \times N$ matrix $A$ and an $N \times P$
matrix $B$ with ball entries $A = [A_m \pm A_r]$, $B = [B_m \pm B_r]$
must satisfy $\tilde A \tilde B \in C$ for all matrices $\tilde A \in A$,
$\tilde B \in B$. Since
$$\tilde A \tilde B - A_m B_m = (\tilde A - A_m) \tilde B + A_m (\tilde B - B_m),$$
the entrywise bound
$$|\tilde A \tilde B - A_m B_m| \le A_r (|B_m| + B_r) + |A_m| B_r$$
holds (with absolute values and inequalities taken entrywise), so it is
enough to compute an enclosure of the midpoint product $A_m B_m$ and an
upper bound for the radius product on the right. The classical algorithm
computes each entry as a dot product in ball arithmetic (`arb_dot`, see
[](algorithms/ball-arithmetic/dot-products)), which does exactly this
entry by entry with $O(N)$ multiplications of $p$-bit numbers. For large
matrices, both parts can instead be done with fast matrix
multiplication, as observed by Rump for midpoint-radius interval
matrices in floating-point arithmetic [@Rum1999]; FLINT does this in
`arb_mat_mul_block`, with exact integer products for the midpoints.

## Midpoints: exact integer matrix products {#midpoints}

Every midpoint is a binary floating-point number $a = u \cdot 2^e$ with
an odd integer $u$ of at most $p$ bits. If all entries in row $i$ of
$A_m$ are multiples of $2^{s_i}$ and have absolute value less than
$2^{s_i + h}$, then $2^{-s_i}$ times the row is a row of integers of at
most $h$ bits, and similarly for the columns of $B_m$ with exponents
$t_j$. The scaled product of integer matrices is then exact, and
$(A_m B_m)_{ij} = 2^{s_i + t_j} (\hat A \hat B)_{ij}$.

::: algorithm Block multiplication of midpoints {#alg-block}
Input: $A_m$, $B_m$, precision $p$. Let $b$ be the largest number of bits
of any midpoint and $H = 1.25 \min(p, b) + 192$.

1. Split the inner index range $0 \le k < N$ into consecutive blocks
   $K$, each as long as possible subject to the condition that in each
   row $i$ of $A_m$ (restricted to the columns in $K$) and in each
   column $j$ of $B_m$ (restricted to the rows in $K$), the nonzero
   entries span at most $H$ bits, from the lowest set bit of any of them
   to the highest.
2. For a block with at least 30 indices, convert the block to integer
   matrices $\hat A$, $\hat B$ with the row and column scalings above,
   compute $\hat A \hat B$ exactly with `fmpz_mat_mul`, round each entry
   of $2^{s_i + t_j} (\hat A \hat B)_{ij}$ to $p$ bits and add it to
   $C_{ij}$ in ball arithmetic.
3. For a shorter block (at least 30 indices are taken anyway), add the
   dot products of the midpoints over the block with `arb_dot`.
:::

For matrices whose entries in each row and column have similar
magnitudes, which is the common case, there is a single block, and the
computed midpoint of each entry has a single rounding error of at most
half an ulp at $p$ bits: the result is as accurate as the classical
algorithm, which also rounds each dot product once. With $k$ blocks,
each entry accumulates at most $k$ roundings. Entries of very different
magnitude in the same row (for example, a matrix with a huge first
column) create several blocks, and a matrix with huge dynamic range
degenerates to classical multiplication in blocks of 30. The bound $H$
of about $1.25p$ bits allows the integer entries to be somewhat longer
than the precision, so that ordinary dynamic range within a row does not
split blocks; the constant 192 accommodates low-precision inputs.

The integer matrix multiplication `fmpz_mat_mul` chooses its algorithm
from the dimension and the entry size ([](fmpz_mat)). At precisions of a
few hundred bits, its multimodular method, which replaces the $p$-bit
products by matrix products modulo word-size primes, gives most of the
speedup over the classical algorithm.

## Radii: products in double precision {#radii}

The radius term $A_r(|B_m| + B_r) + |A_m| B_r$ needs only an upper bound
of moderate accuracy, so it is computed with matrices of nonnegative
`double` numbers. The magnitudes $|A_m| + A_r$, $|B_m|$, $A_r$, $B_r$
are converted to `mag_t` upper bounds (30-bit mantissas), and
`_arb_mat_addmul_rad_mag_fast` again splits the inner dimension into
blocks in which the exponents of each row and column span at most 900
bits. In a block, each row and each column is scaled by a power of two
so that all entries are representable as normal `double` numbers without
rounding, and the product is computed with ordinary `double` arithmetic
in cache-friendly blocks.

::: lemma {#lem-double}
Let $x_k, y_k \ge 0$ be `double` numbers and let $s$ be the result of
summing the products $x_k y_k$, $k = 1, \ldots, n$, in any order with
rounding to nearest, without overflow or underflow. Then $s (1 + 2.01(n
+ 1) \cdot 2^{-52}) \ge \sum_k x_k y_k$.
:::

::: proof
With unit roundoff $u = 2^{-53}$, each computed product and each partial
sum of nonnegative terms has relative error at most $u$, and the
standard error analysis of summation gives $\sum_k x_k y_k \le s (1 -
u)^{-n}$, since each term passes through at most $n$ roundings (one
multiplication and at most $n - 1$ additions) [@Hig2002]. For $nu <
0.01$ (larger $n$ would need matrices with about $10^{14}$ columns),
$(1 - u)^{-n} \le 1 + nu/(1 - nu) \le 1 + 1.02 nu$, which is well below
the stated factor; the margin also covers the rounding of the final
multiplication by the factor.
:::

The scaled `double` products are multiplied by this factor and added,
with rounding upwards, to the radii of $C$; blocks of at most 5 indices
use `mag_t` arithmetic directly. The radius computation costs $O(MNP)$
floating-point operations in the worst case, but these are much cheaper
than $p$-bit multiplications, and for exact inputs (all radii zero) it is
skipped. If only one factor is exact, only one radius product is
computed.

::: theorem {#thm-matmul}
The matrix $C$ computed by `arb_mat_mul_block` contains $\tilde A \tilde
B$ for all $\tilde A \in A$ and $\tilde B \in B$.
:::

::: proof
The midpoint part encloses $A_m B_m$: in each block, the integer product
is exact and the rounding to $p$ bits and the addition are performed in
ball arithmetic, which adds the rounding errors to the radius. The
radius part adds an upper bound for $A_r(|B_m| + B_r) + |A_m| B_r$ by
[](#lem-double) and the upward rounding of the conversions between
`mag_t` and `double`, which bounds $|\tilde A \tilde B - A_m B_m|$ by the
identity at the beginning of this section.
:::

## Algorithm selection {#selection}

`arb_mat_mul` uses the classical algorithm (threaded if several threads
are available and the matrices are not tiny) when one of the dimensions
is at most 60 (for $p \le 128$), 50 (for $p \le 512$) or 40, and
`arb_mat_mul_block` otherwise. The block algorithm itself falls back to
classical multiplication if some midpoint has more than 8000 bits, if
both matrices have less than 10% nonzero entries and $H > 1024$, or if
some exponent is too large for the fast integer exponent arithmetic.

For complex matrices, `acb_mat_mul_reorder` writes $A = A_1 + i A_2$
and $B = B_1 + i B_2$ and computes the real and imaginary parts $A_1 B_1
- A_2 B_2$ and $A_1 B_2 + A_2 B_1$ with four real block
multiplications (two if one of the matrices is real, one if both are).
`acb_mat_mul` uses it when the inner dimension $N$ and the number of
columns $P$ of $B$ satisfy $\min(N, P) \ge \max(20, 5 + b/64)$, where
$b < 8000$ is the minimum of $p$ and the largest number of bits of a
midpoint; otherwise it multiplies classically (threaded if possible).

## Use in FLINT {#use}

- `arb_mat_mul`, `arb_mat_mul_classical`, `arb_mat_mul_threaded`,
  `arb_mat_mul_block`, `arb_mat_sqr`, `_arb_mat_addmul_rad_mag_fast`
  ([](arb_mat)); `acb_mat_mul`, `acb_mat_mul_reorder` ([](acb_mat)).
- `gr_mat_mul` over the contexts `gr_ctx_init_real_arb` and
  `gr_ctx_init_complex_acb` calls `arb_mat_mul` and `acb_mat_mul`
  ([](gr_mat)).
- The same idea of exact integer products of blocks of similar
  magnitude is used for polynomials in
  [](algorithms/ball-arithmetic/polynomial-multiplication).
- Preconditioned linear solving ([](algorithms/ball-arithmetic/linear-systems))
  and the eigenvalue enclosures
  ([](algorithms/ball-arithmetic/eigenvalues)) spend most of their time
  in these products.
