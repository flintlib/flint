---
title: Matrix arithmetic in multi-double formats
summary: How dfloat multiplies matrices of double-word to quad-word numbers with a SIMD kernel that accumulates the order sums of the products without renormalising each term, uses BLAS for a single double, and builds LU decompositions and triangular solving on the product.
complexity: $mkn$ terms for an $m \times k$ by $k \times n$ product, each $O(N^2)$ double operations (about 140 for $N = 4$); the block recursive LU decomposition and triangular solving cost a constant times the product.
order: 86
---

The plain formats of `dfloat` ([](dfloat)) implement the matrix methods
of the generic rings, so that `gr_mat_mul`, `gr_mat_lu`,
`gr_mat_nonsingular_solve`, `gr_mat_inv` and the eigensolvers of
`gr_mat` run on them with the generic block algorithms and fast base
cases. This page describes the product kernel and the base cases; the
expansions themselves and their error-free transformations are
described in [](algorithms/floating-point/double-word).

## One double {#d1}

For $N = 1$ the product is `flint_dgemm`, an optimised BLAS-like double
precision product. A complex product $C = AB$ is a single real product
of twice the size: the $m \times k$ complex matrix $A$, read as an
$m \times 2k$ real matrix of parts, times the $2k \times 2n$ real matrix
$$\begin{pmatrix} \operatorname{Re} B & \operatorname{Im} B \\ -\operatorname{Im} B & \operatorname{Re} B \end{pmatrix}$$
(with interleaved rows and columns, a copy of $B$), whose result is the
$m \times 2n$ matrix of the parts of $C$.

## The product kernel {#kernel}

For $N \ge 2$ the product is a SIMD kernel. The columns of $B$ are
packed, a block of $W$ columns at a time ($W = 4$ with AVX2 or NEON and
8 with AVX-512), in structure-of-arrays form: component $c$ of the
entries $b_{l,j_0}, \ldots, b_{l,j_0+W-1}$ is one vector register. Each
row of $A$ then runs along a block, its entries broadcast to the lanes,
so that the $W$ lanes compute $W$ entries of $C$ in parallel.

The terms $a_{il} b_{lj}$ of a dot product are not renormalised one by
one. The product of two expansions $a = \sum_i a_i$ and
$b = \sum_j b_j$ splits into order sums $T_0, \ldots, T_{N-1}$, where
$T_r$ collects the products $a_i b_j$ with $i + j = r$ and the error
terms of the products of the order below (as in the scalar
multiplication kernel; exact except for the plain additions of the last
order, which are relative to the term). These go into accumulators
$S_0, \ldots, S_N$ of the order sums of the whole dot product: each
addition into $S_r$ is a TwoSum whose error term goes on to $S_{r+1}$,
and only the additions into the extra accumulator $S_N$ round. The
accumulators are renormalised every 32 terms, before their magnitudes
can drift apart, and at the end. For `d4` this costs about 140 double
operations per term, a third fewer than a renormalised product followed
by a renormalised sum, and it is more accurate than a dot product of
rounded products: the error is that of one rounding of the exact dot
product plus the terms dropped from each product. For $k = 64$ it was
measured at most $0.4 \cdot 2^{-53N} \sum_l |a_{il}| |b_{lj}|$ for
entries of uniform magnitude and $3.2 \cdot 2^{-53N} \sum_l |a_{il}|
|b_{lj}|$ for magnitudes spread over $2^{\pm 20}$, against 1.0 and 5.1
for the classical algorithm with the scalar operations. (Without the
extra accumulator $S_N$ the kernel lost five to eight bits.)

Only exact operations and fused multiply-subtracts (for the error terms
of products) are used, so the results do not depend on the vector
width. A complex term is four real terms; with real entries the zero
imaginary parts cost nothing extra, since the additions of an exact
zero are exact. Two rows of $A$ are processed together, sharing the
loads of $B$. The rows of $A$ are split into panels small enough for the
panel and a packed block of $B$ to stay in the L2 cache, and the panels
run in parallel with the threads FLINT is set to use. An entry whose
result is not finite is recomputed with the scalar operations.

The kernel is used from dimension 2. Its cost is cubic, while the block
algorithms of `arb_mat_approx_mul` and `acb_mat_approx_mul` (exact
integer matrix products through BLAS) are cheaper per term for large
matrices; from dimension 1600, 1000 and 600 (`d2`, `d3`, `d4`; complex:
1500, 900, 500) the entries are converted to `arf_t` and multiplied at
$53N + 4$ bits, with an accuracy in the norm of the blocks rather than
entrywise. This makes `dfloat` depend on `arb_mat` and `acb_mat`, which
it otherwise sits below; it is a provisional solution until a fast
arbitrary-precision floating-point matrix product is available at a
lower level.

## LU decomposition and triangular solving {#lu}

`gr_mat_lu` is the block recursive algorithm (`gr_mat_lu_recursive`)
above 16 rows and columns, which spends most of its time in matrix
products $A_{11} - A_{10} A_{01}$. It leaves to the base case the
factorization of tall panels of a few columns, where the generic
classical algorithm, working along the short rows, pays a function call
per entry. The base case has the same pivoting (the entry of largest
magnitude, judged by the head of the expansion) and the same layout of
$L$ as `gr_mat_lu_classical`, but eliminates column by column, four rows
at a time in the SIMD lanes: the multipliers $e_j = a_{j,c}/p$ for the
rows below the pivot $p$, then $x_{j,c'} \gets x_{j,c'} - e_j\, p_{c'}$
with $p_{c'}$ broadcast. The results are those of the scalar operations.

The triangular solvers are the block recursive algorithms above 16. In
the base case, for four or more right-hand sides the rows of the
solution are computed as vectors, $X_j = (B_j - \sum_i L_{ji} X_i) /
L_{jj}$, with the SIMD products of a vector by a scalar; for fewer, by
the classical algorithm with dot products along the rows of $L$.

The recursion subtracts each product from the block it updates in a
separate pass, for want of a fused operation $C \gets C - AB$ in
`gr_mat`; such an operation, as a ring method tuned for each type,
would save a temporary matrix and a pass over it.

## Accuracy {#accuracy}

The decompositions are those of Gaussian elimination with partial
pivoting in the arithmetic of the format, backward stable in the usual
sense with the unit roundoff replaced by the effective accuracy of the
format, a small multiple of $2^{-53N}$. Measured on random linear
systems, the solutions are as accurate as with `arf_t` at $53N$ bits to
within a bit; see [](algorithms/ball-arithmetic/linear-systems#backends)
for how the approximate functions of `arb_mat` use this.

## Use in FLINT {#use}

- The methods `GR_METHOD_MAT_MUL`, `GR_METHOD_MAT_LU`,
  `GR_METHOD_MAT_NONSINGULAR_SOLVE_TRIL` and
  `GR_METHOD_MAT_NONSINGULAR_SOLVE_TRIU` of the plain real and complex
  `dfloat` rings; timings in [](dfloat/performance#linear-algebra).
- The approximate functions of `arb_mat` and `acb_mat` (LU
  decomposition, solving, inversion, eigenvalues) convert to these
  rings when they are faster.
