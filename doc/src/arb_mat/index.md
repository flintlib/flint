---
title: arb_mat.h
header: arb_mat.h
description: Dense matrices with real ball entries: rigorous arithmetic, verified solving of linear systems and inversion, determinants, Cholesky and LDL decompositions, the matrix exponential, and approximate floating-point linear algebra for preconditioning.
category: Real and complex numbers
order: 9
---

This module implements dense matrices whose entries are real balls
(`arb_t`). A matrix of balls represents the set of all real matrices
whose entries lie in the balls, and the functions return enclosures of
the exact results for every matrix in that set. Besides the usual
arithmetic, the module provides verified linear algebra: solving linear
systems, inversion and determinants by Gaussian elimination in ball
arithmetic or by preconditioning with an approximate inverse, Cholesky
and $LDL^T$ decompositions of symmetric positive definite matrices, the
characteristic polynomial and the matrix exponential. The functions
with `approx` in their names ignore the radii and compute
floating-point approximations without error bounds; they serve to
compute preconditioners and starting values for the verified
algorithms.

Large products are computed by splitting the matrices into blocks of
similarly scaled entries and multiplying the midpoints as integer
matrices with `fmpz_mat_mul`, while the radii are bounded by a product
of nonnegative matrices in double precision
([](algorithms/ball-arithmetic/matrix-multiplication)). The solvers
follow [](algorithms/ball-arithmetic/linear-systems). Complex matrices,
including eigenvalues and eigenvectors, are in [](acb_mat). Several
functions are wrappers around the generic implementations in
[](gr_mat) with the real ball context, and a matrix over the generic
real ball context uses this module for multiplication, determinants and
the exponential. The guide [](arb_mat/getting-started) introduces the
module, and [](arb_mat/solve) describes the solvers in detail.

::: summary algorithms/ball-arithmetic/linear-systems :::

## Types and conventions {#conventions}

::: type
arb_mat_struct
arb_mat_t

A matrix with real ball entries. The structure has the fields `entries`
(a pointer to the entry in the top left corner), `r` and `c` (the
numbers of rows and columns) and `stride` (the distance in `arb_struct`
between the starts of consecutive rows). An `arb_mat_t` is an array of
one `arb_mat_struct`, so it is passed by reference.
:::

::: macro
arb_mat_entry(mat, i, j)

Evaluates to a pointer to the entry in row $i$ and column $j$ of `mat`
(no bounds check).
:::

::: function
arb_ptr arb_mat_entry_ptr(arb_mat_t mat, slong i, slong j)

Returns a pointer to the entry in row $i$ and column $j$ of `mat`; the
same as `arb_mat_entry`, as a function.
:::

::: macro
arb_mat_nrows(mat)
arb_mat_ncols(mat)

Evaluate to the number of rows and columns of `mat`.
:::

- The dimensions of a matrix are fixed when it is initialised, and the
  caller must pass operands of compatible dimensions. Multiplication,
  solving and the decompositions throw an exception (`FLINT_ERROR`)
  otherwise; entrywise functions such as `arb_mat_add` do not check.
  Matrices may have zero rows or columns.
- Outputs may be the same object as inputs (for square matrices where
  this makes sense) unless stated otherwise; the functions use
  temporary storage when needed.
- Functions that can fail numerically (`arb_mat_solve`, `arb_mat_inv`,
  `arb_mat_lu`, `arb_mat_cho`, ...) return an `int`: nonzero means that
  the result is certified (for example, that the matrix is certainly
  invertible and the output contains the exact solution); zero means
  that the matrix is singular, or that the input or the working
  precision is not accurate enough to decide. The output is then
  undefined.
- Predicates such as `arb_mat_eq` return 1 if the property certainly
  holds for all matrices represented, and 0 otherwise.
- The functions with `approx` in their name use only the midpoints and
  return matrices of floating-point numbers (zero radii, unless stated
  otherwise) without any error bound.
- Functions described as wrappers of a `gr_mat` function call it with
  `gr_ctx_init_real_arb(ctx, prec)` (see [](gr_domains)) on the same
  data.

## Memory management {#memory}

::: function
void arb_mat_init(arb_mat_t mat, slong r, slong c)

Initialises `mat` as the zero matrix with `r` rows and `c` columns.
:::

::: function
void arb_mat_clear(arb_mat_t mat)

Frees the memory used by `mat`.
:::

::: function
void arb_mat_swap(arb_mat_t mat1, arb_mat_t mat2)
void arb_mat_swap_entrywise(arb_mat_t mat1, arb_mat_t mat2)

Swap two matrices. `arb_mat_swap` swaps the structures (also when the
dimensions differ); `arb_mat_swap_entrywise` swaps the entries one by
one and requires equal dimensions, which also works for windows.
:::

::: function
slong arb_mat_allocated_bytes(const arb_mat_t x)

Returns the number of bytes allocated on the heap by `x`, excluding
`sizeof(arb_mat_struct)`.
:::

::: function
void arb_mat_window_init(arb_mat_t window, const arb_mat_t mat, slong r1, slong c1, slong r2, slong c2)
void arb_mat_window_clear(arb_mat_t window)

Initialise `window` as a view of the submatrix of `mat` formed by rows
$r_1 \le i < r_2$ and columns $c_1 \le j < c_2$, sharing its entries,
and release it again (which frees nothing). Require $0 \le r_1 \le r_2
\le$ the number of rows and likewise for the columns. A window must not
be used after `mat` has been cleared.
:::

## Conversions {#conversions}

::: function
void arb_mat_set(arb_mat_t dest, const arb_mat_t src)
void arb_mat_set_fmpz_mat(arb_mat_t dest, const fmpz_mat_t src)
void arb_mat_set_round_fmpz_mat(arb_mat_t dest, const fmpz_mat_t src, slong prec)
void arb_mat_set_fmpq_mat(arb_mat_t dest, const fmpq_mat_t src, slong prec)

Set `dest` to `src`, which must have the same dimensions: exactly, or
with the entries rounded to `prec` bits.
:::

## Random generation {#random}

::: function
void arb_mat_randtest(arb_mat_t mat, flint_rand_t state, slong prec, slong mag_bits)

Sets `mat` to a random matrix of random density whose nonzero entries
are generated by `arb_randtest` or `arb_randtest_precise` with
midpoints of up to `prec` bits and exponents of up to `mag_bits` bits.
:::

::: function
void arb_mat_randtest_cho(arb_mat_t mat, flint_rand_t state, slong prec, slong mag_bits)
void arb_mat_randtest_spd(arb_mat_t mat, flint_rand_t state, slong prec, slong mag_bits)

Set the square matrix `mat` to a random lower triangular matrix with
exact entries and a positive diagonal (a Cholesky factor), or to the
symmetric positive definite matrix $L L^T$ for such a random $L$.
:::

## Input and output {#io}

::: function
void arb_mat_printd(const arb_mat_t mat, slong digits)
void arb_mat_fprintd(FILE * file, const arb_mat_t mat, slong digits)

Print `mat` to standard output or to `file`, one row per line as `[a,
b, ...]`, each entry printed with `arb_fprintd` to `digits` decimal
digits.
:::

## Comparisons {#comparisons}

::: function
int arb_mat_equal(const arb_mat_t mat1, const arb_mat_t mat2)

Returns whether the matrices have the same dimensions and identical
balls as entries.
:::

::: function
int arb_mat_overlaps(const arb_mat_t mat1, const arb_mat_t mat2)

Returns whether the matrices have the same dimensions and each entry of
`mat1` overlaps the corresponding entry of `mat2`.
:::

::: function
int arb_mat_contains(const arb_mat_t mat1, const arb_mat_t mat2)
int arb_mat_contains_fmpz_mat(const arb_mat_t mat1, const fmpz_mat_t mat2)
int arb_mat_contains_fmpq_mat(const arb_mat_t mat1, const fmpq_mat_t mat2)

Return whether the matrices have the same dimensions and each entry of
`mat2` is contained in the corresponding entry of `mat1`.
:::

::: function
int arb_mat_eq(const arb_mat_t mat1, const arb_mat_t mat2)
int arb_mat_ne(const arb_mat_t mat1, const arb_mat_t mat2)

Return whether `mat1` and `mat2` certainly represent the same matrix
(same dimensions and all entries exact and equal), respectively
certainly represent different matrices (different dimensions, or some
pair of entries that do not overlap).
:::

::: function
int arb_mat_is_empty(const arb_mat_t mat)
int arb_mat_is_square(const arb_mat_t mat)

Return whether `mat` has zero rows or zero columns, respectively as
many rows as columns.
:::

::: function
int arb_mat_is_exact(const arb_mat_t mat)
int arb_mat_is_zero(const arb_mat_t mat)
int arb_mat_is_finite(const arb_mat_t mat)

Return whether all entries of `mat` have zero radius, are exactly zero,
or are finite.
:::

::: function
int arb_mat_is_triu(const arb_mat_t mat)
int arb_mat_is_tril(const arb_mat_t mat)
int arb_mat_is_diag(const arb_mat_t mat)

Return whether all entries below, above, respectively off the main
diagonal of `mat` are exactly zero. The matrix need not be square.
:::

## Special matrices {#special-matrices}

::: function
void arb_mat_zero(arb_mat_t mat)
void arb_mat_one(arb_mat_t mat)
void arb_mat_ones(arb_mat_t mat)
void arb_mat_indeterminate(arb_mat_t mat)

Set all entries of `mat` to zero; the diagonal entries to one and the
others to zero; all entries to one; or all entries to indeterminate
balls $[\operatorname{NaN} \pm \infty]$.
:::

::: function
void arb_mat_hilbert(arb_mat_t mat, slong prec)

Sets `mat` to the Hilbert matrix with entries $1/(j + k + 1)$, rounded
to `prec` bits. A wrapper of `gr_mat_hilbert`.
:::

::: function
void arb_mat_pascal(arb_mat_t mat, int triangular, slong prec)

Sets `mat` to a Pascal matrix of binomial coefficients: for `triangular`
0 the symmetric matrix with entries $\binom{j+k}{j}$ (the rows of
Pascal's triangle as antidiagonals), for 1 the upper triangular matrix
with entries $\binom{k}{j}$, for $-1$ the lower triangular matrix with
entries $\binom{j}{k}$. A wrapper of `gr_mat_pascal`; the entries are
computed by recurrences and are exact as long as they fit in `prec`
bits.
:::

::: function
void arb_mat_stirling(arb_mat_t mat, int kind, slong prec)

Sets `mat` to the matrix of Stirling numbers $s(j, k)$: unsigned of the
first kind for `kind` 0, signed of the first kind for 1, of the second
kind for 2. A wrapper of `gr_mat_stirling`; the entries are computed by
recurrences.
:::

::: function
void arb_mat_dct(arb_mat_t mat, int type, slong prec)

Sets `mat` to the orthogonal DCT-II matrix of order $n$, the smaller
dimension of `mat`,
$$A_{j,k} = \sqrt{\frac{2 - \delta_{j,0}}{n}} \cos\left(\frac{\pi j}{n} \left(k + \frac{1}{2}\right)\right),$$
which satisfies $A^{-1} = A^T$ (if `mat` is not square, the formula is
used for all entries). `type` is ignored and should be 0.
:::

## Transpose {#transpose}

::: function
void arb_mat_transpose(arb_mat_t dest, const arb_mat_t src)

Sets `dest` to the transpose of `src`, exactly. Aliasing is allowed for
square matrices.
:::

## Norms {#norms}

::: function
void arb_mat_bound_inf_norm(mag_t b, const arb_mat_t A)

Sets `b` to an upper bound for the infinity norm $\max_i \sum_j
|a_{ij}|$ of every matrix in `A`.
:::

::: function
void arb_mat_frobenius_norm(arb_t res, const arb_mat_t A, slong prec)
void arb_mat_bound_frobenius_norm(mag_t res, const arb_mat_t A)

Set `res` to an enclosure, respectively an upper bound, of the
Frobenius norm $(\sum_{i,j} a_{ij}^2)^{1/2}$.
:::

## Arithmetic {#arithmetic}

::: function
void arb_mat_neg(arb_mat_t dest, const arb_mat_t src)

Sets `dest` to $-\mathit{src}$, exactly.
:::

::: function
void arb_mat_add(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, slong prec)
void arb_mat_sub(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, slong prec)

Set `res` to the sum or difference of matrices of equal dimensions.
:::

::: function
void arb_mat_mul_classical(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, slong prec)
void arb_mat_mul_threaded(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, slong prec)
void arb_mat_mul_block(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, slong prec)
void arb_mat_mul(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, slong prec)

Set `C` (`res`) to the matrix product $AB$. The *classical* version
computes each entry with `arb_dot`; the *threaded* version does the
same with the rows of the output distributed over the threads of
`flint_get_num_threads()`; the *block* version multiplies the midpoints
exactly as integer matrices and bounds the radii separately. The
default chooses automatically.

::: details
`arb_mat_mul` uses the block version if all three dimensions exceed a
cutoff (60 for $\mathit{prec} \le 128$, 50 for $\mathit{prec} \le 512$,
40 otherwise), and the classical (or, with several threads and more
than $10^5$ bit operations per thread unit, the threaded) version
otherwise.

The *block* version groups consecutive columns of $A$ and rows of $B$
into blocks in which the entries of each row of $A$ (each column of $B$)
have exponents within about $1.25 \min(\mathit{prec}, b) + 192$ bits,
where $b$ is the largest midpoint size; converts each block to an
integer matrix scaled per row (per column); multiplies the integer
matrices with `fmpz_mat_mul` (multimodular or BLAS-like algorithms,
see [](algorithms/linear-algebra/integer-matrix-multiplication)); and
rounds and accumulates the results. Blocks with fewer than 30 columns
are multiplied by dot products instead. The radii are bounded by
$|A_m| B_r + A_r(|B_m| + B_r)$, computed as products of nonnegative
matrices with `_arb_mat_addmul_rad_mag_fast`. Matrices with entries
whose exponents are not small, or whose midpoints have more than 8000
bits, or that are very sparse at high precision, are multiplied
classically ([](algorithms/ball-arithmetic/matrix-multiplication)).
:::
:::

::: function
void arb_mat_mul_entrywise(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, slong prec)

Sets `C` to the entrywise (Hadamard) product of matrices of equal
dimensions.
:::

::: function
void arb_mat_sqr_classical(arb_mat_t B, const arb_mat_t A, slong prec)
void arb_mat_sqr(arb_mat_t res, const arb_mat_t mat, slong prec)

Set the output to $A^2$ for a square matrix $A$. The classical version
uses the symmetry of the products $a_{ik} a_{kj}$ to save
multiplications.
:::

::: function
void arb_mat_pow_ui(arb_mat_t res, const arb_mat_t mat, ulong exp, slong prec)

Sets `res` to $A^e$ for a square matrix $A$, by binary exponentiation.
A wrapper of `gr_mat_pow_ui`.
:::

::: function
void arb_mat_approx_mul(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, slong prec)

Approximate matrix multiplication: the input radii are ignored and the
output is a floating-point approximation of the product of the
midpoints. The radii of the output are not necessarily zero (they are
whatever the algorithm used leaves), and carry no meaning.
:::

## Scalar arithmetic {#scalar}

::: function
void arb_mat_scalar_mul_2exp_si(arb_mat_t B, const arb_mat_t A, slong c)

Sets `B` to $2^c A$, exactly.
:::

::: function
void arb_mat_scalar_addmul_si(arb_mat_t B, const arb_mat_t A, slong c, slong prec)
void arb_mat_scalar_addmul_fmpz(arb_mat_t B, const arb_mat_t A, const fmpz_t c, slong prec)
void arb_mat_scalar_addmul_arb(arb_mat_t B, const arb_mat_t A, const arb_t c, slong prec)

Set `B` to $B + cA$.
:::

::: function
void arb_mat_scalar_mul_si(arb_mat_t B, const arb_mat_t A, slong c, slong prec)
void arb_mat_scalar_mul_fmpz(arb_mat_t B, const arb_mat_t A, const fmpz_t c, slong prec)
void arb_mat_scalar_mul_arb(arb_mat_t B, const arb_mat_t A, const arb_t c, slong prec)

Set `B` to $cA$.
:::

::: function
void arb_mat_scalar_div_si(arb_mat_t B, const arb_mat_t A, slong c, slong prec)
void arb_mat_scalar_div_fmpz(arb_mat_t B, const arb_mat_t A, const fmpz_t c, slong prec)
void arb_mat_scalar_div_arb(arb_mat_t B, const arb_mat_t A, const arb_t c, slong prec)

Set `B` to $A / c$.
:::

## Matrix-vector products {#vector}

::: function
void _arb_mat_vector_mul_row(arb_ptr res, arb_srcptr v, const arb_mat_t A, slong prec)
void _arb_mat_vector_mul_col(arb_ptr res, const arb_mat_t A, arb_srcptr v, slong prec)
void arb_mat_vector_mul_row(arb_ptr res, arb_srcptr v, const arb_mat_t A, slong prec)
void arb_mat_vector_mul_col(arb_ptr res, const arb_mat_t A, arb_srcptr v, slong prec)

Set `res` to the product $vA$ of the row vector $v$ (of length the
number of rows of $A$) by $A$, respectively the product $Av$ of $A$ by
the column vector $v$, each entry computed by `arb_dot`. The underscore
versions do not allow `res` to alias `v`.
:::

## LU decomposition and solving {#solving}

::: function
int arb_mat_lu_classical(slong * perm, arb_mat_t LU, const arb_mat_t A, slong prec)
int arb_mat_lu_recursive(slong * perm, arb_mat_t LU, const arb_mat_t A, slong prec)
int arb_mat_lu(slong * perm, arb_mat_t LU, const arb_mat_t A, slong prec)

Compute an LU decomposition $PA = LU$ of an $n \times n$ matrix $A$ by
Gaussian elimination with partial pivoting in ball arithmetic. Return
nonzero if $n$ pivots that do not contain zero are found, which proves
that $A$ is invertible; then $L$ (unit lower triangular, without its
diagonal) and $U$ are stored in `LU`, and `perm[i]` is the row of $A$
that became row $i$. Return zero, leaving `perm` and `LU` undefined,
otherwise. `LU` may be the same object as `A`.

::: details
The pivot in each column is the entry not containing zero with the
largest midpoint. The *classical* version eliminates row by row; the
*recursive* version splits the columns in two halves, recursively
decomposes the left half and updates the right half with a matrix
product, so that it benefits from fast multiplication. The default uses
the classical version if either dimension is below 8. In ball
arithmetic the radii of the entries typically grow exponentially with
$n$ even for well-conditioned matrices, so this decomposition is mainly
useful at high precision or for small matrices; see [](arb_mat/solve).
:::
:::

::: function
void arb_mat_solve_tril_classical(arb_mat_t X, const arb_mat_t L, const arb_mat_t B, int unit, slong prec)
void arb_mat_solve_tril_recursive(arb_mat_t X, const arb_mat_t L, const arb_mat_t B, int unit, slong prec)
void arb_mat_solve_tril(arb_mat_t X, const arb_mat_t L, const arb_mat_t B, int unit, slong prec)
void arb_mat_solve_triu_classical(arb_mat_t X, const arb_mat_t U, const arb_mat_t B, int unit, slong prec)
void arb_mat_solve_triu_recursive(arb_mat_t X, const arb_mat_t U, const arb_mat_t B, int unit, slong prec)
void arb_mat_solve_triu(arb_mat_t X, const arb_mat_t U, const arb_mat_t B, int unit, slong prec)

Solve $LX = B$ for a lower triangular $L$, respectively $UX = B$ for an
upper triangular $U$, by substitution. Only the relevant triangle is
read; the diagonal must not contain zero (not checked; otherwise the
output is not finite). If `unit` is nonzero, the diagonal is taken to be
all ones and is not read. The *recursive* versions use block
substitution with matrix products; the defaults use them when $B$ has
at least 40 rows and 40 columns.
:::

::: function
void arb_mat_solve_lu_precomp(arb_mat_t X, const slong * perm, const arb_mat_t LU, const arb_mat_t B, slong prec)

Solves $AX = B$ given an LU decomposition $PA = LU$ computed by
`arb_mat_lu` (assumed valid, not checked). `X` may be the same object
as `B` but not as `LU`.
:::

::: function
int arb_mat_solve(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec)
int arb_mat_solve_lu(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec)
int arb_mat_solve_precond(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec)

Solve $AX = B$ for a square matrix $A$ and $X$, $B$ with $m$ columns.
Return nonzero if $A$ is proved invertible, in which case `X` contains
the exact solution for every pair of matrices in the balls $A$ and $B$;
return zero (with `X` undefined) if $m > 0$ and this could not be
proved. The *lu* version uses Gaussian elimination in ball arithmetic,
the *precond* version preconditions with an approximate inverse
[@HS1967], and the default uses the *lu* version if $n \le 4$ or
$\mathit{prec} > 10n$ and the *precond* version otherwise. See
[](arb_mat/solve).
:::

::: function
int arb_mat_solve_preapprox(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, const arb_mat_t R, const arb_mat_t T, slong prec)

Solves $AX = B$ given an approximate inverse $R$ of $A$ and an
approximate solution $T$, by bounding the error of $T$ (Theorem 10.2 of
[@Rum2010]). Returns nonzero if $\|I - RA\|_\infty < 1$ could be proved,
which shows that $A$ is invertible; `X` then contains the exact
solution. Returns zero (with `X` undefined) otherwise, for example if
$R$ is not a good enough approximation.
:::

::: function
int arb_mat_inv(arb_mat_t X, const arb_mat_t A, slong prec)

Sets `X` to $A^{-1}$ by solving $AX = I$ with `arb_mat_solve`. Returns
nonzero if $A$ is proved invertible (and `X` contains the exact inverse)
and zero otherwise, with `X` undefined.
:::

::: function
void arb_mat_det_lu(arb_t det, const arb_mat_t A, slong prec)
void arb_mat_det_precond(arb_t det, const arb_mat_t A, slong prec)
void arb_mat_det(arb_t det, const arb_mat_t A, slong prec)

Set `det` to an enclosure of the determinant of the square matrix $A$.
The *lu* version uses Gaussian elimination and bounds the remaining
submatrix by Hadamard's inequality when no more pivots can be found;
the *precond* version multiplies $A$ by approximate inverses of the
factors of an approximate LU decomposition and encloses the determinant
of the result, which is close to the identity [@Rum2010]. The default
uses direct formulas for $n \le 3$ and triangular matrices, the *lu*
version if $n \le 10$ or $\mathit{prec} > 10n$, and the *precond*
version otherwise; it throws an exception for a non-square matrix.
:::

## Approximate solving {#approx}

These functions perform floating-point linear algebra on the midpoints
of the input, *without any error control*. The input radii are ignored,
and the output matrices hold the approximate results with zero radii.
They use the same algorithms as the rigorous versions (partial
pivoting, block recursion for large matrices).

::: function
void arb_mat_approx_solve_triu(arb_mat_t X, const arb_mat_t U, const arb_mat_t B, int unit, slong prec)
void arb_mat_approx_solve_tril(arb_mat_t X, const arb_mat_t L, const arb_mat_t B, int unit, slong prec)

Approximately solve $UX = B$ or $LX = B$ for a triangular matrix with
nonzero diagonal (not checked); with `unit` set, the diagonal is taken
to be ones and not read.
:::

::: function
int arb_mat_approx_lu(slong * P, arb_mat_t LU, const arb_mat_t A, slong prec)
void arb_mat_approx_solve_lu_precomp(arb_mat_t X, const slong * perm, const arb_mat_t A, const arb_mat_t B, slong prec)

Compute an approximate LU decomposition with partial pivoting, returning
zero if an exactly zero pivot column is encountered, and solve with such
a decomposition (assumed to come from `arb_mat_approx_lu`).
:::

::: function
int arb_mat_approx_solve(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec)
int arb_mat_approx_inv(arb_mat_t X, const arb_mat_t A, slong prec)

Approximately solve $AX = B$ or invert $A$. Return zero if the
approximate LU decomposition fails, nonzero otherwise; a nonzero return
value says nothing about the accuracy of the result.
:::

## Cholesky and LDL decompositions {#cholesky}

The functions in this section read only the lower triangular part of
their symmetric input matrix $A$ (including the diagonal).

::: function
int _arb_mat_cholesky_banachiewicz(arb_mat_t A, slong prec)
int arb_mat_cho(arb_mat_t L, const arb_mat_t A, slong prec)

Compute the Cholesky decomposition $A = LL^T$ with $L$ lower triangular,
by the Cholesky–Banachiewicz algorithm. Return nonzero if all pivots
are certainly positive, which proves that every symmetric matrix
represented by the lower triangle of $A$ is positive definite; return
zero if $A$ is not positive definite or the precision is insufficient.
The underscore version works in place and leaves the strict upper
triangle undefined.
:::

::: function
void arb_mat_solve_cho_precomp(arb_mat_t X, const arb_mat_t L, const arb_mat_t B, slong prec)

Solves $AX = B$ given the Cholesky factor $L$ of $A$ computed by
`arb_mat_cho` (assumed valid, not checked). `X` may be the same object
as `B` but not as `L`.
:::

::: function
int arb_mat_spd_solve(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec)

Solves $AX = B$ for a symmetric positive definite matrix $A$ using
`arb_mat_cho`. Returns nonzero if the decomposition succeeds, which
proves that $A$ is positive definite and that `X` contains the exact
solution, and zero (if $B$ has at least one column) otherwise, with `X`
undefined.
:::

::: function
void arb_mat_inv_cho_precomp(arb_mat_t X, const arb_mat_t L, slong prec)
int arb_mat_spd_inv(arb_mat_t X, const arb_mat_t A, slong prec)

Set `X` to $A^{-1}$ for a symmetric positive definite $A$, given its
Cholesky factor $L$ (assumed valid) or computing it. The inverse is
obtained by the method of [@Kri2013], which computes fewer intermediate
results than solving $AX = I$. `arb_mat_spd_inv` returns zero (with `X`
undefined) if the Cholesky decomposition fails.
:::

::: function
int _arb_mat_ldl_inplace(arb_mat_t A, slong prec)
int _arb_mat_ldl_golub_and_van_loan(arb_mat_t A, slong prec)
int arb_mat_ldl(arb_mat_t res, const arb_mat_t A, slong prec)

Compute the decomposition $A = LDL^T$ with $L$ unit lower triangular and
$D$ diagonal, stored in `res`: the strict lower triangle holds $L$ and
the diagonal holds $D$. Return nonzero if all entries of $D$ are
certainly positive, which proves that $A$ is positive definite, and
zero otherwise. The underscore versions work in place and leave the
strict upper triangle undefined; the default uses Algorithm 4.1.2 of
[@GVL1996] (`_arb_mat_ldl_golub_and_van_loan`).
:::

::: function
void arb_mat_solve_ldl_precomp(arb_mat_t X, const arb_mat_t L, const arb_mat_t B, slong prec)
void arb_mat_inv_ldl_precomp(arb_mat_t X, const arb_mat_t L, slong prec)

Solve $AX = B$, respectively compute $A^{-1}$ (by the method of
[@Kri2013]), given the $LDL^T$ decomposition of $A$ encoded in `L` as
computed by `arb_mat_ldl` (assumed valid). `X` may be the same object
as `B` but not as `L`.
:::

## Characteristic polynomial and companion matrix {#charpoly}

::: function
void _arb_mat_charpoly(arb_ptr poly, const arb_mat_t mat, slong prec)
void arb_mat_charpoly(arb_poly_t poly, const arb_mat_t mat, slong prec)

Set `poly` to the characteristic polynomial $\det(xI - A)$ of the
square $n \times n$ matrix $A$ (an exception is raised otherwise). The
underscore version writes $n + 1$ coefficients. Wrappers of
`_gr_mat_charpoly_berkowitz`, the division-free algorithm of Berkowitz
with $O(n^4)$ operations, which avoids divisions by balls that contain
zero ([](algorithms/linear-algebra/charpoly)). A non-finite matrix
gives indeterminate coefficients.
:::

::: function
void _arb_mat_companion(arb_mat_t mat, arb_srcptr poly, slong prec)
void arb_mat_companion(arb_mat_t mat, const arb_poly_t poly, slong prec)

Set the $n \times n$ matrix `mat` to the companion matrix of the
polynomial of degree $n$ (the underscore version reads $n + 1$
coefficients). If the leading coefficient cannot be inverted, the
result is indeterminate. Wrappers of `_gr_mat_companion` and
`gr_mat_companion`.
:::

## Matrix functions {#functions}

::: function
void arb_mat_exp_taylor_sum(arb_mat_t S, const arb_mat_t A, slong N, slong prec)

Sets `S` to the truncated exponential series $\sum_{k=0}^{N-1} A^k /
k!$, using rectangular splitting with $O(\sqrt{N})$ matrix
multiplications. The factorials are generated by their recurrence so
that the scalars are small integers, and the divisions are postponed by
computing partial factorials of length $O(\sqrt{N})$ [@Joh2014b].
:::

::: function
void arb_mat_exp(arb_mat_t B, const arb_mat_t A, slong prec)

Sets `B` to the exponential $\exp(A) = \sum_{k \ge 0} A^k/k!$ of the
square matrix $A$.

::: details
The function uses scaling and squaring: it evaluates
$\exp(2^{-r} A)$ by `arb_mat_exp_taylor_sum` with $N$ terms and squares
the result $r$ times, at a working precision of $\mathit{prec} + 3
\log_2 \mathit{prec}$ bits, with $r$ chosen so that $\|2^{-r}
A\|_\infty \approx 2^{-\mathit{prec}^{1/4}}$. The truncation error in
each entry is bounded by
$$\left\|\sum_{k=N}^{\infty} \frac{(2^{-r} A)^k}{k!}\right\|_{\infty}
\le \sum_{k=N}^{\infty} \frac{(2^{-r} \|A\|_{\infty})^k}{k!},$$
computed by `mag_exp_tail`. The sparsity structure of $A$ is used: a
diagonal matrix is exponentiated entrywise, a nilpotent structure
limits the number of terms, and no truncation error is added to entries
that are zero for structural reasons. The structure is found with the
graph functions of [](bool_mat/graphs) applied to the pattern of nonzero
entries. If $\|A\|_\infty > 2^{2w}$ for
the working precision $w$, the result is indeterminate.
:::
:::

::: function
void arb_mat_trace(arb_t trace, const arb_mat_t mat, slong prec)

Sets `trace` to the sum of the diagonal entries of the square matrix
`mat`.
:::

::: function
void _arb_mat_diag_prod(arb_t res, const arb_mat_t mat, slong a, slong b, slong prec)
void arb_mat_diag_prod(arb_t res, const arb_mat_t mat, slong prec)

Set `res` to the product of the diagonal entries of `mat`; the
underscore version takes the product of the entries with index $a \le
i < b$ (which must be in range; the empty product is 1).
:::

## Sparsity structure {#sparsity}

::: function
slong arb_mat_count_is_zero(const arb_mat_t mat)
slong arb_mat_count_not_is_zero(const arb_mat_t mat)

Return the number of entries of `mat` that are exactly zero
(`arb_is_zero`), respectively not exactly zero.
:::

## Component and error operations {#components}

::: function
void arb_mat_get_mid(arb_mat_t B, const arb_mat_t A)

Sets `B` to the matrix of midpoints of `A`, with zero radii.
:::

::: function
void arb_mat_add_error_mag(arb_mat_t mat, const mag_t err)

Adds `err` to the radius of every entry of `mat`.
:::

## Eigenvalues {#eigenvalues}

There are no eigenvalue functions for real matrices in this module.
Convert the matrix to an `acb_mat_t` with `acb_mat_set_arb_mat` and use
the functions in [](acb_mat/eigenvalues).

## LLL reduction of quadratic forms {#lll}

::: function
int arb_mat_spd_get_fmpz_mat(fmpz_mat_t B, const arb_mat_t A, slong prec)

Attempts to set `B` to the symmetric positive definite integer matrix
obtained by rounding $2^{\mathit{prec}} A$ entrywise, using the
midpoint of the intersection of $a_{jk}$ and $a_{kj}$ for both entries.
Returns 1 on success, and 0 (with `B` undefined) if some pair $a_{jk}$,
$a_{kj}$ does not overlap or the rounded matrix is not positive
definite. As for `arf_get_fmpz`, huge exponents can exhaust memory.
:::

::: function
void arb_mat_spd_lll_reduce(fmpz_mat_t U, const arb_mat_t A, double delta, double eta, slong prec)

Given a symmetric positive definite matrix $A$, sets `U` to an
invertible integer matrix such that $U A U^T$ is close to being
LLL-reduced with the parameters $\delta$, $\eta$: if
`arb_mat_spd_get_fmpz_mat` succeeds at precision `prec`, `fmpz_lll` is
applied to the rounded Gram matrix, and otherwise `U` is the identity.
If $A$ is not finite, `U` is not modified.
:::

::: function
int arb_mat_spd_is_lll_reduced(const arb_mat_t A, double delta, double eta, slong prec)

Given a symmetric positive definite matrix $A$ (only the lower triangle
is read; neither property is checked), returns nonzero if $A$ is
certainly LLL-reduced for the parameters $\delta$, $\eta$ as a Gram
matrix: the Gram–Schmidt coefficients satisfy $|\mu_{j,k}| \le \eta$
and $\delta \|b_{k-1}^*\|^2 \le \|b_k^*\|^2 + \mu_{k,k-1}^2
\|b_{k-1}^*\|^2$.
:::

## Internal functions {#internal}

::: function
void arb_mat_swap_rows(arb_mat_t mat, slong * perm, slong r, slong s)

Swaps rows $r$ and $s$ of `mat` and, if `perm` is not `NULL`, the
entries $r$ and $s$ of `perm`.
:::

::: function
slong arb_mat_find_pivot_partial(const arb_mat_t mat, slong start_row, slong end_row, slong c)

Returns the index $i$ with $\mathit{start\_row} \le i <
\mathit{end\_row}$ of the entry in column $c$ that does not contain
zero and has the largest midpoint in absolute value, or $-1$ if every
entry in that range contains zero.
:::

::: function
void _arb_mat_addmul_rad_mag_fast(arb_mat_t C, mag_srcptr A, mag_srcptr B, slong ar, slong ac, slong bc)

Adds to the radii of `C` an upper bound for the product of the
nonnegative matrices $A$ ($\mathit{ar} \times \mathit{ac}$, row-major
array) and $B$ ($\mathit{ac} \times \mathit{bc}$, given as a column-major
array). The exponents of all entries must be small (this is not
checked), which `arb_mat_mul_block` ensures. Products are computed in
double precision when the exponent range permits and the result is
multiplied by $1 + 2.01(\mathit{ac} + 1) 2^{-52}$ to account for
rounding errors.
:::
