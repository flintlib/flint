---
title: acb_mat.h
header: acb_mat.h
description: Dense matrices with complex ball entries: rigorous arithmetic, verified solving and inversion, determinants, the matrix exponential, and certified eigenvalues and eigenvectors.
category: Real and complex numbers
order: 10
---

This module implements dense matrices whose entries are complex balls
(`acb_t`), with the same interface as [](arb_mat) for arithmetic,
solving linear systems, inversion, determinants, the characteristic
polynomial and the exponential. A matrix of balls represents the set of
complex matrices whose entries lie in the balls, and the functions
return enclosures of the exact results for every matrix in that set.
In addition, the module computes eigenvalues and eigenvectors: the QR
algorithm gives floating-point approximations, and the methods of Rump,
Miyajima and van der Hoeven and Mourrain turn them into certified
enclosures of the individual eigenvalues, of clusters of eigenvalues,
or of the whole spectrum ([](acb_mat/eigenvalues)). Real matrices use
these functions after conversion with `acb_mat_set_arb_mat`.

Large products of complex matrices are computed by `acb_mat_mul_reorder`
as products of real matrices with `arb_mat_mul`
([](algorithms/ball-arithmetic/matrix-multiplication)). The solvers
follow [](algorithms/ball-arithmetic/linear-systems), and the
eigenvalue algorithms [](algorithms/ball-arithmetic/eigenvalues). The
characteristic polynomial, the companion matrix and powers are computed
by the generic implementations in [](gr_mat) with the complex ball
context.

::: summary algorithms/ball-arithmetic/eigenvalues :::

## Types and conventions {#conventions}

::: type
acb_mat_struct
acb_mat_t

A matrix with complex ball entries. The structure has the fields
`entries` (a pointer to the entry in the top left corner), `r` and `c`
(the numbers of rows and columns) and `stride` (the distance in
`acb_struct` between the starts of consecutive rows). An `acb_mat_t` is
an array of one `acb_mat_struct`, so it is passed by reference.
:::

::: macro
acb_mat_entry(mat, i, j)

Evaluates to a pointer to the entry in row $i$ and column $j$ of `mat`
(no bounds check).
:::

::: function
acb_ptr acb_mat_entry_ptr(acb_mat_t mat, slong i, slong j)

Returns a pointer to the entry in row $i$ and column $j$ of `mat`; the
same as `acb_mat_entry`, as a function.
:::

::: macro
acb_mat_nrows(mat)
acb_mat_ncols(mat)

Evaluate to the number of rows and columns of `mat`.
:::

The conventions are those of [](arb_mat#conventions): dimensions are
fixed at initialisation and must be compatible; outputs may alias
inputs; functions that can fail numerically return nonzero only when
the result is certified; predicates return 1 only if the property
certainly holds; the `approx` functions compute floating-point results
without error bounds. `acb_mat.h` does not include `acb.h`, which
programs that work with the entries must include.

## Memory management {#memory}

::: function
void acb_mat_init(acb_mat_t mat, slong r, slong c)

Initialises `mat` as the zero matrix with `r` rows and `c` columns.
:::

::: function
void acb_mat_clear(acb_mat_t mat)

Frees the memory used by `mat`.
:::

::: function
void acb_mat_swap(acb_mat_t mat1, acb_mat_t mat2)
void acb_mat_swap_entrywise(acb_mat_t mat1, acb_mat_t mat2)

Swap two matrices. `acb_mat_swap` swaps the structures (also when the
dimensions differ); `acb_mat_swap_entrywise` swaps the entries one by
one and requires equal dimensions, which also works for windows.
:::

::: function
slong acb_mat_allocated_bytes(const acb_mat_t x)

Returns the number of bytes allocated on the heap by `x`, excluding
`sizeof(acb_mat_struct)`.
:::

::: function
void acb_mat_window_init(acb_mat_t window, const acb_mat_t mat, slong r1, slong c1, slong r2, slong c2)
void acb_mat_window_clear(acb_mat_t window)

Initialise `window` as a view of the submatrix of `mat` formed by rows
$r_1 \le i < r_2$ and columns $c_1 \le j < c_2$, sharing its entries,
and release it again (which frees nothing). Require $0 \le r_1 \le r_2
\le$ the number of rows and likewise for the columns. A window must not
be used after `mat` has been cleared.
:::

## Conversions {#conversions}

::: function
void acb_mat_set(acb_mat_t dest, const acb_mat_t src)
void acb_mat_set_fmpz_mat(acb_mat_t dest, const fmpz_mat_t src)
void acb_mat_set_round_fmpz_mat(acb_mat_t dest, const fmpz_mat_t src, slong prec)
void acb_mat_set_fmpq_mat(acb_mat_t dest, const fmpq_mat_t src, slong prec)
void acb_mat_set_arb_mat(acb_mat_t dest, const arb_mat_t src)
void acb_mat_set_round_arb_mat(acb_mat_t dest, const arb_mat_t src, slong prec)

Set `dest` to `src`, which must have the same dimensions: exactly, or
with the entries rounded to `prec` bits. Real input gives entries with
imaginary part exactly zero.
:::

::: function
void acb_mat_get_real(arb_mat_t re, const acb_mat_t mat)
void acb_mat_get_imag(arb_mat_t im, const acb_mat_t mat)
void acb_mat_set_real_imag(acb_mat_t mat, const arb_mat_t re, const arb_mat_t im)

Extract the real or imaginary parts of the entries of `mat`, or set
`mat` to $\mathit{re} + i \cdot \mathit{im}$, exactly.
:::

## Random generation {#random}

::: function
void acb_mat_randtest(acb_mat_t mat, flint_rand_t state, slong prec, slong mag_bits)

Sets `mat` to a random matrix of random density whose nonzero entries
have midpoints of up to `prec` bits and exponents of up to `mag_bits`
bits.
:::

::: function
void acb_mat_randtest_eig(acb_mat_t mat, flint_rand_t state, acb_srcptr E, slong prec)

Sets the square matrix `mat` to a random matrix with the eigenvalues
given by the vector $E$ (with their multiplicities), constructed as $Q
U Q^{*}$ with $Q$ the exponential of a random skew-Hermitian matrix
(unitary up to rounding) and $U$ a random upper triangular matrix with
diagonal $E$. The result is computed in ball arithmetic and contains a
matrix with exactly these eigenvalues. For testing.
:::

## Input and output {#io}

::: function
void acb_mat_printd(const acb_mat_t mat, slong digits)
void acb_mat_fprintd(FILE * file, const acb_mat_t mat, slong digits)

Print `mat` to standard output or to `file`, one row per line as `[a,
b, ...]`, each entry printed with `acb_fprintd` to `digits` decimal
digits.
:::

## Comparisons {#comparisons}

::: function
int acb_mat_equal(const acb_mat_t mat1, const acb_mat_t mat2)

Returns whether the matrices have the same dimensions and identical
balls as entries.
:::

::: function
int acb_mat_overlaps(const acb_mat_t mat1, const acb_mat_t mat2)

Returns whether the matrices have the same dimensions and each entry of
`mat1` overlaps the corresponding entry of `mat2`.
:::

::: function
int acb_mat_contains(const acb_mat_t mat1, const acb_mat_t mat2)
int acb_mat_contains_fmpz_mat(const acb_mat_t mat1, const fmpz_mat_t mat2)
int acb_mat_contains_fmpq_mat(const acb_mat_t mat1, const fmpq_mat_t mat2)

Return whether the matrices have the same dimensions and each entry of
`mat2` is contained in the corresponding entry of `mat1`.
:::

::: function
int acb_mat_eq(const acb_mat_t mat1, const acb_mat_t mat2)
int acb_mat_ne(const acb_mat_t mat1, const acb_mat_t mat2)

Return whether `mat1` and `mat2` certainly represent the same matrix
(same dimensions and all entries exact and equal), respectively
certainly represent different matrices (different dimensions, or some
pair of entries that do not overlap).
:::

::: function
int acb_mat_is_empty(const acb_mat_t mat)
int acb_mat_is_square(const acb_mat_t mat)

Return whether `mat` has zero rows or zero columns, respectively as
many rows as columns.
:::

::: function
int acb_mat_is_exact(const acb_mat_t mat)
int acb_mat_is_zero(const acb_mat_t mat)
int acb_mat_is_finite(const acb_mat_t mat)

Return whether all entries of `mat` have zero radius, are exactly zero,
or are finite.
:::

::: function
int acb_mat_is_real(const acb_mat_t mat)

Returns whether all entries of `mat` have imaginary part exactly zero.
:::

::: function
int acb_mat_is_triu(const acb_mat_t mat)
int acb_mat_is_tril(const acb_mat_t mat)
int acb_mat_is_diag(const acb_mat_t mat)

Return whether all entries below, above, respectively off the main
diagonal of `mat` are exactly zero. The matrix need not be square.
:::

## Special matrices {#special-matrices}

::: function
void acb_mat_zero(acb_mat_t mat)
void acb_mat_one(acb_mat_t mat)
void acb_mat_ones(acb_mat_t mat)
void acb_mat_onei(acb_mat_t mat)
void acb_mat_indeterminate(acb_mat_t mat)

Set all entries of `mat` to zero; the diagonal entries to one and the
others to zero; all entries to one; the diagonal entries to $i$ and the
others to zero; or all entries to indeterminate balls.
:::

::: function
void acb_mat_dft(acb_mat_t mat, int type, slong prec)

Sets `mat` to the unitary DFT matrix of order $n$, the smaller dimension
of `mat`,
$$A_{j,k} = \frac{1}{\sqrt{n}} \exp\left(-\frac{2\pi i}{n} jk\right),$$
which satisfies $A^{-1} = A^{*}$ (if `mat` is not square, the formula
is used for all entries). `type` is ignored and should be 0.
:::

## Transpose and conjugate {#transpose}

::: function
void acb_mat_transpose(acb_mat_t dest, const acb_mat_t src)
void acb_mat_conjugate_transpose(acb_mat_t dest, const acb_mat_t src)

Set `dest` to the transpose, respectively the conjugate transpose, of
`src`, exactly. Aliasing is allowed for square matrices.
:::

::: function
void acb_mat_conjugate(acb_mat_t dest, const acb_mat_t src)

Sets `dest` to the entrywise complex conjugate of `src`.
:::

## Norms {#norms}

::: function
void acb_mat_bound_inf_norm(mag_t b, const acb_mat_t A)

Sets `b` to an upper bound for the infinity norm $\max_i \sum_j
|a_{ij}|$ of every matrix in `A`.
:::

::: function
void acb_mat_frobenius_norm(arb_t res, const acb_mat_t A, slong prec)
void acb_mat_bound_frobenius_norm(mag_t res, const acb_mat_t A)

Set `res` to an enclosure, respectively an upper bound, of the
Frobenius norm $(\sum_{i,j} |a_{ij}|^2)^{1/2}$.
:::

## Arithmetic {#arithmetic}

::: function
void acb_mat_neg(acb_mat_t dest, const acb_mat_t src)

Sets `dest` to $-\mathit{src}$, exactly.
:::

::: function
void acb_mat_add(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, slong prec)
void acb_mat_sub(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, slong prec)

Set `res` to the sum or difference of matrices of equal dimensions.
:::

::: function
void acb_mat_mul_classical(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, slong prec)
void acb_mat_mul_threaded(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, slong prec)
void acb_mat_mul_reorder(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, slong prec)
void acb_mat_mul(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, slong prec)

Set `res` to the matrix product. The *classical* version computes each
entry with `acb_dot`, and the *threaded* version does the same with the
rows of the output distributed over `flint_get_num_threads()` threads.
The *reorder* version writes the operands as $A = A_r + iA_i$, $B = B_r
+ iB_i$ with real matrices and computes the product with real matrix
multiplications by `arb_mat_mul` (one if both operands are real, two if
one is, otherwise four, $A_r B_r - A_i B_i$ and $A_r B_i + A_i B_r$),
which benefits from the block algorithm of `arb_mat_mul_block`.

::: details
`acb_mat_mul` uses the *reorder* version if the inner dimension and the
number of columns of $B$ are at least 20, the midpoints have fewer than
8000 bits and these dimensions are at least $5 + b/64$, where $b$ is
the smaller of `prec` and the largest midpoint size of the entries;
otherwise the threaded version if several threads are available and
the work is large enough, and the classical version otherwise.
:::
:::

::: function
void acb_mat_mul_entrywise(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, slong prec)

Sets `res` to the entrywise product of matrices of equal dimensions.
:::

::: function
void acb_mat_sqr_classical(acb_mat_t res, const acb_mat_t mat, slong prec)
void acb_mat_sqr(acb_mat_t res, const acb_mat_t mat, slong prec)

Set `res` to $A^2$ for a square matrix $A$.
:::

::: function
void acb_mat_pow_ui(acb_mat_t res, const acb_mat_t mat, ulong exp, slong prec)

Sets `res` to $A^e$ for a square matrix $A$, by binary exponentiation.
A wrapper of `gr_mat_pow_ui`.
:::

::: function
void acb_mat_approx_mul(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, slong prec)

Approximate matrix multiplication: the input radii are ignored and the
output is a floating-point approximation of the product of the
midpoints. The output radii carry no meaning.
:::

## Scalar arithmetic {#scalar}

::: function
void acb_mat_scalar_mul_2exp_si(acb_mat_t B, const acb_mat_t A, slong c)

Sets `B` to $2^c A$, exactly.
:::

::: function
void acb_mat_scalar_addmul_si(acb_mat_t B, const acb_mat_t A, slong c, slong prec)
void acb_mat_scalar_addmul_fmpz(acb_mat_t B, const acb_mat_t A, const fmpz_t c, slong prec)
void acb_mat_scalar_addmul_arb(acb_mat_t B, const acb_mat_t A, const arb_t c, slong prec)
void acb_mat_scalar_addmul_acb(acb_mat_t B, const acb_mat_t A, const acb_t c, slong prec)

Set `B` to $B + cA$.
:::

::: function
void acb_mat_scalar_mul_si(acb_mat_t B, const acb_mat_t A, slong c, slong prec)
void acb_mat_scalar_mul_fmpz(acb_mat_t B, const acb_mat_t A, const fmpz_t c, slong prec)
void acb_mat_scalar_mul_arb(acb_mat_t B, const acb_mat_t A, const arb_t c, slong prec)
void acb_mat_scalar_mul_acb(acb_mat_t B, const acb_mat_t A, const acb_t c, slong prec)

Set `B` to $cA$.
:::

::: function
void acb_mat_scalar_div_si(acb_mat_t B, const acb_mat_t A, slong c, slong prec)
void acb_mat_scalar_div_fmpz(acb_mat_t B, const acb_mat_t A, const fmpz_t c, slong prec)
void acb_mat_scalar_div_arb(acb_mat_t B, const acb_mat_t A, const arb_t c, slong prec)
void acb_mat_scalar_div_acb(acb_mat_t B, const acb_mat_t A, const acb_t c, slong prec)

Set `B` to $A / c$.
:::

## Matrix-vector products {#vector}

::: function
void _acb_mat_vector_mul_row(acb_ptr res, acb_srcptr v, const acb_mat_t A, slong prec)
void _acb_mat_vector_mul_col(acb_ptr res, const acb_mat_t A, acb_srcptr v, slong prec)
void acb_mat_vector_mul_row(acb_ptr res, acb_srcptr v, const acb_mat_t A, slong prec)
void acb_mat_vector_mul_col(acb_ptr res, const acb_mat_t A, acb_srcptr v, slong prec)

Set `res` to the product $vA$ of the row vector $v$ (of length the
number of rows of $A$) by $A$, respectively the product $Av$ of $A$ by
the column vector $v$, each entry computed by `acb_dot`. The underscore
versions do not allow `res` to alias `v`.
:::

## LU decomposition and solving {#solving}

::: function
int acb_mat_lu_classical(slong * perm, acb_mat_t LU, const acb_mat_t A, slong prec)
int acb_mat_lu_recursive(slong * perm, acb_mat_t LU, const acb_mat_t A, slong prec)
int acb_mat_lu(slong * perm, acb_mat_t LU, const acb_mat_t A, slong prec)

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
void acb_mat_solve_tril_classical(acb_mat_t X, const acb_mat_t L, const acb_mat_t B, int unit, slong prec)
void acb_mat_solve_tril_recursive(acb_mat_t X, const acb_mat_t L, const acb_mat_t B, int unit, slong prec)
void acb_mat_solve_tril(acb_mat_t X, const acb_mat_t L, const acb_mat_t B, int unit, slong prec)
void acb_mat_solve_triu_classical(acb_mat_t X, const acb_mat_t U, const acb_mat_t B, int unit, slong prec)
void acb_mat_solve_triu_recursive(acb_mat_t X, const acb_mat_t U, const acb_mat_t B, int unit, slong prec)
void acb_mat_solve_triu(acb_mat_t X, const acb_mat_t U, const acb_mat_t B, int unit, slong prec)

Solve $LX = B$ for a lower triangular $L$, respectively $UX = B$ for an
upper triangular $U$, by substitution. Only the relevant triangle is
read; the diagonal must not contain zero (not checked; otherwise the
output is not finite). If `unit` is nonzero, the diagonal is taken to be
all ones and is not read. The *recursive* versions use block
substitution with matrix products; the defaults use them when $B$ has
at least 40 rows and 40 columns.
:::

::: function
void acb_mat_solve_lu_precomp(acb_mat_t X, const slong * perm, const acb_mat_t LU, const acb_mat_t B, slong prec)

Solves $AX = B$ given an LU decomposition $PA = LU$ computed by
`acb_mat_lu` (assumed valid, not checked). `X` may be the same object
as `B` but not as `LU`.
:::

::: function
int acb_mat_solve(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, slong prec)
int acb_mat_solve_lu(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, slong prec)
int acb_mat_solve_precond(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, slong prec)

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
int acb_mat_inv(acb_mat_t X, const acb_mat_t A, slong prec)

Sets `X` to $A^{-1}$ by solving $AX = I$ with `acb_mat_solve`. Returns
nonzero if $A$ is proved invertible (and `X` contains the exact inverse)
and zero otherwise, with `X` undefined.
:::

::: function
void acb_mat_det_lu(acb_t det, const acb_mat_t A, slong prec)
void acb_mat_det_precond(acb_t det, const acb_mat_t A, slong prec)
void acb_mat_det(acb_t det, const acb_mat_t A, slong prec)

Set `det` to an enclosure of the determinant of the square matrix $A$.
The *lu* version uses Gaussian elimination and bounds the remaining
submatrix by Hadamard's inequality when no more pivots can be found;
the *precond* version multiplies $A$ by approximate inverses of the
factors of an approximate LU decomposition and encloses the determinant
of the result, which is close to the identity [@Rum2010]. The default
uses direct formulas for $n \le 3$ and triangular matrices, the *lu*
version if $n \le 14$ or $\mathit{prec} > 10n$, and the *precond*
version otherwise; it throws an exception for a non-square matrix.
:::

## Approximate solving {#approx}

These functions perform floating-point linear algebra on the midpoints
of the input, *without any error control*. The input radii are ignored,
and the output matrices hold the approximate results with zero radii.
They use the same algorithms as the rigorous versions (partial
pivoting, block recursion for large matrices).

::: function
void acb_mat_approx_solve_triu(acb_mat_t X, const acb_mat_t U, const acb_mat_t B, int unit, slong prec)
void acb_mat_approx_solve_tril(acb_mat_t X, const acb_mat_t L, const acb_mat_t B, int unit, slong prec)

Approximately solve $UX = B$ or $LX = B$ for a triangular matrix with
nonzero diagonal (not checked); with `unit` set, the diagonal is taken
to be ones and not read.
:::

::: function
int acb_mat_approx_lu(slong * P, acb_mat_t LU, const acb_mat_t A, slong prec)
void acb_mat_approx_solve_lu_precomp(acb_mat_t X, const slong * perm, const acb_mat_t A, const acb_mat_t B, slong prec)

Compute an approximate LU decomposition with partial pivoting, returning
zero if an exactly zero pivot column is encountered, and solve with such
a decomposition (assumed to come from `acb_mat_approx_lu`).
:::

::: function
int acb_mat_approx_solve(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, slong prec)
int acb_mat_approx_inv(acb_mat_t X, const acb_mat_t A, slong prec)

Approximately solve $AX = B$ or invert $A$. Return zero if the
approximate LU decomposition fails, nonzero otherwise; a nonzero return
value says nothing about the accuracy of the result.
:::

## Characteristic polynomial and companion matrix {#charpoly}

::: function
void _acb_mat_charpoly(acb_ptr poly, const acb_mat_t mat, slong prec)
void acb_mat_charpoly(acb_poly_t poly, const acb_mat_t mat, slong prec)

Set `poly` to the characteristic polynomial $\det(xI - A)$ of the
square $n \times n$ matrix $A$ (an exception is raised otherwise). The
underscore version writes $n + 1$ coefficients. Wrappers of
`_gr_mat_charpoly_berkowitz`, the division-free algorithm of Berkowitz
with $O(n^4)$ operations, which avoids divisions by balls that contain
zero ([](algorithms/linear-algebra/charpoly)). A non-finite matrix
gives indeterminate coefficients.
:::

::: function
void _acb_mat_companion(acb_mat_t mat, acb_srcptr poly, slong prec)
void acb_mat_companion(acb_mat_t mat, const acb_poly_t poly, slong prec)

Set the $n \times n$ matrix `mat` to the companion matrix of the
polynomial of degree $n$ (the underscore version reads $n + 1$
coefficients). If the leading coefficient cannot be inverted, the
result is indeterminate. Wrappers of `_gr_mat_companion` and
`gr_mat_companion`.
:::

## Matrix functions {#functions}

::: function
void acb_mat_exp_taylor_sum(acb_mat_t S, const acb_mat_t A, slong N, slong prec)

Sets `S` to the truncated exponential series $\sum_{k=0}^{N-1} A^k /
k!$, using rectangular splitting with $O(\sqrt{N})$ matrix
multiplications. The factorials are generated by their recurrence so
that the scalars are small integers, and the divisions are postponed by
computing partial factorials of length $O(\sqrt{N})$ [@Joh2014b].
:::

::: function
void acb_mat_exp(acb_mat_t B, const acb_mat_t A, slong prec)

Sets `B` to the exponential $\exp(A) = \sum_{k \ge 0} A^k/k!$ of the
square matrix $A$.

::: details
The function uses scaling and squaring: it evaluates
$\exp(2^{-r} A)$ by `acb_mat_exp_taylor_sum` with $N$ terms and squares
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
the working precision $w$, the result is indeterminate. A matrix with
real entries is exponentiated with `arb_mat_exp`.
:::
:::

::: function
void acb_mat_trace(acb_t trace, const acb_mat_t mat, slong prec)

Sets `trace` to the sum of the diagonal entries of the square matrix
`mat`.
:::

::: function
void _acb_mat_diag_prod(acb_t res, const acb_mat_t mat, slong a, slong b, slong prec)
void acb_mat_diag_prod(acb_t res, const acb_mat_t mat, slong prec)

Set `res` to the product of the diagonal entries of `mat`; the
underscore version takes the product of the entries with index $a \le
i < b$ (which must be in range; the empty product is 1).
:::

## Component and error operations {#components}

::: function
void acb_mat_get_mid(acb_mat_t B, const acb_mat_t A)

Sets `B` to the matrix of midpoints of `A`, with zero radii.
:::

::: function
void acb_mat_add_error_mag(acb_mat_t mat, const mag_t err)

Adds `err` to the radius of every entry of `mat`.
:::

## Eigenvalues and eigenvectors {#eigenvalues}

The functions in this section are described in more detail, with
examples, in [](acb_mat/eigenvalues). They are experimental: there are
matrices for which the algorithms fail to converge even as `prec` is
increased, or for which the error bounds are much worse than necessary.
A similarity transformation $T^{-1} A T$ can help: balancing, if $A$ is
badly scaled, or a random $T$, if an algorithm fails although $A$ is
well-scaled (for example with multiple eigenvalues).

::: function
int acb_mat_approx_eig_qr(acb_ptr E, acb_mat_t L, acb_mat_t R, const acb_mat_t A, const mag_t tol, slong maxiter, slong prec)

Computes floating-point approximations of the $n$ eigenvalues of the $n
\times n$ matrix $A$ (in no particular order) in `E`, and optionally of
the corresponding left eigenvectors in the rows of `L` and right
eigenvectors in the columns of `R` (either may be `NULL`), without error
bounds. Returns nonzero if the QR iteration converged numerically; this
is only a heuristic, and the output may be accurate even if zero is
returned.

::: details
The matrix is reduced to Hessenberg form by Householder reflections and
the eigenvalues are computed by the implicitly shifted QR algorithm,
working on the midpoints at `prec` bits. `tol` (default
$2^{-\mathit{prec}} / (100n)$ if `NULL`) is the tolerance for
negligible subdiagonal entries, and `maxiter` (default $14n +
\mathit{prec}/10$ if 0 or negative) bounds the number of iterations.
The eigenvectors are computed from the triangular Schur form by
back-substitution and transformed back with the accumulated unitary
transformation.
:::
:::

::: function
void acb_mat_eig_global_enclosure(mag_t eps, const acb_mat_t A, acb_srcptr E, const acb_mat_t R, slong prec)

Given approximations $E$ of the eigenvalues and $R$ of the right
eigenvectors of the $n \times n$ matrix $A$, computes $\varepsilon$
such that every eigenvalue of $A$ lies within distance $\varepsilon$ of
some $\hat\lambda_k \in E$, by Theorem 1 of [@Miy2010], an extension of
the Bauer–Fike theorem. If $\|YR - I\|_\infty \ge 1$ for the
approximate inverse $Y$ of $R$, `eps` is set to infinity.

::: details
The union of the disks $\{z : |z - \hat\lambda_k| \le \varepsilon\}$
encloses the spectrum, but a single disk need not contain exactly one
eigenvalue; when the approximations are good, the disks usually
correspond to the eigenvalues one to one. The bound is
$$\varepsilon = \frac{\|Y(AR - RD)\|_\infty}{1 - \|YR - I\|_\infty},$$
with $D = \operatorname{diag}(E)$ and $Y \approx R^{-1}$. It certifies
regions of the complex plane (half-planes, disks, annuli) that contain
the whole spectrum, cheaply and without assumptions on $A$.
:::
:::

::: function
void acb_mat_eig_enclosure_rump(acb_t lambda, acb_mat_t J, acb_mat_t R, const acb_mat_t A, const acb_t lambda_approx, const acb_mat_t R_approx, slong prec)

Given an approximate eigenvalue `lambda_approx` of $A$ and an $n \times
k$ matrix `R_approx` of approximate eigenvectors (a basis of the
invariant subspace of a presumed cluster of $k$ eigenvalues), computes a
ball `lambda` that contains at least $k$ eigenvalues of $A$ (counted
with multiplicity) and an $n \times k$ matrix `R` that contains a basis
of the corresponding invariant subspace. If `J` is not `NULL`, it is set
to a $k \times k$ matrix $J \approx \lambda I$ such that $AR - RJ$
contains zero; the Jordan form of $A$ then has a $k \times k$ block
equal to the Jordan form of a matrix in $J$. For $k = 1$, $AR - R
\lambda$ contains zero. On failure, the outputs are set to
indeterminate values.

::: details
The algorithm is that of Section 13.4 of [@Rum2010] (the INTLAB
function `verifyeig`): with a preconditioner $R_0 \approx ((A -
\tilde\lambda I) U U^T - \tilde X V^T)^{-1}$, where the rows selected
by $V$ (for $k = 1$ the row of $\tilde X$ of largest absolute value,
otherwise the last $k$ rows) are frozen, a fixed point iteration with
$\varepsilon$-inflation is run for at most $5 + b$ steps, where $b$ is
the bit length of `prec`, until the image of the candidate set lies in
its interior, as Rump's theorem requires (each entry is tested with
`acb_contains_interior`). See [](algorithms/ball-arithmetic/eigenvalues).
:::
:::

::: function
int acb_mat_eig_simple_rump(acb_ptr E, acb_mat_t L, acb_mat_t R, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, slong prec)
int acb_mat_eig_simple_vdhoeven_mourrain(acb_ptr E, acb_mat_t L, acb_mat_t R, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, slong prec)
int acb_mat_eig_simple(acb_ptr E, acb_mat_t L, acb_mat_t R, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, slong prec)

Attempt to prove that the $n \times n$ matrix $A$ has $n$ simple
eigenvalues, given approximations of the eigenvalues and right
eigenvectors. Return 1 on success, writing pairwise disjoint enclosures
of the eigenvalues to `E` (in no particular order), and if `L` or `R` is
not `NULL`, enclosures of the corresponding left eigenvectors (rows of
`L`) and right eigenvectors (columns of `R`), normalised so that $L =
R^{-1}$, which gives a diagonalisation $LAR = D$. Return 0 otherwise.

::: details
The *rump* version certifies the eigenpairs one by one with
`acb_mat_eig_enclosure_rump` and fails as soon as an enclosure overlaps
a previous one; $L$ is obtained by inverting $R$; the cost is
$O(n^4)$. The *vdhoeven_mourrain* version certifies all eigenpairs at
once with the algorithm of [@HM2017] at a cost of $O(n^3)$, and is the
default. These functions do not try to compute clusters; use
`acb_mat_eig_multiple` for matrices that may have multiple or
clustered eigenvalues.
:::
:::

::: function
int acb_mat_eig_multiple_rump(acb_ptr E, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, slong prec)
int acb_mat_eig_multiple(acb_ptr E, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, slong prec)

Compute enclosures of all $n$ eigenvalues of $A$ with multiplicities,
given approximations of the eigenvalues and right eigenvectors. On
success (return value 1), the $n$ entries of `E` are either disjoint or
identical, identical entries are consecutive, and a run of $k$
identical entries is a ball that contains exactly $k$ eigenvalues
(a cluster that could not be separated at this precision) and is
disjoint from the other $n - k$ eigenvalues. Return 0 otherwise.

::: details
The *rump* version groups the approximate eigenvalues that lie within
$M 2^{-3\mathit{prec}/4}$ of each other, where $M$ is the largest
approximate eigenvalue (and, in a second attempt, within $M
2^{-\mathit{prec}/2}$), validates each group with
`acb_mat_eig_enclosure_rump` and checks that the enclosures are
disjoint; for $m$ clusters it costs $O(m n^3)$. The default
`acb_mat_eig_multiple` first tries `acb_mat_eig_simple_vdhoeven_mourrain`
and falls back to the *rump* version.
:::
:::

## Internal functions {#internal}

::: function
void acb_mat_swap_rows(acb_mat_t mat, slong * perm, slong r, slong s)

Swaps rows $r$ and $s$ of `mat` and, if `perm` is not `NULL`, the
entries $r$ and $s$ of `perm`.
:::

::: function
slong acb_mat_find_pivot_partial(const acb_mat_t mat, slong start_row, slong end_row, slong c)

Returns the index $i$ with $\mathit{start\_row} \le i <
\mathit{end\_row}$ of an entry in column $c$ that does not contain zero
and has the largest midpoint (comparing the larger of the absolute values
of the real and imaginary parts), or $-1$ if every entry in that
range contains zero.
:::
