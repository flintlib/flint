---
title: nmod_poly_mat.h
header: nmod_poly_mat.h
description: Dense matrices of polynomials over $\Z/n\Z$ for word-size $n$, with multiplication by Kronecker substitution, evaluation and interpolation or the small-prime FFT, determinants, and exact solving by fraction-free elimination.
category: Matrices
order: 7
---

This module implements dense matrices whose entries are polynomials over
$\Z/n\Z$ (`nmod_poly_t`) for a word-size modulus $n$. Such matrices can
also be seen as polynomials with matrix coefficients, and the module
converts between the two views. Multiplication uses the classical
algorithm, Kronecker substitution into an integer matrix product,
evaluation and interpolation at many points (for a large prime modulus),
or the small-prime FFT with each entry transformed once. Determinants,
ranks, inverses, solutions and nullspaces are computed by fraction-free
elimination over $(\Z/n\Z)[x]$ (or, for determinants, by evaluation and
interpolation), which requires $n$ to be prime. The algorithms are
described in [](algorithms/linear-algebra/polynomial-matrices).

::: summary algorithms/linear-algebra/polynomial-matrices :::

## Types and conventions {#types}

::: type
nmod_poly_mat_struct
nmod_poly_mat_t

An `nmod_poly_mat_struct` has the fields `entries` (a pointer to the
entry in row 0, column 0), `r` and `c` (the numbers of rows and columns),
`stride` (the distance between the starts of consecutive rows) and
`modulus`. An `nmod_poly_mat_t` is an array of one
`nmod_poly_mat_struct`.
:::

All matrices passed to one function must have the same modulus. The
entries are `nmod_poly_t` with that modulus and can be manipulated with
the functions of [](nmod_poly). Outputs must be initialised with the
dimensions of the result. The arithmetic functions accept any modulus
$n \ge 1$ (except `nmod_poly_mat_mul_interpolate`); elimination,
determinants for dimension at least 3, ranks, solving, inverses and
nullspaces divide by polynomials and require $n$ to be prime.

## Memory management {#memory}

::: function
void nmod_poly_mat_init(nmod_poly_mat_t mat, slong rows, slong cols, ulong n)

Initialises `mat` as a `rows` by `cols` matrix of zero polynomials modulo
$n$. Throws if $\mathit{rows} \cdot \mathit{cols}$ overflows an `slong`.
:::

::: function
void nmod_poly_mat_init_set(nmod_poly_mat_t mat, const nmod_poly_mat_t src)

Initialises `mat` as a copy of `src`, with the same modulus.
:::

::: function
void nmod_poly_mat_clear(nmod_poly_mat_t mat)

Frees the memory used by `mat`. Must not be called on a window.
:::

::: function
void nmod_poly_mat_swap(nmod_poly_mat_t mat1, nmod_poly_mat_t mat2)
void nmod_poly_mat_swap_entrywise(nmod_poly_mat_t mat1, nmod_poly_mat_t mat2)

Swaps two matrices in constant time by swapping their structs (not for
windows), respectively swaps the entries of two matrices of the same
dimensions one by one.
:::

## Entries and dimensions {#entries}

::: function
nmod_poly_struct * nmod_poly_mat_entry(const nmod_poly_mat_t mat, slong i, slong j)

Returns a pointer to the entry in row $i$ and column $j$.
:::

::: function
slong nmod_poly_mat_nrows(const nmod_poly_mat_t mat)
slong nmod_poly_mat_ncols(const nmod_poly_mat_t mat)

Returns the number of rows, respectively columns, of `mat`.
:::

::: function
ulong nmod_poly_mat_modulus(const nmod_poly_mat_t mat)

Returns the modulus $n$ of `mat`.
:::

::: function
int nmod_poly_mat_is_empty(const nmod_poly_mat_t mat)
int nmod_poly_mat_is_square(const nmod_poly_mat_t mat)

Returns 1 if `mat` has no rows or no columns, respectively as many rows as
columns, and 0 otherwise.
:::

::: function
slong nmod_poly_mat_max_length(const nmod_poly_mat_t A)
slong nmod_poly_mat_degree(const nmod_poly_mat_t pmat)

Returns the largest length of an entry (0 if all entries are zero),
respectively the largest degree, which is the largest length minus one
($-1$ for the zero matrix).
:::

## Polynomial-matrix view {#coefficients}

::: function
void nmod_poly_mat_get_coeff_mat(nmod_mat_t coeff, const nmod_poly_mat_t pmat, slong deg)
void nmod_poly_mat_set_coeff_mat(nmod_poly_mat_t pmat, const nmod_mat_t coeff, slong deg)

Sets the matrix `coeff` to the coefficients of $x^{\mathit{deg}}$ of the
entries of `pmat`, respectively sets these coefficients of `pmat` to the
entries of `coeff`. The dimensions must agree.
:::

::: function
void nmod_poly_mat_set_nmod_mat(nmod_poly_mat_t pmat, const nmod_mat_t cmat)

Sets `pmat` to the constant matrix `cmat`, which must have the same
dimensions.
:::

::: function
int nmod_poly_mat_equal_nmod_mat(const nmod_poly_mat_t pmat, const nmod_mat_t cmat)

Returns 1 if `pmat` is the constant matrix `cmat` (same dimensions, all
entries constant and equal), and 0 otherwise.
:::

::: function
void nmod_poly_mat_set_trunc(nmod_poly_mat_t res, const nmod_poly_mat_t pmat, long len)
void nmod_poly_mat_truncate(nmod_poly_mat_t pmat, long len)

Sets `res` to `pmat` with every entry truncated to length `len`
(reduced modulo $x^{\mathit{len}}$), respectively truncates `pmat` in
place.
:::

::: function
void nmod_poly_mat_shift_left(nmod_poly_mat_t res, const nmod_poly_mat_t pmat, slong k)
void nmod_poly_mat_shift_right(nmod_poly_mat_t res, const nmod_poly_mat_t pmat, slong k)

Sets `res` to `pmat` multiplied by $x^k$, respectively divided by $x^k$
with the coefficients of degree less than $k$ discarded. Aliasing is
allowed.
:::

## Windows and concatenation {#windows}

::: function
void nmod_poly_mat_window_init(nmod_poly_mat_t window, const nmod_poly_mat_t mat, slong r1, slong c1, slong r2, slong c2)

Initialises `window` as the submatrix of `mat` formed by rows $r_1$ to
$r_2 - 1$ and columns $c_1$ to $c_2 - 1$, sharing its entries and modulus
with `mat`.
:::

::: function
void nmod_poly_mat_window_clear(nmod_poly_mat_t window)

Releases a window. This does nothing at present.
:::

::: function
void nmod_poly_mat_concat_horizontal(nmod_poly_mat_t res, const nmod_poly_mat_t mat1, const nmod_poly_mat_t mat2)
void nmod_poly_mat_concat_vertical(nmod_poly_mat_t res, const nmod_poly_mat_t mat1, const nmod_poly_mat_t mat2)

Sets `res` to `mat1` placed to the left of, respectively above, `mat2`,
whose dimensions must be compatible.
:::

## Input and output {#io}

::: function
void nmod_poly_mat_print(const nmod_poly_mat_t mat, const char * x)

Prints `mat` to `stdout`: a header `<r x c matrix over Z/nZ[x]>` (with
the variable name `x`), then each row in square brackets with the entries
printed by `nmod_poly_print` (length, modulus and coefficients) and
separated by commas, then an empty line.
:::

## Random matrices {#random}

::: function
void nmod_poly_mat_rand(nmod_poly_mat_t A, flint_rand_t state, slong len)
void nmod_poly_mat_randtest(nmod_poly_mat_t mat, flint_rand_t state, slong len)

Sets each entry to a random polynomial of length up to `len`, with
uniformly random coefficients (`nmod_poly_rand`), respectively with a bias
towards special values (`nmod_poly_randtest`).
:::

::: function
void nmod_poly_mat_randtest_sparse(nmod_poly_mat_t A, flint_rand_t state, slong len, float density)

Sets each entry, with probability `density`, to a random polynomial of
length up to `len`, and to zero otherwise.
:::

## Comparison and special matrices {#comparison}

::: function
int nmod_poly_mat_equal(const nmod_poly_mat_t mat1, const nmod_poly_mat_t mat2)

Returns 1 if the matrices have the same dimensions and entries, and 0
otherwise.
:::

::: function
int nmod_poly_mat_is_zero(const nmod_poly_mat_t mat)
int nmod_poly_mat_is_one(const nmod_poly_mat_t mat)

Returns 1 if `mat` is zero, respectively has the constant 1 on the main
diagonal and zeros elsewhere, and 0 otherwise.
:::

::: function
void nmod_poly_mat_zero(nmod_poly_mat_t mat)
void nmod_poly_mat_one(nmod_poly_mat_t mat)

Sets `mat` to the zero matrix, respectively to the matrix with the
constant 1 on the main diagonal and zeros elsewhere.
:::

## Basic arithmetic {#arithmetic}

::: function
void nmod_poly_mat_set(nmod_poly_mat_t mat1, const nmod_poly_mat_t mat2)

Copies `mat2` to `mat1`, which must have the same dimensions.
:::

::: function
void nmod_poly_mat_transpose(nmod_poly_mat_t B, const nmod_poly_mat_t A)

Sets $B$ to the transpose of $A$; $B$ may be $A$ if it is square. Throws
if the dimensions are incompatible.
:::

::: function
void nmod_poly_mat_evaluate_nmod(nmod_mat_t B, const nmod_poly_mat_t A, ulong x)

Sets the matrix $B$ over $\Z/n\Z$ to $A$ evaluated at $x$.
:::

::: function
void nmod_poly_mat_add(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B)
void nmod_poly_mat_sub(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B)
void nmod_poly_mat_neg(nmod_poly_mat_t B, const nmod_poly_mat_t A)

Sets $C = A + B$, $C = A - B$, respectively $B = -A$. Aliasing is allowed.
:::

::: function
void nmod_poly_mat_scalar_mul_nmod_poly(nmod_poly_mat_t B, const nmod_poly_mat_t A, const nmod_poly_t c)
void nmod_poly_mat_scalar_mul_nmod(nmod_poly_mat_t B, const nmod_poly_mat_t A, ulong c)

Sets $B = cA$ for a polynomial, respectively a reduced constant, $c$.
Aliasing is allowed.
:::

::: function
void nmod_poly_mat_trace(nmod_poly_t trace, const nmod_poly_mat_t mat)

Sets `trace` to the sum of the diagonal entries of the square matrix
`mat`.
:::

## Multiplication {#mul}

::: function
void nmod_poly_mat_mul(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B)

Sets $C = AB$, choosing between the methods below from the dimensions,
the lengths of the entries and the modulus. $C$ may be $A$ or $B$. See
[](algorithms/linear-algebra/polynomial-matrices#choice).
:::

::: function
void nmod_poly_mat_mul_classical(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B)
void nmod_poly_mat_mul_KS(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B)

Sets $C = AB$ by sums of polynomial products, respectively by Kronecker
substitution: the entries are packed into integers with
`nmod_poly_bit_pack`, multiplied with `fmpz_mat_mul` and unpacked. $C$
may be $A$ or $B$.
:::

::: function
void nmod_poly_mat_mul_interpolate(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B)

Sets $C = AB$ by evaluating $A$ and $B$ at $L = \ell_A + \ell_B - 1$
points (where $\ell_A$, $\ell_B$ are the largest lengths of the entries),
multiplying the $L$ pairs of matrices over $\Z/n\Z$ with `nmod_mat_mul`
and interpolating. Requires $n$ to be a prime with $n \ge L$; throws if
$n < L$.
:::

::: function
int nmod_poly_mat_mulmid_fft_small(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B, slong zl, slong zh)

Sets each entry of $C$ to the coefficients of degrees $z_l$ to $z_h - 1$
of the corresponding entry of $AB$ (the product shifted right by $z_l$
and truncated to length $z_h - z_l$), and returns 1; returns 0 if the
small-prime FFT is not available or no transform plan admits the size of
the accumulated products, in which case $C$ may be partially written.
$C$ may be $A$ or $B$.

::: details
Each entry of $A$ and $B$ is transformed once in the transformed
polynomial ring of `gr` (`gr_ctx_init_gr_poly_transformed_repr`), the
$k$ pointwise products of each output entry are accumulated in the
transformed domain, and only the requested window is reconstructed by
Chinese remaindering. The work is spread over threads. If the transforms
of all inputs do not fit in the storage budget of [](fft_small), the
product is computed in blocks of rows and columns. With $z_l = 0$ and
$z_h = \ell_A + \ell_B - 1$ this is the full product, which is how
`nmod_poly_mat_mul` uses it for long entries.
:::
:::

::: function
void nmod_poly_mat_sqr(nmod_poly_mat_t B, const nmod_poly_mat_t A)
void nmod_poly_mat_sqr_classical(nmod_poly_mat_t B, const nmod_poly_mat_t A)
void nmod_poly_mat_sqr_KS(nmod_poly_mat_t B, const nmod_poly_mat_t A)
void nmod_poly_mat_sqr_interpolate(nmod_poly_mat_t B, const nmod_poly_mat_t A)

Sets $B = A^2$ for a square matrix $A$; $B$ may be $A$. The classical
version uses five polynomial products for $2 \times 2$ matrices; the
others are the squaring versions of the multiplication methods.
`nmod_poly_mat_sqr` uses the classical algorithm for dimension less than
10. Above, it uses evaluation and interpolation for a prime modulus of
more than $\beta/4$ bits that is at least $2\ell - 1$, where $\ell$ is
the largest length of an entry, and dimension above $80 + \sqrt{\ell}$;
otherwise Kronecker substitution, or the classical algorithm for entries
of length above 128.
:::

::: function
void nmod_poly_mat_pow(nmod_poly_mat_t B, const nmod_poly_mat_t A, ulong exp)

Sets $B = A^{\mathit{exp}}$ for a square matrix $A$ by binary powering.
$A^0$ is the identity. $B$ may be $A$ but must not be a window.
:::

## Row reduction {#rref}

These functions require a prime modulus.

::: function
slong nmod_poly_mat_find_pivot_any(const nmod_poly_mat_t mat, slong start_row, slong end_row, slong c)
slong nmod_poly_mat_find_pivot_partial(const nmod_poly_mat_t mat, slong start_row, slong end_row, slong c)

Returns the index of a row $r$ with `start_row` $\le r <$ `end_row`
whose entry in column $c$ is nonzero, or $-1$ if there is none:
the first such row, respectively one whose entry has least length.
:::

::: function
slong nmod_poly_mat_fflu(nmod_poly_mat_t B, nmod_poly_t den, slong * perm, const nmod_poly_mat_t A, int rank_check)

Computes a fraction-free LU decomposition of $A$ in $B$ and returns the
rank $r$ of $A$ over the field of rational functions. `den` is set to the
last pivot, a nonzero $r \times r$ minor of $A$ (1 if $r = 0$). Row
exchanges are applied to `perm` if it is not `NULL`; pivots are chosen
with `nmod_poly_mat_find_pivot_partial`. If `rank_check` is set, the
function returns 0 with `den` = 0 as soon as a column without pivot is
found. $B$ may be $A$. See [](algorithms/linear-algebra/fraction-free).
:::

::: function
slong nmod_poly_mat_rref(nmod_poly_mat_t B, nmod_poly_t den, const nmod_poly_mat_t A)

Sets ($B$, `den`) to the reduced row echelon form of $A$ over the field
of rational functions and returns the rank: $B/\mathit{den}$ is in
reduced row echelon form, with `den` on the pivot positions. $B$ may be
$A$.
:::

## Determinant and rank {#det}

::: function
void nmod_poly_mat_det(nmod_poly_t det, const nmod_poly_mat_t A)
void nmod_poly_mat_det_fflu(nmod_poly_t det, const nmod_poly_mat_t A)
void nmod_poly_mat_det_interpolate(nmod_poly_t det, const nmod_poly_mat_t A)

Sets `det` to the determinant of the square matrix $A$.
`nmod_poly_mat_det_fflu` uses fraction-free elimination;
`nmod_poly_mat_det_interpolate` evaluates $A$ at the points
$0, 1, \ldots, L - 1$ with $L = n(\ell - 1) + 1$, where $\ell$ is the
largest length of an entry, computes the determinants with
`nmod_mat_det` and interpolates, falling back to fraction-free
elimination if $L$ exceeds the modulus. `nmod_poly_mat_det` uses
explicit formulas for $n \le 2$, fraction-free elimination for $n < 15$
and interpolation otherwise. The modulus must be prime for $n \ge 3$.
:::

::: function
slong nmod_poly_mat_rank(const nmod_poly_mat_t A)

Returns the rank of $A$ over the field of rational functions, computed
by fraction-free elimination. Requires a prime modulus.
:::

## Solving and inversion {#solve}

These functions require a prime modulus.

::: function
int nmod_poly_mat_solve(nmod_poly_mat_t X, nmod_poly_t den, const nmod_poly_mat_t A, const nmod_poly_mat_t B)
int nmod_poly_mat_solve_fflu(nmod_poly_mat_t X, nmod_poly_t den, const nmod_poly_mat_t A, const nmod_poly_mat_t B)

For a square $A$ of dimension $n$ and an $n \times k$ matrix $B$, sets
($X$, `den`) with $AX = \mathit{den} \cdot B$ and $\mathit{den} =
\det A$, and returns 1 if $A$ is nonsingular; returns 0 with `den` = 0 if
$A$ is singular. Both use fraction-free elimination. $X$ may be $A$ but
not $B$.
:::

::: function
void nmod_poly_mat_solve_fflu_precomp(nmod_poly_mat_t X, const slong * perm, const nmod_poly_mat_t FFLU, const nmod_poly_mat_t B)

Solves $AX = \mathit{den} \cdot B$ given the output `FFLU`, `den`, `perm`
of `nmod_poly_mat_fflu` for a nonsingular $A$, with fraction-free forward
and back substitution. The sign of the permutation is not taken into
account. `perm` must not be `NULL`, and $X$ must not be $B$.
:::

::: function
int nmod_poly_mat_inv(nmod_poly_mat_t Ainv, nmod_poly_t den, const nmod_poly_mat_t A)

Sets (`Ainv`, `den`) to the inverse of the square matrix $A$, so that
$A \cdot \mathit{Ainv} = \mathit{den} \cdot I$ with
$\mathit{den} = \det A$, and returns 1 if $A$ is nonsingular; returns 0
if it is singular. `Ainv` may be $A$.
:::

::: function
slong nmod_poly_mat_nullspace(nmod_poly_mat_t res, const nmod_poly_mat_t mat)

Sets the first $n - r$ columns of the $n \times n$ matrix `res` to a
basis of the right nullspace of the $m \times n$ matrix `mat` over the
field of rational functions, with polynomial entries, and returns
$n - r$; the other columns are set to zero.
:::
