---
title: fmpz_poly_mat.h
header: fmpz_poly_mat.h
description: Dense matrices of polynomials over the integers, with multiplication by Kronecker substitution, determinants by fraction-free elimination or interpolation, and exact solving with polynomial denominators.
category: Matrices
order: 6
---

This module implements dense matrices whose entries are polynomials in
$\Z[x]$ (`fmpz_poly_t`). Besides arithmetic, it provides multiplication
by Kronecker substitution, which turns a product of polynomial matrices
into one product of integer matrices handled by `fmpz_mat_mul`, and
linear algebra over the integral domain $\Z[x]$ by fraction-free
elimination: determinants (also by evaluation and interpolation), ranks,
inverses, solutions of linear systems and nullspaces, with results
expressed over a common polynomial denominator. The algorithms are
described in [](algorithms/linear-algebra/polynomial-matrices) and
[](algorithms/linear-algebra/fraction-free).

::: summary algorithms/linear-algebra/polynomial-matrices :::

## Types and conventions {#types}

::: type
fmpz_poly_mat_struct
fmpz_poly_mat_t

An `fmpz_poly_mat_struct` has the fields `entries` (a pointer to the
entry in row 0, column 0), `r` and `c` (the numbers of rows and columns)
and `stride` (the distance between the starts of consecutive rows). An
`fmpz_poly_mat_t` is an array of one `fmpz_poly_mat_struct`.
:::

Matrices with zero rows or columns are allowed. Outputs must be
initialised with the dimensions of the result; dimensions are not
checked. Unless stated otherwise, outputs may be aliased with inputs only
for entrywise operations.

## Memory management {#memory}

::: function
void fmpz_poly_mat_init(fmpz_poly_mat_t mat, slong rows, slong cols)

Initialises `mat` as a `rows` by `cols` matrix of zero polynomials.
Throws if $\mathit{rows} \cdot \mathit{cols}$ overflows an `slong`.
:::

::: function
void fmpz_poly_mat_init_set(fmpz_poly_mat_t mat, const fmpz_poly_mat_t src)

Initialises `mat` as a copy of `src`.
:::

::: function
void fmpz_poly_mat_clear(fmpz_poly_mat_t mat)

Frees the memory used by `mat`. Must not be called on a window.
:::

::: function
void fmpz_poly_mat_swap(fmpz_poly_mat_t mat1, fmpz_poly_mat_t mat2)
void fmpz_poly_mat_swap_entrywise(fmpz_poly_mat_t mat1, fmpz_poly_mat_t mat2)

Swaps two matrices in constant time by swapping their structs (not for
windows), respectively swaps the entries of two matrices of the same
dimensions one by one.
:::

## Entries and dimensions {#entries}

::: function
fmpz_poly_struct * fmpz_poly_mat_entry(const fmpz_poly_mat_t mat, slong i, slong j)

Returns a pointer to the entry in row $i$ and column $j$, which can be
used with the functions of [](fmpz_poly).
:::

::: function
slong fmpz_poly_mat_nrows(const fmpz_poly_mat_t mat)
slong fmpz_poly_mat_ncols(const fmpz_poly_mat_t mat)

Returns the number of rows, respectively columns, of `mat`.
:::

::: function
int fmpz_poly_mat_is_empty(const fmpz_poly_mat_t mat)
int fmpz_poly_mat_is_square(const fmpz_poly_mat_t mat)

Returns 1 if `mat` has no rows or no columns, respectively as many rows as
columns, and 0 otherwise.
:::

::: function
slong fmpz_poly_mat_max_bits(const fmpz_poly_mat_t A)

Returns $\pm b$, where $b$ is the largest bit size of the absolute value
of a coefficient of an entry of $A$, negative if some coefficient is
negative.
:::

::: function
slong fmpz_poly_mat_max_length(const fmpz_poly_mat_t A)

Returns the largest length of an entry of $A$ (0 if all entries are
zero).
:::

## Windows and concatenation {#windows}

::: function
void fmpz_poly_mat_window_init(fmpz_poly_mat_t window, const fmpz_poly_mat_t mat, slong r1, slong c1, slong r2, slong c2)

Initialises `window` as the submatrix of `mat` formed by rows $r_1$ to
$r_2 - 1$ and columns $c_1$ to $c_2 - 1$, sharing its entries with `mat`.
The bounds are checked only by `FLINT_ASSERT`.
:::

::: function
void fmpz_poly_mat_window_clear(fmpz_poly_mat_t window)

Releases a window. This does nothing at present.
:::

::: function
void fmpz_poly_mat_concat_horizontal(fmpz_poly_mat_t res, const fmpz_poly_mat_t mat1, const fmpz_poly_mat_t mat2)
void fmpz_poly_mat_concat_vertical(fmpz_poly_mat_t res, const fmpz_poly_mat_t mat1, const fmpz_poly_mat_t mat2)

Sets `res` to `mat1` placed to the left of, respectively above, `mat2`,
whose dimensions must be compatible; `res` must have the dimensions of
the result.
:::

## Input and output {#io}

::: function
void fmpz_poly_mat_print(const fmpz_poly_mat_t mat, const char * x)

Prints `mat` to `stdout`: a header `<r x c matrix over Z[x]>` (with the
variable name `x`), then each row in square brackets with the entries
printed by `fmpz_poly_print_pretty` and separated by commas, then an
empty line.
:::

## Random matrices {#random}

::: function
void fmpz_poly_mat_randtest(fmpz_poly_mat_t mat, flint_rand_t state, slong len, flint_bitcnt_t bits)
void fmpz_poly_mat_randtest_unsigned(fmpz_poly_mat_t mat, flint_rand_t state, slong len, flint_bitcnt_t bits)

Sets each entry to a random polynomial of length up to `len` with signed,
respectively nonnegative, coefficients of up to `bits` bits.
:::

::: function
void fmpz_poly_mat_randtest_sparse(fmpz_poly_mat_t A, flint_rand_t state, slong len, flint_bitcnt_t bits, float density)

Sets each entry, with probability `density`, to a random nonzero
polynomial of length up to `len` with coefficients of up to `bits` bits,
and to zero otherwise.
:::

## Comparison and special matrices {#comparison}

::: function
int fmpz_poly_mat_equal(const fmpz_poly_mat_t mat1, const fmpz_poly_mat_t mat2)

Returns 1 if the matrices have the same dimensions and entries, and 0
otherwise.
:::

::: function
int fmpz_poly_mat_is_zero(const fmpz_poly_mat_t mat)
int fmpz_poly_mat_is_one(const fmpz_poly_mat_t mat)

Returns 1 if `mat` is zero, respectively has the constant 1 on the main
diagonal and zeros elsewhere (it need not be square), and 0 otherwise.
:::

::: function
void fmpz_poly_mat_zero(fmpz_poly_mat_t mat)
void fmpz_poly_mat_one(fmpz_poly_mat_t mat)

Sets `mat` to the zero matrix, respectively to the matrix with the
constant 1 on the main diagonal and zeros elsewhere.
:::

## Basic arithmetic {#arithmetic}

::: function
void fmpz_poly_mat_set(fmpz_poly_mat_t mat1, const fmpz_poly_mat_t mat2)

Copies `mat2` to `mat1`, which must have the same dimensions.
:::

::: function
void fmpz_poly_mat_transpose(fmpz_poly_mat_t B, const fmpz_poly_mat_t A)

Sets $B$ to the transpose of $A$; $B$ may be $A$ if it is square.
:::

::: function
void fmpz_poly_mat_truncate(fmpz_poly_mat_t A, slong len)

Truncates every entry of $A$ to length `len` (reduces it modulo
$x^{\mathit{len}}$), in place.
:::

::: function
void fmpz_poly_mat_evaluate_fmpz(fmpz_mat_t B, const fmpz_poly_mat_t A, const fmpz_t x)

Sets the integer matrix $B$ to $A$ evaluated at $x$.
:::

::: function
void fmpz_poly_mat_add(fmpz_poly_mat_t C, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B)
void fmpz_poly_mat_sub(fmpz_poly_mat_t C, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B)
void fmpz_poly_mat_neg(fmpz_poly_mat_t B, const fmpz_poly_mat_t A)

Sets $C = A + B$, $C = A - B$, respectively $B = -A$. Aliasing is allowed.
:::

::: function
void fmpz_poly_mat_scalar_mul_fmpz_poly(fmpz_poly_mat_t B, const fmpz_poly_mat_t A, const fmpz_poly_t c)
void fmpz_poly_mat_scalar_mul_fmpz(fmpz_poly_mat_t B, const fmpz_poly_mat_t A, const fmpz_t c)

Sets $B = cA$ for a polynomial, respectively integer, $c$. Aliasing is
allowed.
:::

::: function
void fmpz_poly_mat_trace(fmpz_poly_t trace, const fmpz_poly_mat_t mat)

Sets `trace` to the sum of the diagonal entries of the square matrix
`mat`. Throws if `mat` is not square.
:::

## Multiplication {#mul}

::: function
void fmpz_poly_mat_mul(fmpz_poly_mat_t C, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B)
void fmpz_poly_mat_mul_classical(fmpz_poly_mat_t C, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B)
void fmpz_poly_mat_mul_KS(fmpz_poly_mat_t C, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B)

Sets $C = AB$. $C$ may be $A$ or $B$. `fmpz_poly_mat_mul_classical`
computes each entry as a sum of polynomial products;
`fmpz_poly_mat_mul_KS` packs the entries into integers by Kronecker
substitution and multiplies the integer matrices with `fmpz_mat_mul`.
`fmpz_poly_mat_mul` uses the classical algorithm if a dimension is less
than 8 and Kronecker substitution otherwise. See
[](algorithms/linear-algebra/polynomial-matrices#mul).
:::

::: function
void fmpz_poly_mat_mullow(fmpz_poly_mat_t C, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B, slong len)

Sets $C$ to $AB$ with every entry truncated to length `len`, computing
only the needed coefficients of each polynomial product. $C$ may be $A$
or $B$.
:::

::: function
void fmpz_poly_mat_sqr(fmpz_poly_mat_t B, const fmpz_poly_mat_t A)
void fmpz_poly_mat_sqr_classical(fmpz_poly_mat_t B, const fmpz_poly_mat_t A)
void fmpz_poly_mat_sqr_KS(fmpz_poly_mat_t B, const fmpz_poly_mat_t A)

Sets $B = A^2$ for a square matrix $A$; $B$ may be $A$. The classical
version uses five polynomial products for $2 \times 2$ matrices;
`fmpz_poly_mat_sqr` uses it for dimension less than 8 and Kronecker
substitution otherwise.
:::

::: function
void fmpz_poly_mat_sqrlow(fmpz_poly_mat_t B, const fmpz_poly_mat_t A, slong len)

Sets $B$ to $A^2$ with every entry truncated to length `len`.
:::

::: function
void fmpz_poly_mat_pow(fmpz_poly_mat_t B, const fmpz_poly_mat_t A, ulong exp)

Sets $B = A^{\mathit{exp}}$ for a square matrix $A$ by binary powering.
$A^0$ is the identity. $B$ may be $A$ but must not be a window.
:::

::: function
void fmpz_poly_mat_pow_trunc(fmpz_poly_mat_t B, const fmpz_poly_mat_t A, ulong exp, slong len)

Sets $B$ to $A^{\mathit{exp}}$ with every entry truncated to length
`len`, truncating after every multiplication. $B$ may be $A$ but must not
be a window.
:::

::: function
void fmpz_poly_mat_prod(fmpz_poly_mat_t res, fmpz_poly_mat_t * const factors, slong n)

Sets `res` to the product of the $n$ matrices `factors[0]`, ...,
`factors[n-1]`, all of the same square dimension, computed by a balanced
product tree. The empty product is the identity.
:::

## Row reduction {#rref}

::: function
slong fmpz_poly_mat_find_pivot_any(const fmpz_poly_mat_t mat, slong start_row, slong end_row, slong c)

Returns the index of the first row $r$ with `start_row` $\le r <$
`end_row` whose entry in column $c$ is nonzero, or $-1$ if there is none.
:::

::: function
slong fmpz_poly_mat_find_pivot_partial(const fmpz_poly_mat_t mat, slong start_row, slong end_row, slong c)

Returns the index of a row $r$ with `start_row` $\le r <$ `end_row`
whose entry in column $c$ is nonzero, choosing one of least length and,
among those, with coefficients of least bit size; returns $-1$ if all
entries are zero. Small pivots keep the entries of fraction-free
elimination small.
:::

::: function
slong fmpz_poly_mat_fflu(fmpz_poly_mat_t B, fmpz_poly_t den, slong * perm, const fmpz_poly_mat_t A, int rank_check)

Computes a fraction-free LU decomposition of $A$ in $B$ and returns the
rank $r$ of $A$ over $\Q(x)$. `den` is set to the last pivot, a nonzero
$r \times r$ minor of $A$ (1 if $r = 0$). Row exchanges are applied to
`perm` if it is not `NULL`. Pivots are chosen with
`fmpz_poly_mat_find_pivot_partial`. If `rank_check` is set, the function
returns 0 with `den` = 0 as soon as a column without pivot is found. $B$
may be $A$.

::: details
The layout is that of `fmpz_mat_fflu`: the fraction-free echelon form on
and above the pivots, and the multipliers below them; see
[](algorithms/linear-algebra/fraction-free).
:::
:::

::: function
slong fmpz_poly_mat_rref(fmpz_poly_mat_t B, fmpz_poly_t den, const fmpz_poly_mat_t A)

Sets ($B$, `den`) to the reduced row echelon form of $A$ over $\Q(x)$ and
returns the rank: $B/\mathit{den}$ is in reduced row echelon form, with
`den` on the pivot positions. Uses fraction-free elimination followed by
fraction-free back substitution. $B$ may be $A$.
:::

## Determinant, rank and trace {#det}

::: function
void fmpz_poly_mat_det(fmpz_poly_t det, const fmpz_poly_mat_t A)
void fmpz_poly_mat_det_fflu(fmpz_poly_t det, const fmpz_poly_mat_t A)
void fmpz_poly_mat_det_interpolate(fmpz_poly_t det, const fmpz_poly_mat_t A)

Sets `det` to the determinant of the square matrix $A$.
`fmpz_poly_mat_det_fflu` uses fraction-free elimination;
`fmpz_poly_mat_det_interpolate` evaluates $A$ at $L = n(\ell - 1) + 1$
consecutive integers around 0, where $\ell$ is the largest length of an
entry, computes the integer determinants with `fmpz_mat_det`, and
interpolates. `fmpz_poly_mat_det` uses explicit formulas for $n \le 2$,
fraction-free elimination for $n < 15$ and interpolation otherwise.
:::

::: function
slong fmpz_poly_mat_rank(const fmpz_poly_mat_t A)

Returns the rank of $A$ over $\Q(x)$, computed by fraction-free
elimination.
:::

## Solving and inversion {#solve}

::: function
int fmpz_poly_mat_solve(fmpz_poly_mat_t X, fmpz_poly_t den, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B)
int fmpz_poly_mat_solve_fflu(fmpz_poly_mat_t X, fmpz_poly_t den, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B)

For a square $A$ of dimension $n$ and an $n \times k$ matrix $B$, sets
($X$, `den`) with $AX = \mathit{den} \cdot B$ and $\mathit{den} =
\det A$, and returns 1 if $A$ is nonsingular; returns 0 with `den` = 0 if
$A$ is singular. Both use fraction-free elimination
(`fmpz_poly_mat_solve` calls `fmpz_poly_mat_solve_fflu`). $X$ may be $A$
but not $B$.
:::

::: function
void fmpz_poly_mat_solve_fflu_precomp(fmpz_poly_mat_t X, const slong * perm, const fmpz_poly_mat_t FFLU, const fmpz_poly_mat_t B)

Solves $AX = \mathit{den} \cdot B$ given the output `FFLU`, `den`, `perm`
of `fmpz_poly_mat_fflu` for a nonsingular $A$, with fraction-free forward
and back substitution. The sign of the permutation is not taken into
account. `perm` must not be `NULL`, and $X$ must not be $B$; both raise an
exception.
:::

::: function
int fmpz_poly_mat_inv(fmpz_poly_mat_t Ainv, fmpz_poly_t den, const fmpz_poly_mat_t A)

Sets (`Ainv`, `den`) to the inverse of the square matrix $A$, so that
$A \cdot \mathit{Ainv} = \mathit{den} \cdot I$ with $\mathit{den} =
\det A$, and returns 1 if $A$ is nonsingular; returns 0 if it is
singular. `Ainv` may be $A$.
:::

::: function
slong fmpz_poly_mat_nullspace(fmpz_poly_mat_t res, const fmpz_poly_mat_t mat)

Sets the first $n - r$ columns of the $n \times n$ matrix `res` to a basis
of the right nullspace of the $m \times n$ matrix `mat` over $\Q(x)$,
with polynomial entries, and returns $n - r$; the other columns are set
to zero. The basis is read off from `fmpz_poly_mat_rref`.
:::
