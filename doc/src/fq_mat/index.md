---
title: fq_mat.h
header: fq_mat.h
template-headers: fq_mat_templates
description: Dense matrices over finite fields $\mathbb{F}_q$ of arbitrary characteristic, with Kronecker substitution multiplication and linear algebra through the generic `gr_mat` algorithms.
category: Finite fields
order: 7
---

This module implements dense matrices over a finite field $\mathbb{F}_q$
in the representation of [](fq). Multiplication packs the entries into
integers ([Kronecker substitution](algorithms/polynomials/kronecker)) and
multiplies integer matrices for all but small matrices. The linear
algebra functions (LU decomposition, row echelon forms, solving,
inversion, rank, nullspace, characteristic and minimal polynomial) are
thin wrappers around the generic implementations in [](gr_mat), which
call `fq_mat_mul` for the matrix products. The interface is generated
from `fq_mat_templates.h` and is the same for `fq_nmod_mat` and
`fq_zech_mat`. `fq_default_mat` provides the same functions for any
finite field.

::: summary algorithms/linear-algebra/gaussian-elimination :::

## Types and conventions {#types}

::: type
fq_mat_struct
fq_mat_t

A dense matrix with `r` rows and `c` columns. The entries are stored in
the array `entries` row by row, with `stride` elements between the
starts of consecutive rows (equal to `c` except for windows). `fq_mat_t`
is an array of one structure.
:::

- Matrices must be initialised with the dimensions they are used with;
  no function resizes its output, and dimensions are checked only by
  assertions.
- Unless stated otherwise, outputs may be aliased with inputs of the same
  dimensions.
- Entries are elements of $\mathbb{F}_q$ as in [](fq): reduced and
  normalised.

## Memory management {#memory}

::: function
void fq_mat_init(fq_mat_t mat, slong rows, slong cols, const fq_ctx_t ctx)

Initialises `mat` as a `rows` by `cols` zero matrix.
:::

::: function
void fq_mat_init_set(fq_mat_t mat, const fq_mat_t src, const fq_ctx_t ctx)

Initialises `mat` as a copy of `src`, with the same dimensions.
:::

::: function
void fq_mat_clear(fq_mat_t mat, const fq_ctx_t ctx)

Frees the memory used by `mat`.
:::

::: function
void fq_mat_window_init(fq_mat_t window, const fq_mat_t mat, slong r1, slong c1, slong r2, slong c2, const fq_ctx_t ctx)

Initialises `window` as the submatrix of `mat` formed by the rows
$r_1 \le i < r_2$ and the columns $c_1 \le j < c_2$. The entries are
shared with `mat`, not copied, so changes to one are visible in the
other. `mat` must not be cleared while the window is in use.
:::

::: function
void fq_mat_window_clear(fq_mat_t window, const fq_ctx_t ctx)

Clears a window. The entries, which belong to the underlying matrix, are
not freed.
:::

## Entries and dimensions {#entries}

::: function
fq_struct * fq_mat_entry(const fq_mat_t mat, slong i, slong j)

Returns a pointer to the entry in row $i$ and column $j$, indexed from
zero. No bounds checking is performed.
:::

::: function
void fq_mat_entry_set(fq_mat_t mat, slong i, slong j, const fq_t x, const fq_ctx_t ctx)

Sets the entry in row $i$ and column $j$ to $x$.
:::

::: function
slong fq_mat_nrows(const fq_mat_t mat, const fq_ctx_t ctx)
slong fq_mat_ncols(const fq_mat_t mat, const fq_ctx_t ctx)

Return the number of rows, respectively columns, of `mat`.
:::

## Assignment and permutations {#assignment}

::: function
void fq_mat_set(fq_mat_t mat1, const fq_mat_t mat2, const fq_ctx_t ctx)

Copies `mat2` into `mat1`, which must have the same dimensions.
:::

::: function
void fq_mat_swap(fq_mat_t mat1, fq_mat_t mat2, const fq_ctx_t ctx)

Swaps the two matrices by swapping their structures; the dimensions may
differ. Not for windows.
:::

::: function
void fq_mat_swap_entrywise(fq_mat_t mat1, fq_mat_t mat2, const fq_ctx_t ctx)

Swaps the entries of two matrices of the same dimensions one by one;
this also works for windows.
:::

::: function
void fq_mat_zero(fq_mat_t A, const fq_ctx_t ctx)
void fq_mat_one(fq_mat_t A, const fq_ctx_t ctx)

Set `A` to the zero matrix, respectively to the matrix with ones on the
main diagonal and zeros elsewhere (also for non-square `A`).
:::

::: function
void fq_mat_swap_rows(fq_mat_t mat, slong * perm, slong r, slong s, const fq_ctx_t ctx)
void fq_mat_swap_cols(fq_mat_t mat, slong * perm, slong r, slong s, const fq_ctx_t ctx)

Swap rows, respectively columns, $r$ and $s$ of `mat`. If `perm` is not
`NULL`, its entries $r$ and $s$ are swapped as well, so that it records
the permutation.
:::

::: function
void fq_mat_invert_rows(fq_mat_t mat, slong * perm, const fq_ctx_t ctx)
void fq_mat_invert_cols(fq_mat_t mat, slong * perm, const fq_ctx_t ctx)

Reverse the order of the rows, respectively columns, of `mat`, by
swapping row $i$ and row $r - 1 - i$ for $0 \le i < \lfloor r/2
\rfloor$ (and likewise for columns). If `perm` is not `NULL`, the same
permutation is applied to it.
:::

## Conversion and concatenation {#conversion}

::: function
void fq_mat_set_nmod_mat(fq_mat_t mat1, const nmod_mat_t mat2, const fq_ctx_t ctx)
void fq_mat_set_fmpz_mod_mat(fq_mat_t mat1, const fmpz_mod_mat_t mat2, const fq_ctx_t ctx)

Set `mat1` to the image of the matrix `mat2` over $\mathbb{F}_p$ (or over
$\Z/n\Z$, whose entries are mapped to $\mathbb{F}_p$ by reduction) in
$\mathbb{F}_q$. The dimensions must agree.
:::

::: function
void fq_mat_concat_vertical(fq_mat_t res, const fq_mat_t mat1, const fq_mat_t mat2, const fq_ctx_t ctx)
void fq_mat_concat_horizontal(fq_mat_t res, const fq_mat_t mat1, const fq_mat_t mat2, const fq_ctx_t ctx)

Set `res` to the block matrix $\begin{pmatrix} A \\ B \end{pmatrix}$,
respectively $(A \ B)$, for $A$ = `mat1` and $B$ = `mat2`. For vertical
concatenation of an $m \times n$ and a $k \times n$ matrix, `res` must be
$(m + k) \times n$; for horizontal concatenation of an $m \times n$ and
an $m \times k$ matrix, $m \times (n + k)$.
:::

## Comparison {#comparison}

::: function
int fq_mat_equal(const fq_mat_t mat1, const fq_mat_t mat2, const fq_ctx_t ctx)

Returns whether the matrices have the same dimensions and entries.
:::

::: function
int fq_mat_is_zero(const fq_mat_t mat, const fq_ctx_t ctx)
int fq_mat_is_one(const fq_mat_t mat, const fq_ctx_t ctx)

Return whether all entries are zero, respectively whether the diagonal
entries are one and the others zero. An empty matrix is both zero and
one.
:::

::: function
int fq_mat_is_empty(const fq_mat_t mat, const fq_ctx_t ctx)
int fq_mat_is_square(const fq_mat_t mat, const fq_ctx_t ctx)

Return whether `mat` has zero rows or zero columns, respectively whether
it has as many rows as columns.
:::

## Random matrices {#random}

::: function
void fq_mat_randtest(fq_mat_t mat, flint_rand_t state, const fq_ctx_t ctx)

Sets the entries of `mat` to random elements for testing (see
`fq_randtest`).
:::

::: function
int fq_mat_randpermdiag(fq_mat_t mat, flint_rand_t state, fq_struct * diag, slong n, const fq_ctx_t ctx)

Sets `mat` to a random row and column permutation of the diagonal matrix
whose first $n$ diagonal entries are `diag` and whose other entries are
zero; $n$ must not exceed the number of rows or columns. Returns the
parity of the permutation (0 for even, 1 for odd).
:::

::: function
void fq_mat_randrank(fq_mat_t mat, flint_rand_t state, slong rank, const fq_ctx_t ctx)

Sets `mat` to a random sparse matrix of the given rank, with exactly
`rank` nonzero entries, which are random nonzero elements.
`fq_mat_randops` turns it into a dense matrix of the same rank.
:::

::: function
void fq_mat_randops(fq_mat_t mat, flint_rand_t state, slong count, const fq_ctx_t ctx)

Applies `count` random elementary row or column operations (additions or
subtractions of a row or column to another) to `mat`. This preserves the
rank and, for square matrices, the determinant.
:::

::: function
void fq_mat_randtril(fq_mat_t mat, flint_rand_t state, int unit, const fq_ctx_t ctx)
void fq_mat_randtriu(fq_mat_t mat, flint_rand_t state, int unit, const fq_ctx_t ctx)

Set `mat` to a random lower, respectively upper, triangular matrix,
with ones on the diagonal if `unit` is nonzero and random nonzero
diagonal entries otherwise.
:::

## Input and output {#io}

::: function
int fq_mat_fprint(FILE * file, const fq_mat_t mat, const fq_ctx_t ctx)
int fq_mat_print(const fq_mat_t mat, const fq_ctx_t ctx)

Print `mat` to `file`, respectively `stdout`: the number of rows and
columns separated by a space, two spaces, and the entries in row-major
order in the format of `fq_fprint`, separated by spaces. Return a
positive value on success and a non-positive value on failure.
:::

::: function
int fq_mat_fprint_pretty(FILE * file, const fq_mat_t mat, const fq_ctx_t ctx)
int fq_mat_print_pretty(const fq_mat_t mat, const fq_ctx_t ctx)

Print `mat` to `file`, respectively `stdout`, as a bracketed list of
bracketed rows, with the entries in the format of `fq_fprint_pretty`.
Return a positive value on success.
:::

## Arithmetic {#arithmetic}

::: function
void fq_mat_transpose(fq_mat_t B, const fq_mat_t A, const fq_ctx_t ctx)

Sets $B$ to $A^T$. $B$ must have the transposed dimensions of $A$;
aliasing is allowed for square matrices.
:::

::: function
void fq_mat_add(fq_mat_t C, const fq_mat_t A, const fq_mat_t B, const fq_ctx_t ctx)
void fq_mat_sub(fq_mat_t C, const fq_mat_t A, const fq_mat_t B, const fq_ctx_t ctx)

Set $C$ to $A + B$, respectively $A - B$. All three matrices must have
the same dimensions.
:::

::: function
void fq_mat_neg(fq_mat_t B, const fq_mat_t A, const fq_ctx_t ctx)

Sets $B$ to $-A$.
:::

::: function
void fq_mat_scalar_mul(fq_mat_t B, const fq_mat_t A, const fq_t c, const fq_ctx_t ctx)
void fq_mat_scalar_mul_fmpz(fq_mat_t B, const fq_mat_t A, const fmpz_t c, const fq_ctx_t ctx)
void fq_mat_scalar_mul_si(fq_mat_t B, const fq_mat_t A, const slong c, const fq_ctx_t ctx)
void fq_mat_scalar_mul_ui(fq_mat_t B, const fq_mat_t A, const ulong c, const fq_ctx_t ctx)

Set $B$ to $cA$, where $c$ is an element of $\mathbb{F}_q$ or an integer
mapped to $\mathbb{F}_p$.
:::

::: function
void fq_mat_mul(fq_mat_t C, const fq_mat_t A, const fq_mat_t B, const fq_ctx_t ctx)

Sets $C = AB$, where $A$ is $m \times k$, $B$ is $k \times n$ and $C$ is
$m \times n$. Aliasing is allowed. The function uses `fq_mat_mul_KS` if
`FQ_MAT_MUL_KS_CUTOFF` returns 1 for $m$ and $n$, and
`fq_mat_mul_classical` otherwise.
:::

::: function
int FQ_MAT_MUL_KS_CUTOFF(slong r, slong c, const fq_ctx_t ctx)

Returns 1 if `fq_mat_mul` should use Kronecker substitution for an
$r \times k$ times $k \times c$ product, which is the case if
$5 \min(r, c) > 8d + 29$. (`FQ_NMOD_MAT_MUL_KS_CUTOFF` and
`FQ_ZECH_MAT_MUL_KS_CUTOFF` are the corresponding tuned functions of the
other representations.)
:::

::: function
void fq_mat_mul_classical(fq_mat_t C, const fq_mat_t A, const fq_mat_t B, const fq_ctx_t ctx)

Sets $C = AB$ by computing each entry as a dot product with
`_fq_vec_dot`, which reduces modulo the defining polynomial only once per
entry. Aliasing is allowed.
:::

::: function
void fq_mat_mul_KS(fq_mat_t C, const fq_mat_t A, const fq_mat_t B, const fq_ctx_t ctx)

Sets $C = AB$ by Kronecker substitution. $C$ must not be aliased with
$A$ or $B$.

::: details
Every entry, a polynomial of degree less than $d$ with coefficients in
$[0, p)$, is packed into an integer with `fq_bit_pack`, using fields of
$b = \operatorname{bits}((p - 1)^2 k d) + 1$ bits for an inner dimension
$k$, enough to hold every coefficient of a sum of $k$ products of two
such polynomials. The integer matrices are multiplied with
`fmpz_mat_mul`, which uses multimodular and BLAS-like methods (see
[](algorithms/linear-algebra/matrix-multiplication)), and the entries of
the product are unpacked with `fq_bit_unpack`, which reduces modulo $p$
and the defining polynomial. One reduction per entry is needed, and the
product of the packed matrices has entries of about $2db$ bits.
:::
:::

::: function
void fq_mat_submul(fq_mat_t D, const fq_mat_t C, const fq_mat_t A, const fq_mat_t B, const fq_ctx_t ctx)

Sets $D = C - AB$. $D$ may be aliased with $C$, $A$ or $B$.
:::

::: function
void fq_mat_mul_vec(fq_struct * c, const fq_mat_t A, const fq_struct * b, slong blen, const fq_ctx_t ctx)
void fq_mat_mul_vec_ptr(fq_struct * const * c, const fq_mat_t A, const fq_struct * const * b, slong blen, const fq_ctx_t ctx)

Set the vector $c$ to the matrix-vector product $Ab$. The vector
`(b, blen)` is truncated or extended by zeros to the number of columns
of $A$, and $c$ receives as many entries as $A$ has rows. The `_ptr`
version takes arrays of pointers to the entries. $c$ must not alias $b$.
:::

::: function
void fq_mat_vec_mul(fq_struct * c, const fq_struct * a, slong alen, const fq_mat_t B, const fq_ctx_t ctx)
void fq_mat_vec_mul_ptr(fq_struct * const * c, const fq_struct * const * a, slong alen, const fq_mat_t B, const fq_ctx_t ctx)

Set the vector $c$ to the vector-matrix product $aB$. The vector
`(a, alen)` is truncated or extended by zeros to the number of rows of
$B$, and $c$ receives as many entries as $B$ has columns. $c$ must not
alias $a$.
:::

## Gaussian elimination {#elimination}

These functions use recursive (block) LU decomposition and triangular
solving for large matrices
([](algorithms/linear-algebra/gaussian-elimination)), so that their cost
is dominated by matrix multiplications with `fq_mat_mul`.

::: function
slong fq_mat_lu(slong * P, fq_mat_t A, int rank_check, const fq_ctx_t ctx)

Computes a generalised LU decomposition $PA = LU$ of $A$ in place, and
returns the rank of $A$. `P` must have room for as many entries as $A$
has rows; it receives the row permutation.

::: details
If $A$ is a nonsingular square matrix, it is overwritten with the unit
lower triangular $L$ (without its diagonal) and the upper triangular
$U$. If $A$ has rank $r$, then $U$ is in row echelon form with $r$
nonzero rows, and $L$ is lower triangular with $r$ columns and implicit
ones on its diagonal. If `rank_check` is nonzero and $A$ is found to be
rank-deficient, the function stops early and returns 0, leaving $A$ in an
undefined state.
:::
:::

::: function
slong fq_mat_rref(fq_mat_t B, const fq_mat_t A, const fq_ctx_t ctx)

Sets $B$ to the reduced row echelon form of $A$ and returns the rank of
$A$. The form is computed from an LU decomposition. $B$ may be aliased
with $A$. If $A$ is zero, $B$
is set to zero and 0 is returned.
:::

::: function
slong fq_mat_reduce_row(fq_mat_t A, slong * P, slong * L, slong m, const fq_ctx_t ctx)

Reduces row $m$ of $A$ against the earlier rows, which must already be in
reduced form (pivot entries 1, not necessarily in order), and returns
the column of the new pivot, or $-1$ if row $m$ becomes zero. `P[i]` is
the row with its pivot in column $i$, or $-1$; it is updated when a new
pivot is found, and the new pivot row is scaled to make the pivot 1.
`L[i]` is the number of possibly nonzero columns of row $i$, which allows
faster reduction of matrices that are zero on the right; the entries of
`L` must be nondecreasing, and can all be set to the number of columns.
The arithmetic is done without intermediate reductions modulo the
defining polynomial.
:::

::: function
slong fq_mat_rank(const fq_mat_t A, const fq_ctx_t ctx)

Returns the rank of $A$, computed by LU decomposition.
:::

::: function
int fq_mat_inv(fq_mat_t B, fq_mat_t A, const fq_ctx_t ctx)

If the square matrix $A$ is invertible, sets $B = A^{-1}$ and returns 1;
otherwise returns 0 and the entries of $B$ are undefined. $A$ is not
modified (although not declared `const`), and $B$ may be aliased with
$A$.
:::

::: function
void fq_mat_solve_tril(fq_mat_t X, const fq_mat_t L, const fq_mat_t B, int unit, const fq_ctx_t ctx)
void fq_mat_solve_triu(fq_mat_t X, const fq_mat_t U, const fq_mat_t B, int unit, const fq_ctx_t ctx)

Set $X = L^{-1} B$, respectively $X = U^{-1} B$, for a square lower,
respectively upper, triangular matrix with nonzero diagonal entries. If
`unit` is nonzero, the diagonal is assumed to consist of ones and is not
read. Only the triangular part of the matrix is read. Aliasing of $X$ and
$B$ is allowed. Classical substitution is used if $B$ has fewer than 64
rows or columns, and block recursion otherwise.
:::

::: function
int fq_mat_solve(fq_mat_t X, const fq_mat_t A, const fq_mat_t C, const fq_ctx_t ctx)

Solves $AX = C$ for a square matrix $A$. Returns 1 and sets $X$ if $A$ is
nonsingular; otherwise returns 0 and $X$ is undefined.
:::

::: function
int fq_mat_can_solve(fq_mat_t X, const fq_mat_t A, const fq_mat_t B, const fq_ctx_t ctx)

Solves $AX = B$ for a matrix $A$ of any shape: if a solution exists, sets
$X$ to one and returns 1, otherwise returns 0 and $X$ is undefined. For
an $m \times n$ matrix $A$ and an $m \times k$ matrix $B$, $X$ must be
$n \times k$.
:::

::: function
slong fq_mat_nullspace(fq_mat_t X, const fq_mat_t A, const fq_ctx_t ctx)

Computes a basis of the right nullspace $\{v : Av = 0\}$ of $A$ and
returns its dimension $k$ (the nullity). $X$ must have as many rows and
at least as many columns as $A$ has columns; its first $k$ columns
receive the basis and the remaining columns are set to zero.
:::

## Characteristic and minimal polynomial {#charpoly}

::: function
void fq_mat_similarity(fq_mat_t A, slong r, fq_t d, const fq_ctx_t ctx)

Replaces the square matrix $A$ by $P^{-1} A P$, where $P$ is the
identity matrix with the entries in row $r$ and every column
$j \notin \{r - 1, r\}$ replaced by $d$. Such transforms preserve the
characteristic and minimal polynomials.
:::

The characteristic and minimal polynomials of a matrix, `fq_mat_charpoly`
and `fq_mat_minpoly`, are declared in `fq_poly.h` and documented with
[](fq_poly).
