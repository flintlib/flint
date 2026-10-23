---
title: fq_nmod_mat.h
header: fq_nmod_mat.h
template-headers: fq_mat_templates
description: Dense matrices over finite fields $\mathbb{F}_q$ in the word-size characteristic.
category: Finite fields
order: 13
---

This module implements dense matrices over a finite field $\mathbb{F}_q$ in the representation of [](fq_nmod) (a prime $p < 2^{\bits}$, elements stored as `nmod_poly`). It is generated from the same template as [](fq_mat) and has the same functions, with `fq_nmod` in place of `fq`. As there, the linear algebra functions are wrappers of the generic `gr_mat` algorithms, and multiplication uses Kronecker substitution for all but small matrices.

::: summary algorithms/linear-algebra/gaussian-elimination :::

## Types and conventions {#types}

::: type
fq_nmod_mat_struct
fq_nmod_mat_t

A dense matrix with `r` rows and `c` columns. The entries are stored in
the array `entries` row by row, with `stride` elements between the
starts of consecutive rows (equal to `c` except for windows). `fq_nmod_mat_t`
is an array of one structure.
:::

## Memory management {#memory}

::: function
void fq_nmod_mat_init(fq_nmod_mat_t mat, slong rows, slong cols, const fq_nmod_ctx_t ctx)

Initialises `mat` as a `rows` by `cols` zero matrix.

::: details
Generated from the same template as `fq_mat_init`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_init_set(fq_nmod_mat_t mat, const fq_nmod_mat_t src, const fq_nmod_ctx_t ctx)

Initialises `mat` as a copy of `src`, with the same dimensions.

::: details
Generated from the same template as `fq_mat_init_set`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_clear(fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)

Frees the memory used by `mat`.

::: details
Generated from the same template as `fq_mat_clear`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_window_init(fq_nmod_mat_t window, const fq_nmod_mat_t mat, slong r1, slong c1, slong r2, slong c2, const fq_nmod_ctx_t ctx)

Initialises `window` as the submatrix of `mat` formed by the rows
$r_1 \le i < r_2$ and the columns $c_1 \le j < c_2$. The entries are
shared with `mat`, not copied, so changes to one are visible in the
other. `mat` must not be cleared while the window is in use.

::: details
Generated from the same template as `fq_mat_window_init`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_window_clear(fq_nmod_mat_t window, const fq_nmod_ctx_t ctx)

Clears a window. The entries, which belong to the underlying matrix, are
not freed.

::: details
Generated from the same template as `fq_mat_window_clear`, whose documentation has the details.
:::
:::

## Entries and dimensions {#entries}

::: function
fq_nmod_struct * fq_nmod_mat_entry(const fq_nmod_mat_t mat, slong i, slong j)

Returns a pointer to the entry in row $i$ and column $j$, indexed from
zero. No bounds checking is performed.

::: details
Generated from the same template as `fq_mat_entry`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_entry_set(fq_nmod_mat_t mat, slong i, slong j, const fq_nmod_t x, const fq_nmod_ctx_t ctx)

Sets the entry in row $i$ and column $j$ to $x$.

::: details
Generated from the same template as `fq_mat_entry_set`, whose documentation has the details.
:::
:::

::: function
slong fq_nmod_mat_nrows(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
slong fq_nmod_mat_ncols(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)

Return the number of rows, respectively columns, of `mat`.

::: details
Generated from the same template as `fq_mat_nrows`, whose documentation has the details.
:::
:::

## Assignment and permutations {#assignment}

::: function
void fq_nmod_mat_set(fq_nmod_mat_t mat1, const fq_nmod_mat_t mat2, const fq_nmod_ctx_t ctx)

Copies `mat2` into `mat1`, which must have the same dimensions.

::: details
Generated from the same template as `fq_mat_set`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_swap(fq_nmod_mat_t mat1, fq_nmod_mat_t mat2, const fq_nmod_ctx_t ctx)

Swaps the two matrices by swapping their structures; the dimensions may
differ. Not for windows.

::: details
Generated from the same template as `fq_mat_swap`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_swap_entrywise(fq_nmod_mat_t mat1, fq_nmod_mat_t mat2, const fq_nmod_ctx_t ctx)

Swaps the entries of two matrices of the same dimensions one by one;
this also works for windows.

::: details
Generated from the same template as `fq_mat_swap_entrywise`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_zero(fq_nmod_mat_t A, const fq_nmod_ctx_t ctx)
void fq_nmod_mat_one(fq_nmod_mat_t A, const fq_nmod_ctx_t ctx)

Set `A` to the zero matrix, respectively to the matrix with ones on the
main diagonal and zeros elsewhere (also for non-square `A`).

::: details
Generated from the same template as `fq_mat_zero`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_swap_rows(fq_nmod_mat_t mat, slong * perm, slong r, slong s, const fq_nmod_ctx_t ctx)
void fq_nmod_mat_swap_cols(fq_nmod_mat_t mat, slong * perm, slong r, slong s, const fq_nmod_ctx_t ctx)

Swap rows, respectively columns, $r$ and $s$ of `mat`. If `perm` is not
`NULL`, its entries $r$ and $s$ are swapped as well, so that it records
the permutation.

::: details
Generated from the same template as `fq_mat_swap_rows`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_invert_rows(fq_nmod_mat_t mat, slong * perm, const fq_nmod_ctx_t ctx)
void fq_nmod_mat_invert_cols(fq_nmod_mat_t mat, slong * perm, const fq_nmod_ctx_t ctx)

Reverse the order of the rows, respectively columns, of `mat`, by
swapping row $i$ and row $r - 1 - i$ for $0 \le i < \lfloor r/2
\rfloor$ (and likewise for columns). If `perm` is not `NULL`, the same
permutation is applied to it.

::: details
Generated from the same template as `fq_mat_invert_rows`, whose documentation has the details.
:::
:::

## Conversion and concatenation {#conversion}

::: function
void fq_nmod_mat_set_nmod_mat(fq_nmod_mat_t mat1, const nmod_mat_t mat2, const fq_nmod_ctx_t ctx)
void fq_nmod_mat_set_fmpz_mod_mat(fq_nmod_mat_t mat1, const fmpz_mod_mat_t mat2, const fq_nmod_ctx_t ctx)

Set `mat1` to the image of the matrix `mat2` over $\mathbb{F}_p$ (or over
$\Z/n\Z$, whose entries are mapped to $\mathbb{F}_p$ by reduction) in
$\mathbb{F}_q$. The dimensions must agree.

::: details
Generated from the same template as `fq_mat_set_nmod_mat`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_concat_vertical(fq_nmod_mat_t res, const fq_nmod_mat_t mat1, const fq_nmod_mat_t mat2, const fq_nmod_ctx_t ctx)
void fq_nmod_mat_concat_horizontal(fq_nmod_mat_t res, const fq_nmod_mat_t mat1, const fq_nmod_mat_t mat2, const fq_nmod_ctx_t ctx)

Set `res` to the block matrix $\begin{pmatrix} A \\ B \end{pmatrix}$,
respectively $(A \ B)$, for $A$ = `mat1` and $B$ = `mat2`. For vertical
concatenation of an $m \times n$ and a $k \times n$ matrix, `res` must be
$(m + k) \times n$; for horizontal concatenation of an $m \times n$ and
an $m \times k$ matrix, $m \times (n + k)$.

::: details
Generated from the same template as `fq_mat_concat_vertical`, whose documentation has the details.
:::
:::

## Comparison {#comparison}

::: function
int fq_nmod_mat_equal(const fq_nmod_mat_t mat1, const fq_nmod_mat_t mat2, const fq_nmod_ctx_t ctx)

Returns whether the matrices have the same dimensions and entries.

::: details
Generated from the same template as `fq_mat_equal`, whose documentation has the details.
:::
:::

::: function
int fq_nmod_mat_is_zero(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
int fq_nmod_mat_is_one(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)

Return whether all entries are zero, respectively whether the diagonal
entries are one and the others zero. An empty matrix is both zero and
one.

::: details
Generated from the same template as `fq_mat_is_zero`, whose documentation has the details.
:::
:::

::: function
int fq_nmod_mat_is_empty(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
int fq_nmod_mat_is_square(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)

Return whether `mat` has zero rows or zero columns, respectively whether
it has as many rows as columns.

::: details
Generated from the same template as `fq_mat_is_empty`, whose documentation has the details.
:::
:::

## Random matrices {#random}

::: function
void fq_nmod_mat_randtest(fq_nmod_mat_t mat, flint_rand_t state, const fq_nmod_ctx_t ctx)

Sets the entries of `mat` to random elements for testing (see
`fq_nmod_randtest`).

::: details
Generated from the same template as `fq_mat_randtest`, whose documentation has the details.
:::
:::

::: function
int fq_nmod_mat_randpermdiag(fq_nmod_mat_t mat, flint_rand_t state, fq_nmod_struct * diag, slong n, const fq_nmod_ctx_t ctx)

Sets `mat` to a random row and column permutation of the diagonal matrix
whose first $n$ diagonal entries are `diag` and whose other entries are
zero; $n$ must not exceed the number of rows or columns. Returns the
parity of the permutation (0 for even, 1 for odd).

::: details
Generated from the same template as `fq_mat_randpermdiag`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_randrank(fq_nmod_mat_t mat, flint_rand_t state, slong rank, const fq_nmod_ctx_t ctx)

Sets `mat` to a random sparse matrix of the given rank, with exactly
`rank` nonzero entries, which are random nonzero elements.
`fq_nmod_mat_randops` turns it into a dense matrix of the same rank.

::: details
Generated from the same template as `fq_mat_randrank`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_randops(fq_nmod_mat_t mat, flint_rand_t state, slong count, const fq_nmod_ctx_t ctx)

Applies `count` random elementary row or column operations (additions or
subtractions of a row or column to another) to `mat`. This preserves the
rank and, for square matrices, the determinant.

::: details
Generated from the same template as `fq_mat_randops`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_randtril(fq_nmod_mat_t mat, flint_rand_t state, int unit, const fq_nmod_ctx_t ctx)
void fq_nmod_mat_randtriu(fq_nmod_mat_t mat, flint_rand_t state, int unit, const fq_nmod_ctx_t ctx)

Set `mat` to a random lower, respectively upper, triangular matrix,
with ones on the diagonal if `unit` is nonzero and random nonzero
diagonal entries otherwise.

::: details
Generated from the same template as `fq_mat_randtril`, whose documentation has the details.
:::
:::

## Input and output {#io}

::: function
int fq_nmod_mat_fprint(FILE * file, const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
int fq_nmod_mat_print(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)

Print `mat` to `file`, respectively `stdout`: the number of rows and
columns separated by a space, two spaces, and the entries in row-major
order in the format of `fq_nmod_fprint`, separated by spaces. Return a
positive value on success and a non-positive value on failure.

::: details
Generated from the same template as `fq_mat_fprint`, whose documentation has the details.
:::
:::

::: function
int fq_nmod_mat_fprint_pretty(FILE * file, const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
int fq_nmod_mat_print_pretty(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)

Print `mat` to `file`, respectively `stdout`, as a bracketed list of
bracketed rows, with the entries in the format of `fq_nmod_fprint_pretty`.
Return a positive value on success.

::: details
Generated from the same template as `fq_mat_fprint_pretty`, whose documentation has the details.
:::
:::

## Arithmetic {#arithmetic}

::: function
void fq_nmod_mat_transpose(fq_nmod_mat_t B, const fq_nmod_mat_t A, const fq_nmod_ctx_t ctx)

Sets $B$ to $A^T$. $B$ must have the transposed dimensions of $A$;
aliasing is allowed for square matrices.

::: details
Generated from the same template as `fq_mat_transpose`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_add(fq_nmod_mat_t C, const fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)
void fq_nmod_mat_sub(fq_nmod_mat_t C, const fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)

Set $C$ to $A + B$, respectively $A - B$. All three matrices must have
the same dimensions.

::: details
Generated from the same template as `fq_mat_add`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_neg(fq_nmod_mat_t B, const fq_nmod_mat_t A, const fq_nmod_ctx_t ctx)

Sets $B$ to $-A$.

::: details
Generated from the same template as `fq_mat_neg`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_scalar_mul(fq_nmod_mat_t B, const fq_nmod_mat_t A, const fq_nmod_t c, const fq_nmod_ctx_t ctx)
void fq_nmod_mat_scalar_mul_fmpz(fq_nmod_mat_t B, const fq_nmod_mat_t A, const fmpz_t c, const fq_nmod_ctx_t ctx)
void fq_nmod_mat_scalar_mul_si(fq_nmod_mat_t B, const fq_nmod_mat_t A, const slong c, const fq_nmod_ctx_t ctx)
void fq_nmod_mat_scalar_mul_ui(fq_nmod_mat_t B, const fq_nmod_mat_t A, const ulong c, const fq_nmod_ctx_t ctx)

Set $B$ to $cA$, where $c$ is an element of $\mathbb{F}_q$ or an integer
mapped to $\mathbb{F}_p$.

::: details
Generated from the same template as `fq_mat_scalar_mul`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_mul(fq_nmod_mat_t C, const fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)

Sets $C = AB$, where $A$ is $m \times k$, $B$ is $k \times n$ and $C$ is
$m \times n$. Aliasing is allowed. The function uses `fq_nmod_mat_mul_KS` if
`FQ_NMOD_MAT_MUL_KS_CUTOFF` returns 1 for $m$ and $n$, and
`fq_nmod_mat_mul_classical` otherwise.

::: details
Generated from the same template as `fq_mat_mul`, whose documentation has the details.
:::
:::

::: function
int FQ_NMOD_MAT_MUL_KS_CUTOFF(slong r, slong c, const fq_nmod_ctx_t ctx)

Returns 1 if `fq_nmod_mat_mul` should use Kronecker substitution for an
$r \times k$ times $k \times c$ product, which is the case if
$2(r + 1)c > \max(0, 12 - d)^2$, in particular always for $d \ge 12$.
:::

::: function
void fq_nmod_mat_mul_classical(fq_nmod_mat_t C, const fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)

Sets $C = AB$ by computing each entry as a dot product with
`_fq_nmod_vec_dot`, which reduces modulo the defining polynomial only once per
entry. Aliasing is allowed.

::: details
Generated from the same template as `fq_mat_mul_classical`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_mul_KS(fq_nmod_mat_t C, const fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)

Sets $C = AB$ by Kronecker substitution. $C$ must not be aliased with
$A$ or $B$.

::: details
Generated from the same template as `fq_mat_mul_KS`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_submul(fq_nmod_mat_t D, const fq_nmod_mat_t C, const fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)

Sets $D = C - AB$. $D$ may be aliased with $C$, $A$ or $B$.

::: details
Generated from the same template as `fq_mat_submul`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_mul_vec(fq_nmod_struct * c, const fq_nmod_mat_t A, const fq_nmod_struct * b, slong blen, const fq_nmod_ctx_t ctx)
void fq_nmod_mat_mul_vec_ptr(fq_nmod_struct * const * c, const fq_nmod_mat_t A, const fq_nmod_struct * const * b, slong blen, const fq_nmod_ctx_t ctx)

Set the vector $c$ to the matrix-vector product $Ab$. The vector
`(b, blen)` is truncated or extended by zeros to the number of columns
of $A$, and $c$ receives as many entries as $A$ has rows. The `_ptr`
version takes arrays of pointers to the entries. $c$ must not alias $b$.

::: details
Generated from the same template as `fq_mat_mul_vec`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_vec_mul(fq_nmod_struct * c, const fq_nmod_struct * a, slong alen, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)
void fq_nmod_mat_vec_mul_ptr(fq_nmod_struct * const * c, const fq_nmod_struct * const * a, slong alen, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)

Set the vector $c$ to the vector-matrix product $aB$. The vector
`(a, alen)` is truncated or extended by zeros to the number of rows of
$B$, and $c$ receives as many entries as $B$ has columns. $c$ must not
alias $a$.

::: details
Generated from the same template as `fq_mat_vec_mul`, whose documentation has the details.
:::
:::

## Gaussian elimination {#elimination}

::: function
slong fq_nmod_mat_lu(slong * P, fq_nmod_mat_t A, int rank_check, const fq_nmod_ctx_t ctx)

Computes a generalised LU decomposition $PA = LU$ of $A$ in place, and
returns the rank of $A$. `P` must have room for as many entries as $A$
has rows; it receives the row permutation.

::: details
Generated from the same template as `fq_mat_lu`, whose documentation has the details.
:::
:::

::: function
slong fq_nmod_mat_rref(fq_nmod_mat_t B, const fq_nmod_mat_t A, const fq_nmod_ctx_t ctx)

Sets $B$ to the reduced row echelon form of $A$ and returns the rank of
$A$. The form is computed from an LU decomposition. $B$ may be aliased
with $A$. If $A$ is zero, $B$
is set to zero and 0 is returned.

::: details
Generated from the same template as `fq_mat_rref`, whose documentation has the details.
:::
:::

::: function
slong fq_nmod_mat_reduce_row(fq_nmod_mat_t A, slong * P, slong * L, slong m, const fq_nmod_ctx_t ctx)

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

::: details
Generated from the same template as `fq_mat_reduce_row`, whose documentation has the details.
:::
:::

::: function
slong fq_nmod_mat_rank(const fq_nmod_mat_t A, const fq_nmod_ctx_t ctx)

Returns the rank of $A$, computed by LU decomposition.

::: details
Generated from the same template as `fq_mat_rank`, whose documentation has the details.
:::
:::

::: function
int fq_nmod_mat_inv(fq_nmod_mat_t B, fq_nmod_mat_t A, const fq_nmod_ctx_t ctx)

If the square matrix $A$ is invertible, sets $B = A^{-1}$ and returns 1;
otherwise returns 0 and the entries of $B$ are undefined. $A$ is not
modified (although not declared `const`), and $B$ may be aliased with
$A$.

::: details
Generated from the same template as `fq_mat_inv`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_solve_tril(fq_nmod_mat_t X, const fq_nmod_mat_t L, const fq_nmod_mat_t B, int unit, const fq_nmod_ctx_t ctx)
void fq_nmod_mat_solve_triu(fq_nmod_mat_t X, const fq_nmod_mat_t U, const fq_nmod_mat_t B, int unit, const fq_nmod_ctx_t ctx)

Set $X = L^{-1} B$, respectively $X = U^{-1} B$, for a square lower,
respectively upper, triangular matrix with nonzero diagonal entries. If
`unit` is nonzero, the diagonal is assumed to consist of ones and is not
read. Only the triangular part of the matrix is read. Aliasing of $X$ and
$B$ is allowed. Classical substitution is used if $B$ has fewer than 64
rows or columns, and block recursion otherwise.

::: details
Generated from the same template as `fq_mat_solve_tril`, whose documentation has the details.
:::
:::

::: function
int fq_nmod_mat_solve(fq_nmod_mat_t X, const fq_nmod_mat_t A, const fq_nmod_mat_t C, const fq_nmod_ctx_t ctx)

Solves $AX = C$ for a square matrix $A$. Returns 1 and sets $X$ if $A$ is
nonsingular; otherwise returns 0 and $X$ is undefined.

::: details
Generated from the same template as `fq_mat_solve`, whose documentation has the details.
:::
:::

::: function
int fq_nmod_mat_can_solve(fq_nmod_mat_t X, const fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)

Solves $AX = B$ for a matrix $A$ of any shape: if a solution exists, sets
$X$ to one and returns 1, otherwise returns 0 and $X$ is undefined. For
an $m \times n$ matrix $A$ and an $m \times k$ matrix $B$, $X$ must be
$n \times k$.

::: details
Generated from the same template as `fq_mat_can_solve`, whose documentation has the details.
:::
:::

::: function
slong fq_nmod_mat_nullspace(fq_nmod_mat_t X, const fq_nmod_mat_t A, const fq_nmod_ctx_t ctx)

Computes a basis of the right nullspace $\{v : Av = 0\}$ of $A$ and
returns its dimension $k$ (the nullity). $X$ must have as many rows and
at least as many columns as $A$ has columns; its first $k$ columns
receive the basis and the remaining columns are set to zero.

::: details
Generated from the same template as `fq_mat_nullspace`, whose documentation has the details.
:::
:::

## Characteristic and minimal polynomial {#charpoly}

::: function
void fq_nmod_mat_similarity(fq_nmod_mat_t A, slong r, fq_nmod_t d, const fq_nmod_ctx_t ctx)

Replaces the square matrix $A$ by $P^{-1} A P$, where $P$ is the
identity matrix with the entries in row $r$ and every column
$j \notin \{r - 1, r\}$ replaced by $d$. Such transforms preserve the
characteristic and minimal polynomials.

::: details
Generated from the same template as `fq_mat_similarity`, whose documentation has the details.
:::
:::
