---
title: fmpq_mat.h
header: fmpq_mat.h
description: Dense matrices over the rational numbers, with conversions to integer matrices, multiplication, determinants, linear system solving, echelon forms and characteristic polynomials computed over the integers.
category: Matrices
order: 3
---

This module implements dense matrices with entries in $\Q$, stored as
arrays of `fmpq` in canonical form (reduced fractions with positive
denominators). Entries can be read and modified individually, but all
nontrivial operations are carried out by clearing denominators, computing
with integer matrices in [](fmpz_mat), and converting back: products,
determinants and characteristic polynomials use common denominators of
rows, columns or the whole matrix, and the solvers pass integer systems to
the multimodular and $p$-adic solvers of `fmpz_mat`. The conversion
functions are public, so that users can organise such computations
themselves.

::: summary algorithms/linear-algebra/multimodular :::

## Types and conventions {#types}

::: type
fmpq_mat_struct
fmpq_mat_t

An `fmpq_mat_struct` has the fields `entries` (a pointer to the entry in
row 0, column 0), `r` and `c` (the numbers of rows and columns) and
`stride` (the distance between the starts of consecutive rows). An
`fmpq_mat_t` is an array of one `fmpq_mat_struct`.
:::

- *Canonical form.* Functions assume that the input entries are in
  canonical form and produce canonical output. Entries modified directly
  through `fmpq_mat_entry_num` and `fmpq_mat_entry_den` must be put back
  into canonical form with `fmpq_canonicalise`.
- *Dimensions.* Outputs must be initialised with the dimensions of the
  result; dimensions are checked only by some functions.
- *Aliasing.* Entrywise operations allow aliasing; for the others, each
  description states whether the output may be an input.

## Memory management {#memory}

::: function
void fmpq_mat_init(fmpq_mat_t mat, slong rows, slong cols)

Initialises `mat` as a `rows` by `cols` matrix of zeros. Throws if
$\mathit{rows} \cdot \mathit{cols}$ overflows an `slong`.
:::

::: function
void fmpq_mat_init_set(fmpq_mat_t mat1, const fmpq_mat_t mat2)

Initialises `mat1` as a copy of `mat2`.
:::

::: function
void fmpq_mat_clear(fmpq_mat_t mat)

Frees the memory used by `mat`. Must not be called on a window.
:::

::: function
void fmpq_mat_swap(fmpq_mat_t mat1, fmpq_mat_t mat2)
void fmpq_mat_swap_entrywise(fmpq_mat_t mat1, fmpq_mat_t mat2)

Swaps two matrices in constant time by swapping their structs (not for
windows), respectively swaps the entries of two matrices of the same
dimensions one by one.
:::

## Entries and dimensions {#entries}

::: function
fmpq * fmpq_mat_entry(const fmpq_mat_t mat, slong i, slong j)

Returns a pointer to the entry in row $i$ and column $j$.
:::

::: function
fmpz * fmpq_mat_entry_num(const fmpq_mat_t mat, slong i, slong j)
fmpz * fmpq_mat_entry_den(const fmpq_mat_t mat, slong i, slong j)

Returns a pointer to the numerator, respectively denominator, of the entry
in row $i$ and column $j$.
:::

::: function
slong fmpq_mat_nrows(const fmpq_mat_t mat)
slong fmpq_mat_ncols(const fmpq_mat_t mat)

Returns the number of rows, respectively columns, of `mat`.
:::

::: function
int fmpq_mat_is_empty(const fmpq_mat_t mat)
int fmpq_mat_is_square(const fmpq_mat_t mat)

Returns 1 if `mat` has no rows or no columns, respectively as many rows as
columns, and 0 otherwise.
:::

## Windows and concatenation {#windows}

::: function
void fmpq_mat_window_init(fmpq_mat_t window, const fmpq_mat_t mat, slong r1, slong c1, slong r2, slong c2)

Initialises `window` as the submatrix of `mat` formed by rows $r_1$ to
$r_2 - 1$ and columns $c_1$ to $c_2 - 1$, sharing its entries with `mat`.
The bounds are checked only by `FLINT_ASSERT`.
:::

::: function
void fmpq_mat_window_clear(fmpq_mat_t window)

Releases a window. This does nothing at present.
:::

::: function
void fmpq_mat_concat_horizontal(fmpq_mat_t res, const fmpq_mat_t mat1, const fmpq_mat_t mat2)
void fmpq_mat_concat_vertical(fmpq_mat_t res, const fmpq_mat_t mat1, const fmpq_mat_t mat2)

Sets `res` to `mat1` placed to the left of, respectively above, `mat2`,
whose dimensions must be compatible. Wrappers of
`gr_mat_concat_horizontal` and `gr_mat_concat_vertical`.
:::

## Input and output {#io}

::: function
void fmpq_mat_print(const fmpq_mat_t mat)

Prints `mat` to `stdout`: a header `<r x c matrix over Q>`, then each row
in square brackets with entries separated by commas, then an empty line.
:::

## Random and special matrices {#random}

::: function
void fmpq_mat_randbits(fmpq_mat_t mat, flint_rand_t state, flint_bitcnt_t bits)
void fmpq_mat_randtest(fmpq_mat_t mat, flint_rand_t state, flint_bitcnt_t bits)

Sets the entries to random fractions with numerator and denominator of
`bits` bits (`fmpq_randbits`), respectively of at most `bits` bits with a
bias towards special values (`fmpq_randtest`).
:::

::: function
void fmpq_mat_hilbert_matrix(fmpq_mat_t mat)

Sets `mat` to the Hilbert matrix, with entry $(i, j)$ equal to
$1/(i + j + 1)$. `mat` need not be square.
:::

## Comparison {#comparison}

::: function
int fmpq_mat_equal(const fmpq_mat_t mat1, const fmpq_mat_t mat2)

Returns 1 if the matrices have the same dimensions and entries, and 0
otherwise.
:::

::: function
int fmpq_mat_is_zero(const fmpq_mat_t mat)
int fmpq_mat_is_one(const fmpq_mat_t mat)

Returns 1 if `mat` is zero, respectively has ones on the main diagonal
and zeros elsewhere ($\mathit{mat}$ need not be square), and 0 otherwise.
Empty matrices are both zero and one.
:::

::: function
int fmpq_mat_is_integral(const fmpq_mat_t mat)

Returns 1 if all entries have denominator 1, and 0 otherwise.
:::

## Assignment and permutations {#assignment}

::: function
void fmpq_mat_set(fmpq_mat_t dest, const fmpq_mat_t src)

Copies `src` to `dest`, which must have the same dimensions.
:::

::: function
void fmpq_mat_zero(fmpq_mat_t mat)
void fmpq_mat_one(fmpq_mat_t mat)

Sets `mat` to the zero matrix, respectively to the matrix with ones on the
main diagonal and zeros elsewhere.
:::

::: function
void fmpq_mat_transpose(fmpq_mat_t rop, const fmpq_mat_t op)

Sets `rop` to the transpose of `op`, which may be `rop` if it is square.
Wrapper of `gr_mat_transpose`.
:::

::: function
void fmpq_mat_swap_rows(fmpq_mat_t mat, slong * perm, slong r, slong s)
void fmpq_mat_swap_cols(fmpq_mat_t mat, slong * perm, slong r, slong s)

Swaps rows, respectively columns, $r$ and $s$ of `mat`, and entries $r$
and $s$ of `perm` if it is not `NULL`.
:::

::: function
void fmpq_mat_invert_rows(fmpq_mat_t mat, slong * perm)
void fmpq_mat_invert_cols(fmpq_mat_t mat, slong * perm)

Reverses the order of the rows, respectively columns, of `mat`, and of the
entries of `perm` if it is not `NULL`.
:::

## Arithmetic {#arithmetic}

::: function
void fmpq_mat_add(fmpq_mat_t mat, const fmpq_mat_t mat1, const fmpq_mat_t mat2)
void fmpq_mat_sub(fmpq_mat_t mat, const fmpq_mat_t mat1, const fmpq_mat_t mat2)
void fmpq_mat_neg(fmpq_mat_t rop, const fmpq_mat_t op)

Sets `mat` to `mat1` + `mat2`, respectively `mat1` $-$ `mat2`, and `rop`
to $-$`op`. Aliasing is allowed.
:::

::: function
void fmpq_mat_scalar_mul_fmpq(fmpq_mat_t rop, const fmpq_mat_t op, const fmpq_t x)
void fmpq_mat_scalar_mul_fmpz(fmpq_mat_t rop, const fmpq_mat_t op, const fmpz_t x)
void fmpq_mat_scalar_mul_si(fmpq_mat_t rop, const fmpq_mat_t op, const slong x)
void fmpq_mat_scalar_mul_ui(fmpq_mat_t rop, const fmpq_mat_t op, const ulong x)

Sets `rop` to `op` multiplied by $x$. Aliasing is allowed. The `si` and
`ui` versions wrap `gr_mat_mul_si` and `gr_mat_mul_ui`.
:::

::: function
void fmpq_mat_scalar_div_fmpz(fmpq_mat_t rop, const fmpq_mat_t op, const fmpz_t x)
void fmpq_mat_scalar_div_si(fmpq_mat_t rop, const fmpq_mat_t op, const slong x)
void fmpq_mat_scalar_div_ui(fmpq_mat_t rop, const fmpq_mat_t op, const ulong x)

Sets `rop` to `op` divided by $x \ne 0$. Aliasing is allowed.
:::

## Conversion to and from integer matrices {#conversion}

::: function
int fmpq_mat_get_fmpz_mat(fmpz_mat_t dest, const fmpq_mat_t mat)

Sets `dest` to `mat` and returns 1 if all entries are integers; otherwise
returns 0, and `dest` may have been partially written.
:::

::: function
void fmpq_mat_set_fmpz_mat(fmpq_mat_t dest, const fmpz_mat_t src)

Sets `dest` to the integer matrix `src`.
:::

::: function
void fmpq_mat_get_fmpz_mat_entrywise(fmpz_mat_t num, fmpz_mat_t den, const fmpq_mat_t mat)

Sets `num` and `den` to the matrices of numerators and denominators of the
entries of `mat`.
:::

::: function
void fmpq_mat_get_fmpz_mat_matwise(fmpz_mat_t num, fmpz_t den, const fmpq_mat_t mat)

Sets `den` to the least common denominator $d$ of all entries and `num` to
the integer matrix $d \cdot \mathit{mat}$. For an empty matrix, $d = 1$.
:::

::: function
void fmpq_mat_get_fmpz_mat_rowwise(fmpz_mat_t num, fmpz * den, const fmpq_mat_t mat)
void fmpq_mat_get_fmpz_mat_colwise(fmpz_mat_t num, fmpz * den, const fmpq_mat_t mat)

Clears denominators row by row, respectively column by column: `den[i]`
is set to the least common denominator of row (column) $i$, and row
(column) $i$ of `num` to the corresponding integer vector. `den` may be
`NULL` if the denominators are not wanted.
:::

::: function
void fmpq_mat_get_fmpz_mat_rowwise_2(fmpz_mat_t num, fmpz_mat_t num2, fmpz * den, const fmpq_mat_t mat, const fmpq_mat_t mat2)

Clears denominators row by row in the two matrices `mat` and `mat2` with
the same number of rows jointly, as if they were concatenated
horizontally: `den[i]` is the least common denominator of row $i$ of both,
and `num`, `num2` receive the scaled rows. This turns a system
$AX = B$ over $\Q$ into an equivalent integer system. `den` may be
`NULL`.
:::

::: function
void fmpq_mat_set_fmpz_mat_div_fmpz(fmpq_mat_t X, const fmpz_mat_t Xmod, const fmpz_t div)

Sets $X$ to the matrix `Xmod` divided by the nonzero integer `div`,
reducing every entry to canonical form. This converts a pair (`X`, `den`)
returned by the `fmpz_mat` solvers to a rational matrix.
:::

## Modular reduction and rational reconstruction {#modular}

::: function
void fmpq_mat_get_fmpz_mat_mod_fmpz(fmpz_mat_t dest, const fmpq_mat_t mat, const fmpz_t mod)

Sets each entry of `dest` to the corresponding entry of `mat` modulo
`mod`, in $[0, \mathit{mod})$ (see `fmpq_mod_fmpz`). All denominators must
be invertible modulo `mod`; otherwise the corresponding entries are
undefined.
:::

::: function
int fmpq_mat_set_fmpz_mat_mod_fmpz(fmpq_mat_t X, const fmpz_mat_t Xmod, const fmpz_t mod)

Sets $X$ to the rational matrix congruent to `Xmod` modulo `mod` whose
entries are obtained by rational reconstruction with balanced bounds
(`_fmpq_reconstruct_fmpz`,
[](algorithms/number-theory/rational-reconstruction)), using a running
common denominator: each
entry is multiplied by the product of the denominators found so far
before it is reconstructed. Returns 1 on success and 0 if some entry
cannot be reconstructed, in which case $X$ is partially written. The
result is correct if every numerator and every partial common denominator
is less than $\sqrt{\mathit{mod}/2}$.
:::

## Multiplication {#mul}

::: function
void fmpq_mat_mul(fmpq_mat_t C, const fmpq_mat_t A, const fmpq_mat_t B)
void fmpq_mat_mul_cleared(fmpq_mat_t C, const fmpq_mat_t A, const fmpq_mat_t B)
void fmpq_mat_mul_direct(fmpq_mat_t C, const fmpq_mat_t A, const fmpq_mat_t B)

Sets $C = AB$. `fmpq_mat_mul` calls `fmpq_mat_mul_cleared`, which clears
denominators of the rows of $A$ and of the columns of $B$, multiplies the
integer matrices with `fmpz_mat_mul` and divides entry $(i, j)$ by the
product of the $i$-th row and $j$-th column denominators.
`fmpq_mat_mul_direct` computes the dot products in rational arithmetic and
throws if $C$ is $A$ or $B$. `fmpq_mat_mul` and `fmpq_mat_mul_cleared`
allow aliasing.
:::

::: function
void fmpq_mat_mul_fmpz_mat(fmpq_mat_t C, const fmpq_mat_t A, const fmpz_mat_t B)
void fmpq_mat_mul_r_fmpz_mat(fmpq_mat_t C, const fmpz_mat_t A, const fmpq_mat_t B)

Sets $C = AB$ where $B$, respectively $A$, is an integer matrix, by
clearing the denominators of the rows of $A$, respectively the columns of
$B$, and one `fmpz_mat_mul`.
:::

::: function
void fmpq_mat_mul_fmpz_vec(fmpq * c, const fmpq_mat_t A, const fmpz * b, slong blen)
void fmpq_mat_mul_fmpq_vec(fmpq * c, const fmpq_mat_t A, const fmpq * b, slong blen)
void fmpq_mat_mul_fmpz_vec_ptr(fmpq * const * c, const fmpq_mat_t A, const fmpz * const * b, slong blen)
void fmpq_mat_mul_fmpq_vec_ptr(fmpq * const * c, const fmpq_mat_t A, const fmpq * const * b, slong blen)

Sets the column vector $c$ of length `A->r` to $Ab$, where $b$ is
truncated or zero-extended to the number of columns of $A$. The `ptr`
versions take arrays of pointers. The `fmpq_vec` versions clear the
denominator of $b$ first. $c$ must not overlap $b$.
:::

::: function
void fmpq_mat_fmpz_vec_mul(fmpq * c, const fmpz * a, slong alen, const fmpq_mat_t B)
void fmpq_mat_fmpq_vec_mul(fmpq * c, const fmpq * a, slong alen, const fmpq_mat_t B)
void fmpq_mat_fmpz_vec_mul_ptr(fmpq * const * c, const fmpz * const * a, slong alen, const fmpq_mat_t B)
void fmpq_mat_fmpq_vec_mul_ptr(fmpq * const * c, const fmpq * const * a, slong alen, const fmpq_mat_t B)

Sets the row vector $c$ of length `B->c` to $aB$, where $a$ is truncated
or zero-extended to the number of rows of $B$. The `ptr` versions take
arrays of pointers. $c$ must not overlap $a$.
:::

::: function
void fmpq_mat_kronecker_product(fmpq_mat_t C, const fmpq_mat_t A, const fmpq_mat_t B)

Sets $C$ to the Kronecker product $A \otimes B$; for $A$ of size
$m \times n$ and $B$ of size $k \times l$, $C$ must be $mk \times nl$.
$C$ must not be aliased with $A$ or $B$.
:::

## Trace, determinant and permanent {#det}

::: function
void fmpq_mat_trace(fmpq_t trace, const fmpq_mat_t mat)

Sets `trace` to the sum of the diagonal entries of the square matrix
`mat`. Throws if `mat` is not square.
:::

::: function
void fmpq_mat_det(fmpq_t det, const fmpq_mat_t mat)

Sets `det` to the determinant of the square matrix `mat`. For dimensions
up to 2 an explicit formula is used; otherwise the denominators of the
rows are cleared, giving $\operatorname{diag}(d_i) \cdot \mathit{mat}$,
and the determinant of that integer matrix (`fmpz_mat_det`) is divided by
$\prod_i d_i$.
:::

::: function
int fmpq_mat_permanent(fmpq_t res, const fmpq_mat_t A)

Sets `res` to the permanent of the square matrix $A$, computed like the
determinant with `fmpz_mat_permanent` after clearing denominators row by
row, and returns 1.
:::

## Solving linear systems {#solve}

The functions with `fmpz_mat` in their name take integer matrices $A$
and $B$ and return the rational solution $X$; the others take rational
matrices and clear denominators row by row, which turns $AX = B$ into an
equivalent integer system ([](#conversion)). See [](fmpq_mat/solve).

::: function
int fmpq_mat_solve(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)
int fmpq_mat_solve_fmpz_mat(fmpq_mat_t X, const fmpz_mat_t A, const fmpz_mat_t B)
int fmpq_mat_solve_fraction_free(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)
int fmpq_mat_solve_fmpz_mat_fraction_free(fmpq_mat_t X, const fmpz_mat_t A, const fmpz_mat_t B)
int fmpq_mat_solve_dixon(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)
int fmpq_mat_solve_fmpz_mat_dixon(fmpq_mat_t X, const fmpz_mat_t A, const fmpz_mat_t B)
int fmpq_mat_solve_multi_mod(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)
int fmpq_mat_solve_fmpz_mat_multi_mod(fmpq_mat_t X, const fmpz_mat_t A, const fmpz_mat_t B)

For a square $A$, sets $X$ to the solution of $AX = B$ and returns 1 if
$A$ is nonsingular, and returns 0 if $A$ is singular. The functions
without algorithm suffix choose; the others use `fmpz_mat_solve` (which
uses fraction-free elimination up to dimension 15),
`fmpz_mat_solve_dixon_den`, respectively `fmpz_mat_solve_multi_mod_den`.

::: details
`fmpq_mat_solve` and `fmpq_mat_solve_fmpz_mat` use the fraction-free
versions for $n \le 15$, the multimodular versions if $B$ has at least
$n/2$ columns, and the $p$-adic versions otherwise. The functions with
algorithm suffix (except the fraction-free ones) throw if the dimensions
are incompatible. On a return value of 0, $X$ is unchanged.
:::
:::

::: function
int fmpq_mat_can_solve(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)
int fmpq_mat_can_solve_fraction_free(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)
int fmpq_mat_can_solve_dixon(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)
int fmpq_mat_can_solve_multi_mod(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)
int fmpq_mat_can_solve_fmpz_mat_dixon(fmpq_mat_t X, const fmpz_mat_t A, const fmpz_mat_t B)
int fmpq_mat_can_solve_fmpz_mat_multi_mod(fmpq_mat_t X, const fmpz_mat_t A, const fmpz_mat_t B)

For an $m \times n$ matrix $A$ of any rank and an $m \times k$ matrix
$B$, returns 1 and sets the $n \times k$ matrix $X$ to a solution of
$AX = B$ if one exists, and returns 0 if the system is inconsistent. The
answer is proven. The solution has zeros in the rows corresponding to the
non-pivot columns of $A$.

::: details
These call `fmpz_mat_can_solve_fflu`, `fmpz_mat_can_solve_dixon_den` and
`fmpz_mat_can_solve_multi_mod_den` on the integer systems.
`fmpq_mat_can_solve` uses the fraction-free version for $m \le 15$ and
otherwise the modular method of `fmpz_mat_can_solve`. The dimensions are
checked and an exception is thrown if they are incompatible.
:::
:::

::: function
int fmpq_mat_inv(fmpq_mat_t B, const fmpq_mat_t A)

Sets $B$ to the inverse of the square matrix $A$ and returns 1 if $A$ is
nonsingular; returns 0 if $A$ is singular. $B$ may be $A$. For dimension
at least 3, the denominators of the rows of $A$ are cleared,
$\operatorname{diag}(d) A = A'$, and $A' X = \operatorname{diag}(d)$ is
solved with `fmpq_mat_solve_fmpz_mat`.
:::

## Echelon form {#rref}

::: function
int fmpq_mat_pivot(slong * perm, fmpq_mat_t mat, slong r, slong c)

Helper for Gaussian elimination. Returns 1 if the entry $(r, c)$ is
nonzero. Otherwise, if some entry $(j, c)$ with $j > r$ is nonzero, swaps
the first such row with row $r$ (and the entries $r$, $j$ of `perm` if it
is not `NULL`) and returns $-1$; returns 0 if all these entries are zero.
:::

::: function
slong fmpq_mat_rref(fmpq_mat_t B, const fmpq_mat_t A)
slong fmpq_mat_rref_classical(fmpq_mat_t B, const fmpq_mat_t A)
slong fmpq_mat_rref_fraction_free(fmpq_mat_t B, const fmpq_mat_t A)

Sets $B$ to the reduced row echelon form of $A$ and returns the rank. $B$
may be $A$. `fmpq_mat_rref_classical` uses Gauss–Jordan elimination with
rational arithmetic; `fmpq_mat_rref_fraction_free` clears the
denominators of the rows and calls `fmpz_mat_rref`. `fmpq_mat_rref` uses
the classical algorithm if $A$ has at most 2 rows or columns and the
fraction-free one otherwise.
:::

## Gram–Schmidt orthogonalisation {#gso}

::: function
void fmpq_mat_gso(fmpq_mat_t B, const fmpq_mat_t A)

Sets $B$ to the Gram–Schmidt orthogonalisation of the columns of $A$:
column $i$ of $B$ is column $i$ of $A$ minus its projections onto the
previous columns of $B$, without normalisation. Columns that become zero
are skipped in later projections. $B$ may be $A$. Throws if the
dimensions differ.
:::

## Characteristic and minimal polynomials {#charpoly}

::: function
void fmpq_mat_similarity(fmpq_mat_t A, slong r, fmpq_t d)

Replaces the square matrix $A$ by $P^{-1} A P$, where $P$ is the identity
matrix with the entries in row $r$ and every column $j \notin
\{r - 1, r\}$ replaced by $d$ (see `fmpz_mat_similarity`). Wrapper of
`gr_mat_apply_row_similarity`; throws if $A$ is not square or $r$ is out
of range.
:::

::: function
void fmpq_mat_charpoly(fmpq_poly_t pol, const fmpq_mat_t mat)
void _fmpq_mat_charpoly(fmpz * coeffs, fmpz_t den, const fmpq_mat_t mat)

Sets `pol` to the characteristic polynomial of the square matrix `mat`;
the underscore version writes the $n + 1$ coefficients of the numerator
to `coeffs` and the denominator to `den`, not necessarily in canonical
form. Throws if `mat` is not square.

::: details
With $q$ the least common denominator of the entries, the characteristic
polynomial $\sum_k c_k x^k$ of the integer matrix $qA$ is computed with
`_fmpz_mat_charpoly`, and $\chi_A(x) = q^{-n} \chi_{qA}(qx)$, so the
coefficient of $x^k$ is $c_k / q^{n - k}$.
:::
:::

::: function
void fmpq_mat_minpoly(fmpq_poly_t pol, const fmpq_mat_t mat)
slong _fmpq_mat_minpoly(fmpz * coeffs, fmpz_t den, const fmpq_mat_t mat)

Sets `pol` to the minimal polynomial of the square matrix `mat`; the
underscore version writes the coefficients of the numerator to `coeffs`
(space for $n + 1$ of them) and the denominator to `den`, and returns the
length. The minimal polynomial is computed like the characteristic
polynomial, from that of the integer matrix $qA$ with
`_fmpz_mat_minpoly`. Throws if `mat` is not square.
:::
