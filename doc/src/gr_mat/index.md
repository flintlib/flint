---
title: gr_mat.h
header: gr_mat.h
description: Dense matrices over any generic ring, with the generic algorithms for multiplication, LU and fraction-free LU decomposition, solving, determinants, ranks, echelon forms, characteristic polynomials, eigenvalues and matrix functions that FLINT's matrix modules build on.
category: Generic rings
order: 6
---

This module implements dense matrices over a ring $R$ given by a
generic context `gr_ctx_t` (see [](gr)). A `gr_mat_t` stores its entries
row by row in the representation of the elements of $R$, so the same
functions serve integer, rational and modular matrices, matrices over
finite fields and number fields, matrices of real and complex balls or
floating-point numbers, and matrices over polynomial rings. The module
contains the generic algorithms for matrix multiplication (classical,
Strassen, Waksman, Rosowski), Gaussian elimination in its ordinary and
fraction-free forms, triangular and general linear solving,
determinants and permanents, ranks, echelon forms and nullspaces,
characteristic and minimal polynomials, eigenvalues and
diagonalisation, Jordan forms and functions of matrices, and QR
decompositions.

Multiplication, LU decomposition, determinants, characteristic
polynomials, triangular and nonsingular solving, inversion and pivot
selection are methods of the
coefficient context, which a ring can replace by its own
implementation; every algorithm here that multiplies matrices calls
`gr_mat_mul`, so a fast ring-specific product speeds up block LU
decomposition, triangular solving and everything built on them. The
overview [](algorithms/linear-algebra/overview) says which ring uses
which algorithm and how the specialised modules (`nmod_mat`,
`fmpz_mat`, `fmpz_mod_mat`, `arb_mat` and others) relate to the generic
code.

The guide [](gr_mat/exactness) explains how the functions behave over
rings that are not fields and over inexact rings, with complete
programs.

::: summary algorithms/linear-algebra/overview :::

## Types and conventions {#conventions}

::: type
gr_mat_struct
gr_mat_t

A matrix over the ring given by a `gr_ctx_t`. The structure has the
fields `entries` (a pointer to the first entry), `r` and `c` (the
numbers of rows and columns) and `stride` (the distance in elements
between the starts of consecutive rows). The entry in row $i$ and
column $j$ is at index $i \cdot \mathit{stride} + j$ of `entries`. A
matrix created by `gr_mat_init` has $\mathit{stride} = c$; a window
shares the entries and the stride of its parent. A `gr_mat_t` is an
array of one `gr_mat_struct`, so it is passed by reference.
:::

The following conventions hold throughout the module.

- **The context is the coefficient ring.** The `ctx` argument of every
  function describes $R$. A context for a ring of matrices
  (`gr_ctx_init_matrix_ring`, `gr_ctx_init_matrix_space`) is needed
  only to use matrices as elements of another generic structure or with
  the generic element functions such as `gr_mul`.
- **Shapes.** Outputs must be initialised with the shape of the result;
  most functions return `GR_DOMAIN` when the shapes of their arguments
  are incompatible, but not all of them check. Functions that resize
  their output say so (for example `gr_mat_nullspace`).
- **Status codes.** Functions returning `int` return `GR_SUCCESS` or a
  combination of `GR_DOMAIN` (no result in $R$: a singular matrix where
  an invertible one is needed, a non-invertible pivot, an impossible
  shape) and `GR_UNABLE` (the result may exist but cannot be computed,
  typically because a zero test is undecidable or the ring lacks a
  required operation). On failure the outputs are undefined but valid.
  Predicates return a `truth_t`.
- **Rings.** Most functions assume that $R$ is commutative. Elimination
  over a field uses inverses of pivots; fraction-free elimination needs
  exact division and works over integral domains. Algorithms designed
  for integral domains often work over other rings as long as every
  inversion they attempt succeeds, and report a failed inversion with
  `GR_DOMAIN` or `GR_UNABLE`. Rings with undecidable equality (balls,
  exact real numbers) can make a pivot search fail with `GR_UNABLE`;
  see [](gr_mat/exactness).
- **Aliasing.** The output may be the same matrix as an input for the
  entrywise operations and, unless stated otherwise, for the functions
  of this module, which use temporaries where needed. Several outputs
  of one function must be distinct.

The layout of `gr_mat_t` is that of `fmpz_mat_t`, `fmpq_mat_t`,
`fmpz_mod_mat_t`, `fq_mat_t`, `fq_nmod_mat_t`, `fq_zech_mat_t`,
`arb_mat_t`, `acb_mat_t` and `ca_mat_t`, so such matrices can be cast
to a `gr_mat_t` over the matching context. `nmod_mat_t` stores its
modulus in the matrix and must be converted.

## Memory management {#memory}

::: function
void gr_mat_init(gr_mat_t mat, slong rows, slong cols, gr_ctx_t ctx)

Initialises `mat` as a zero matrix with the given numbers of rows and
columns, which may be zero.
:::

::: function
int gr_mat_init_set(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)

Initialises `res` as a copy of `mat`.
:::

::: function
void gr_mat_clear(gr_mat_t mat, gr_ctx_t ctx)

Clears the entries and frees the memory of `mat`. A window must be
released with `gr_mat_window_clear` instead.
:::

::: function
void gr_mat_swap(gr_mat_t mat1, gr_mat_t mat2, gr_ctx_t ctx)

Swaps `mat1` and `mat2` in constant time by exchanging the structures.
:::

::: function
int gr_mat_swap_entrywise(gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)

Swaps the entries of two matrices of the same shape one by one, leaving
the structures in place. This is the way to exchange the contents of
windows.
:::

## Entries and dimensions {#entries}

::: macro
GR_MAT_ENTRY(mat, i, j, sz)

A pointer to the entry in row `i` and column `j` of `mat`, whose
elements have size `sz` bytes (`ctx->sizeof_elem`).
:::

::: function
gr_ptr gr_mat_entry_ptr(gr_mat_t mat, slong i, slong j, gr_ctx_t ctx)
gr_srcptr gr_mat_entry_srcptr(const gr_mat_t mat, slong i, slong j, gr_ctx_t ctx)

Return a pointer to the entry in row `i` and column `j`. The indices
must be in bounds.
:::

::: macro
gr_mat_nrows(mat, ctx)
gr_mat_ncols(mat, ctx)

The number of rows and of columns of `mat`.
:::

::: function
truth_t gr_mat_is_empty(const gr_mat_t mat, gr_ctx_t ctx)
truth_t gr_mat_is_square(const gr_mat_t mat, gr_ctx_t ctx)

Return whether `mat` has no entries (zero rows or zero columns),
respectively as many rows as columns. These predicates depend only on
the shape and always return `T_TRUE` or `T_FALSE`.
:::

## Windows and shallow transposes {#windows}

::: function
void gr_mat_window_init(gr_mat_t window, const gr_mat_t mat, slong r1, slong c1, slong r2, slong c2, gr_ctx_t ctx)
void gr_mat_window_clear(gr_mat_t window, gr_ctx_t ctx)

`gr_mat_window_init` makes `window` the submatrix of `mat` formed by
rows $r_1 \le i < r_2$ and columns $c_1 \le j < c_2$, sharing the
entries of `mat` (writing to the window modifies `mat`). The indices
must be within bounds. `gr_mat_window_clear` releases a window without
touching the entries.
:::

::: macro
GR_MAT_TMP_INIT_SHALLOW_TRANSPOSE(AT, A, ctx)
GR_MAT_SHALLOW_TRANSPOSE(AT, A, ctx)
GR_MAT_TMP_CLEAR_SHALLOW_TRANSPOSE(AT, ctx)

`GR_MAT_TMP_INIT_SHALLOW_TRANSPOSE` makes `AT` a temporary transpose of
`A` whose entries are shallow copies (sharing heap data) of those of
`A`, allocated on the stack when small. It can be read like a window,
and written to provided that the changes are copied back with
`GR_MAT_SHALLOW_TRANSPOSE(A, AT, ctx)`, which copies the transpose of
its second argument shallowly into its first (the two must not be
aliased). `GR_MAT_TMP_CLEAR_SHALLOW_TRANSPOSE` frees the temporary
without clearing the shared entries. The classical multiplication uses
this to turn the columns of the right factor into contiguous vectors
for dot products.
:::

## Input and output {#io}

::: function
int gr_mat_write(gr_stream_t out, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_print(const gr_mat_t mat, gr_ctx_t ctx)

Write `mat` to a stream or to standard output, as a bracketed list of
rows with one row per line.
:::

::: function
int gr_mat_set_str(gr_mat_t mat, const char * s, int resize, gr_ctx_t ctx)

Sets `mat` to the matrix described by `s`, of the form `[[e11, ...,
e1n], ..., [em1, ..., emn]]`, where each row is parsed as by
`gr_vec_set_str` and each entry by `gr_set_str`. Whitespace around
brackets and separators, including the line breaks written by
`gr_mat_write`, is ignored.

::: details
If `resize` is 1, `mat` is resized to the parsed shape; if it is 0,
`GR_DOMAIN` is returned when the parsed shape differs from that of
`mat`. Rows of different lengths also give `GR_DOMAIN`, and a string
that is not a list of rows gives `GR_UNABLE`. The string `[]` has zero
rows and an unspecified number of columns: with `resize` it gives a
$0 \times 0$ matrix, and without it is accepted for any $0 \times c$
matrix. A matrix with $n > 0$ rows and no columns is written `[[], ...,
[]]`.
:::
:::

## Assignment and conversion {#assignment}

::: function
int gr_mat_zero(gr_mat_t res, gr_ctx_t ctx)
int gr_mat_one(gr_mat_t res, gr_ctx_t ctx)

Set `res` to the zero matrix, respectively to the matrix with ones on
the main diagonal and zeros elsewhere (the identity if square).
:::

::: function
int gr_mat_set(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_set_fmpz_mat(gr_mat_t res, const fmpz_mat_t mat, gr_ctx_t ctx)
int gr_mat_set_fmpq_mat(gr_mat_t res, const fmpq_mat_t mat, gr_ctx_t ctx)
int gr_mat_set_gr_mat_other(gr_mat_t res, const gr_mat_t mat, gr_ctx_t mat_ctx, gr_ctx_t ctx)

Set `res` to a copy of `mat`, converting the entries from integers,
rationals or the ring `mat_ctx` into $R$. A conversion that is not
possible gives `GR_DOMAIN` or `GR_UNABLE`. The shapes must agree.
:::

::: function
int gr_mat_set_scalar(gr_mat_t res, gr_srcptr c, gr_ctx_t ctx)
int gr_mat_set_ui(gr_mat_t res, ulong c, gr_ctx_t ctx)
int gr_mat_set_si(gr_mat_t res, slong c, gr_ctx_t ctx)
int gr_mat_set_fmpz(gr_mat_t res, const fmpz_t c, gr_ctx_t ctx)
int gr_mat_set_fmpq(gr_mat_t res, const fmpq_t c, gr_ctx_t ctx)

Set `res` to the matrix with $c$ on the main diagonal and zeros
elsewhere.
:::

## Comparisons and structure {#predicates}

::: function
truth_t gr_mat_equal(const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)

Returns whether the matrices are equal: `T_FALSE` if the shapes
differ, otherwise the conjunction of the entrywise comparisons.
:::

::: function
truth_t gr_mat_is_zero(const gr_mat_t mat, gr_ctx_t ctx)
truth_t gr_mat_is_one(const gr_mat_t mat, gr_ctx_t ctx)
truth_t gr_mat_is_neg_one(const gr_mat_t mat, gr_ctx_t ctx)
truth_t gr_mat_is_scalar(const gr_mat_t mat, gr_ctx_t ctx)

Return whether `mat` is zero, has $1$ or $-1$ on the main diagonal and
zeros elsewhere, or is a scalar matrix (a diagonal matrix with equal
diagonal entries).
:::

::: function
truth_t gr_mat_is_diagonal(const gr_mat_t mat, gr_ctx_t ctx)
truth_t gr_mat_is_upper_triangular(const gr_mat_t mat, gr_ctx_t ctx)
truth_t gr_mat_is_lower_triangular(const gr_mat_t mat, gr_ctx_t ctx)
truth_t gr_mat_is_hessenberg(const gr_mat_t mat, gr_ctx_t ctx)

Return whether all entries off the main diagonal, below it, above it,
or below the first subdiagonal (upper Hessenberg form) are zero. The
matrix need not be square.
:::

## Rows, columns and transposition {#rows}

::: function
int gr_mat_concat_horizontal(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
int gr_mat_concat_vertical(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)

Set `res` to $(A \mid B)$, respectively to $A$ stacked above $B$.
Return `GR_DOMAIN` if the numbers of rows (columns) of the inputs, or
the shape of `res`, do not match.
:::

::: function
int gr_mat_transpose(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)

Sets `res` to the transpose of `mat`. Aliasing is allowed for square
matrices, in which case the entries are swapped in place.
:::

::: function
int gr_mat_swap_rows(gr_mat_t mat, slong * perm, slong r, slong s, gr_ctx_t ctx)
void _gr_mat_swap_rows(gr_mat_t mat, slong * perm, slong r, slong s, gr_ctx_t ctx)
int gr_mat_swap_cols(gr_mat_t mat, slong * perm, slong r, slong s, gr_ctx_t ctx)
void _gr_mat_swap_cols(gr_mat_t mat, slong * perm, slong r, slong s, gr_ctx_t ctx)

Swap rows (columns) `r` and `s` of `mat`, and the entries `r` and `s`
of the permutation `perm` if it is not `NULL`. The non-underscore
versions return `GR_DOMAIN` for indices out of bounds; the underscore
versions do not check them.
:::

::: function
int gr_mat_invert_rows(gr_mat_t mat, slong * perm, gr_ctx_t ctx)
int gr_mat_invert_cols(gr_mat_t mat, slong * perm, gr_ctx_t ctx)

Reverse the order of the rows (columns): swap row $i$ and row $r - 1 -
i$ for $0 \le i < r/2$, and apply the same permutation to `perm` if it
is not `NULL`.
:::

::: function
int gr_mat_permute_rows(gr_mat_t mat, slong * perm_store, const slong * perm_act, gr_ctx_t ctx)
int gr_mat_permute_rows_inv(gr_mat_t mat, slong * perm_store, const slong * perm_act, gr_ctx_t ctx)
int gr_mat_permute_cols(gr_mat_t mat, slong * perm_store, const slong * perm_act, gr_ctx_t ctx)
int gr_mat_permute_cols_inv(gr_mat_t mat, slong * perm_store, const slong * perm_act, gr_ctx_t ctx)

Permute the rows (columns) by `perm_act`, or by its inverse for the
*inv* versions, and apply the same permutation to `perm_store` if it is
not `NULL`. `gr_mat_permute_rows` moves row `perm_act[i]` to row $i$
(and `perm_store[perm_act[i]]` to `perm_store[i]`); the *inv* version
moves row $i$ to row `perm_act[i]`, without computing the inverse
permutation. `perm_store` and `perm_act` may be the same array.

::: details
Swapping rows with `gr_mat_swap_rows` and `gr_mat_invert_rows` while
recording the swaps in a permutation that starts as the identity, and
then calling `gr_mat_permute_rows_inv` with that permutation, restores
the original matrix (and likewise for columns).
:::
:::

::: function
int gr_mat_move_row(gr_mat_t A, slong i, slong new_i, gr_ctx_t ctx)

Moves row $i$ to position $\mathit{new\_i}$, shifting the rows in
between by one: with $i = 4$ and $\mathit{new\_i} = 7$ the rows $A_4,
A_5, A_6, A_7$ become $A_5, A_6, A_7, A_4$, and with $i = 7$ and
$\mathit{new\_i} = 4$ they become $A_7, A_4, A_5, A_6$. Returns
`GR_DOMAIN` if an index is out of bounds and otherwise succeeds.
:::

## Entrywise operations {#entrywise}

::: function
int gr_mat_entrywise_unary_op(gr_mat_t res, gr_method_unary_op f, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_entrywise_binary_op(gr_mat_t res, gr_method_binary_op f, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
int gr_mat_entrywise_binary_op_scalar(gr_mat_t res, gr_method_binary_op f, const gr_mat_t mat, gr_srcptr c, gr_ctx_t ctx)

Set each entry of `res` to $f$ applied to the corresponding entry of
`mat`, to the entries of `mat1` and `mat2`, or to the entry of `mat`
and the scalar $c$. Return `GR_DOMAIN` if the shapes do not match.
These are declared only when `gr.h` is included before `gr_mat.h`.
:::

::: function
truth_t gr_mat_entrywise_unary_predicate_all(gr_method_unary_predicate f, const gr_mat_t mat, gr_ctx_t ctx)
truth_t gr_mat_entrywise_unary_predicate_any(gr_method_unary_predicate f, const gr_mat_t mat, gr_ctx_t ctx)
truth_t gr_mat_entrywise_binary_predicate_all(gr_method_binary_predicate f, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)

Return whether the predicate $f$ holds for all entries, for some entry,
or for all pairs of corresponding entries (`T_FALSE` if the shapes
differ).
:::

## Addition and scalar arithmetic {#arithmetic}

::: function
int gr_mat_neg(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_add(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
int gr_mat_sub(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)

Set `res` to $-A$, $A + B$ or $A - B$.
:::

::: function
int gr_mat_add_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t ctx)
int gr_mat_scalar_add(gr_mat_t res, gr_srcptr x, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_add_ui(gr_mat_t res, const gr_mat_t mat, ulong x, gr_ctx_t ctx)
int gr_mat_add_si(gr_mat_t res, const gr_mat_t mat, slong x, gr_ctx_t ctx)
int gr_mat_add_fmpz(gr_mat_t res, const gr_mat_t mat, const fmpz_t x, gr_ctx_t ctx)
int gr_mat_add_fmpq(gr_mat_t res, const gr_mat_t mat, const fmpq_t x, gr_ctx_t ctx)
int gr_mat_add_scalar_other(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
int gr_mat_scalar_other_add(gr_mat_t res, gr_srcptr x, gr_ctx_t x_ctx, const gr_mat_t mat, gr_ctx_t ctx)

Set `res` to $A + xI$ or $xI + A$, adding $x$ to the main diagonal.
:::

::: function
int gr_mat_sub_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t ctx)
int gr_mat_scalar_sub(gr_mat_t res, gr_srcptr x, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_sub_ui(gr_mat_t res, const gr_mat_t mat, ulong x, gr_ctx_t ctx)
int gr_mat_sub_si(gr_mat_t res, const gr_mat_t mat, slong x, gr_ctx_t ctx)
int gr_mat_sub_fmpz(gr_mat_t res, const gr_mat_t mat, const fmpz_t x, gr_ctx_t ctx)
int gr_mat_sub_fmpq(gr_mat_t res, const gr_mat_t mat, const fmpq_t x, gr_ctx_t ctx)
int gr_mat_sub_scalar_other(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
int gr_mat_scalar_other_sub(gr_mat_t res, gr_srcptr x, gr_ctx_t x_ctx, const gr_mat_t mat, gr_ctx_t ctx)

Set `res` to $A - xI$ or $xI - A$.
:::

::: function
int gr_mat_mul_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t ctx)
int gr_mat_scalar_mul(gr_mat_t res, gr_srcptr x, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_mul_ui(gr_mat_t res, const gr_mat_t mat, ulong x, gr_ctx_t ctx)
int gr_mat_mul_si(gr_mat_t res, const gr_mat_t mat, slong x, gr_ctx_t ctx)
int gr_mat_mul_fmpz(gr_mat_t res, const gr_mat_t mat, const fmpz_t x, gr_ctx_t ctx)
int gr_mat_mul_fmpq(gr_mat_t res, const gr_mat_t mat, const fmpq_t x, gr_ctx_t ctx)
int gr_mat_mul_scalar_other(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
int gr_mat_scalar_other_mul(gr_mat_t res, gr_srcptr x, gr_ctx_t x_ctx, const gr_mat_t mat, gr_ctx_t ctx)

Set `res` to $Ax$ or $xA$ (the order matters over noncommutative
rings).
:::

::: function
int gr_mat_div_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t ctx)
int gr_mat_div_scalar_other(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
int gr_mat_div_ui(gr_mat_t res, const gr_mat_t mat, ulong x, gr_ctx_t ctx)
int gr_mat_div_si(gr_mat_t res, const gr_mat_t mat, slong x, gr_ctx_t ctx)
int gr_mat_div_fmpz(gr_mat_t res, const gr_mat_t mat, const fmpz_t x, gr_ctx_t ctx)
int gr_mat_div_fmpq(gr_mat_t res, const gr_mat_t mat, const fmpq_t x, gr_ctx_t ctx)

Set `res` to $A / x$, dividing each entry with `gr_div` (`GR_DOMAIN` if
a quotient does not exist in $R$).
:::

::: function
int gr_mat_addmul_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr c, gr_ctx_t ctx)
int gr_mat_submul_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr c, gr_ctx_t ctx)

Set `res` to $\mathit{res} + Ac$, respectively $\mathit{res} - Ac$.
:::

::: function
int gr_mat_mul_diag(gr_mat_t res, const gr_mat_t A, const gr_vec_t D, gr_ctx_t ctx)
int gr_mat_diag_mul(gr_mat_t res, const gr_vec_t D, const gr_mat_t A, gr_ctx_t ctx)

Set `res` to $AD$ or $DA$, where $D$ is the diagonal matrix with the
entries of the vector `D` on its diagonal: column $j$ (row $i$) of $A$
is multiplied by $D_j$ ($D_i$).
:::

## Matrix multiplication {#mul}

::: function
int gr_mat_mul(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
int gr_mat_mul_generic(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)

Set `res` to the product $AB$. `gr_mat_mul` calls the ring method
`GR_METHOD_MAT_MUL`; the default `gr_mat_mul_generic` is classical
multiplication. Most rings with a dedicated matrix module override the
method. Returns `GR_DOMAIN` if the shapes are incompatible: `gr_mat_mul`
itself does not check them, so every ring method must (all methods in
FLINT do). See [](gr_mat/mul).
:::

::: function
int gr_mat_mul_classical(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)

Computes the $m \times n$ by $n \times p$ product with the classical
algorithm, as $mp$ dot products of length $n$ (`_gr_vec_dot`), with the
right factor transposed shallowly so that its columns are contiguous.
:::

::: function
int gr_mat_mul_strassen(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)

One level of Strassen's algorithm [@Str1969] with the 7-multiplication,
15-addition sequence of Winograd [@Win1971] in the order given by
Bodrato [@Bodrato2010], which saves further operations when squaring.
The seven products of half size are computed with `gr_mat_mul`, so the
recursion continues only if the ring method dispatches back to this
function above a cutoff; then the complexity is $O(n^{\log_2 7}) =
O(n^{2.81})$. See [](algorithms/linear-algebra/strassen-winograd).
:::

::: function
int gr_mat_mul_waksman(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
int gr_mat_mul_rosowski(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)

Compute $AB$ with about half the multiplications of the classical
algorithm over a commutative ring, at the cost of more additions,
which pays off when multiplications are much more expensive than
additions (large integers, polynomials). Waksman's algorithm
[@Waksman1970] (refining Winograd's inner product trick [@Win1968])
needs exact division by 2, implemented as `gr_mul_2exp_si` by $-1$.
Rosowski's algorithm [@Rosowski2023] needs no division and uses fewer
multiplications than Waksman's when the inner dimension is odd.
:::

## Powers and polynomial evaluation {#powers}

See also [](#functions) for non-integer powers.

::: function
int gr_mat_sqr(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)

Sets `res` to $A^2$ (currently `gr_mat_mul` of $A$ by itself).
:::

::: function
int gr_mat_pow_ui(gr_mat_t res, const gr_mat_t mat, ulong e, gr_ctx_t ctx)
int gr_mat_pow_si(gr_mat_t res, const gr_mat_t mat, slong e, gr_ctx_t ctx)
int gr_mat_pow_fmpz(gr_mat_t res, const gr_mat_t mat, const fmpz_t e, gr_ctx_t ctx)

Set `res` to $A^e$ for a square matrix, by repeated squaring in the
matrix ring (`gr_generic_pow_ui`, `gr_generic_pow_fmpz`); a negative
exponent requires an invertible matrix. $A^0$ is the identity. Returns
`GR_DOMAIN` if the matrices are not square of the same size.
:::

::: function
int _gr_mat_gr_poly_evaluate(gr_mat_t res, gr_srcptr poly, slong len, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_gr_poly_evaluate(gr_mat_t res, const gr_poly_t poly, const gr_mat_t mat, gr_ctx_t ctx)

Set `res` to $f(A)$ for a scalar polynomial $f$ and a square matrix
$A$, by rectangular splitting [@PS1973]: about $2\sqrt{\mathit{len}}$
matrix products and $O(\mathit{len})$ scalar multiplications.
:::

## Norms {#norms}

::: function
int gr_mat_norm_max(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_norm_1(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_norm_inf(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_norm_frobenius(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)

Set `res` to the max norm $\max_{i,j} |a_{ij}|$, the 1-norm
(largest absolute column sum) $\max_j \sum_i |a_{ij}|$, the
$\infty$-norm (largest absolute row sum) $\max_i \sum_j |a_{ij}|$, or
the Frobenius norm $(\sum_{i,j} |a_{ij}|^2)^{1/2}$. The ring must
provide absolute values, comparisons and (for the Frobenius norm) square
roots; the result is an element of $R$.
:::

## Triangular solving {#triangular}

::: function
int gr_mat_nonsingular_solve_tril(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
int gr_mat_nonsingular_solve_tril_generic(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
int gr_mat_nonsingular_solve_tril_classical(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
int gr_mat_nonsingular_solve_tril_recursive(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
int gr_mat_nonsingular_solve_triu(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx)
int gr_mat_nonsingular_solve_triu_generic(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx)
int gr_mat_nonsingular_solve_triu_classical(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx)
int gr_mat_nonsingular_solve_triu_recursive(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx)

Solve $LX = B$ for a lower triangular $L$, or $UX = B$ for an upper
triangular $U$. Only the relevant triangle is read. If `unit` is set,
the diagonal is taken to be all ones and is not read (it may hold other
data, such as the $U$ factor of a packed LU decomposition). Otherwise
the diagonal entries are inverted if possible and divided by otherwise;
a division that fails gives `GR_DOMAIN`, even if the system has a
solution.

::: details
The *classical* versions substitute row by row with one dot product
per entry of $X$. The *recursive* versions split the triangular matrix
into blocks and reduce the work to two half-size triangular solves and
one matrix product, so that they benefit from fast `gr_mat_mul`. The
default methods `GR_METHOD_MAT_NONSINGULAR_SOLVE_TRIL` and `..._TRIU`
use the classical version when $B$ has fewer than 10 rows or columns and
the recursive version otherwise; `fq`, `fq_nmod`, `fq_zech`, `mpn_mod`
and `nfloat` override them. $X$ may be aliased with $B$.
:::
:::

## LU decomposition and solving {#lu}

::: function
int gr_mat_find_nonzero_pivot(slong * pivot_row, gr_mat_t mat, slong start_row, slong end_row, slong column, gr_ctx_t ctx)
int gr_mat_find_nonzero_pivot_generic(slong * pivot_row, gr_mat_t mat, slong start_row, slong end_row, slong column, gr_ctx_t ctx)
int gr_mat_find_nonzero_pivot_large_abs(slong * pivot_row, gr_mat_t mat, slong start_row, slong end_row, slong column, gr_ctx_t ctx)

Look for an entry that is provably nonzero in column `column`, in rows
$\mathit{start\_row} \le i < \mathit{end\_row}$. On success, set
`*pivot_row` and return `GR_SUCCESS`. Return `GR_DOMAIN` if all these
entries are provably zero, and `GR_UNABLE` if none is provably nonzero
but some cannot be decided.

::: details
The *generic* version takes the first provably nonzero entry; the
*large_abs* version takes the one of largest absolute value (partial
pivoting, for floating-point and ball arithmetic). `gr_mat_find_nonzero_pivot`
calls the ring method `GR_METHOD_MAT_FIND_NONZERO_PIVOT`: `arf`, `acf`,
`nfloat` and `nfloat_complex` use partial pivoting, and `ca` tries hard
to decide zero tests.
:::
:::

::: function
int gr_mat_lu(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)
int gr_mat_lu_generic(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)
int gr_mat_lu_classical(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)
int gr_mat_lu_recursive(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)

Compute a generalised LU decomposition $PA = LU$ of an $m \times n$
matrix $A$ with row permutation $P$, writing the rank to `rank`, the
permutation to the array `P` of length $m$, and $L$ and $U$ packed into
`LU`. See [](gr_mat/solve).
:::

::: function
int gr_mat_fflu(slong * rank, slong * P, gr_mat_t LU, gr_ptr den, const gr_mat_t A, int rank_check, gr_ctx_t ctx)

Computes a fraction-free LU decomposition of $A$ by Bareiss' algorithm,
over an integral domain with exact division, writing the final pivot (a
minor of $A$) to `den`. See [](gr_mat/solve).
:::

::: function
int gr_mat_nonsingular_solve(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
int gr_mat_nonsingular_solve_generic(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
int gr_mat_nonsingular_solve_lu(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
int gr_mat_nonsingular_solve_fflu(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)

Solve $AX = B$ for a square matrix $A$. Return `GR_DOMAIN` if $A$ is
singular, even if the system has solutions, and if $X$ does not have
entries in $R$. `gr_mat_nonsingular_solve` calls the ring method
`GR_METHOD_MAT_NONSINGULAR_SOLVE`, whose default
`gr_mat_nonsingular_solve_generic` uses the LU version in fields (or
inexact rings) and the fraction-free LU version otherwise. See
[](gr_mat/solve).
:::

::: function
int gr_mat_nonsingular_solve_lu_precomp(gr_mat_t X, const slong * perm, const gr_mat_t LU, const gr_mat_t B, gr_ctx_t ctx)
int gr_mat_nonsingular_solve_fflu_precomp(gr_mat_t X, const slong * perm, const gr_mat_t LU, const gr_mat_t B, gr_ctx_t ctx)

Solve $AX = B$ given an LU or fraction-free LU decomposition of $A$ with
full rank, as computed by `gr_mat_lu` or `gr_mat_fflu`. The data are
assumed to be valid, which is not checked. The *fflu* version computes
$\mathit{den} \cdot X$, where `den` is the denominator returned by
`gr_mat_fflu`.
:::

::: function
int gr_mat_nonsingular_solve_den(gr_mat_t X, gr_ptr den, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
int gr_mat_nonsingular_solve_den_fflu(gr_mat_t X, gr_ptr den, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)

Solve $AX = B$ over the fraction field of an integral domain $R$:
compute $X$ and $\mathit{den} \in R$ such that $A X = \mathit{den} \cdot
B$, with $\mathit{den} = \pm \det A$ from fraction-free elimination
($\mathit{den} = 1$ for $0 \times 0$ matrices).
Return `GR_DOMAIN` if $A$ is singular. See [](gr_mat/solve).
:::

::: function
int gr_mat_solve_field(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)

Solves $AX = B$ for an arbitrary (possibly rectangular or singular)
$A$ over a field. If the system has several solutions, an arbitrary one
is returned; `GR_DOMAIN` means that it has none. See [](gr_mat/solve).
:::

## Determinant, permanent and trace {#det}

::: function
int gr_mat_det(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_det_generic(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx)
int gr_mat_det_generic_field(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx)
int gr_mat_det_generic_integral_domain(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx)
int gr_mat_det_berkowitz(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_det_cofactor(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_det_lu(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_det_fflu(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)

Set `res` to the determinant of a square matrix (`GR_DOMAIN` otherwise),
with the division-free Berkowitz algorithm, cofactor expansion, LU
decomposition or fraction-free LU decomposition. See [](gr_mat/det).
:::

::: function
int gr_mat_permanent(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx)
int gr_mat_permanent_generic(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx)
int gr_mat_permanent_cofactor(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx)
int gr_mat_permanent_ryser(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx)
int gr_mat_permanent_glynn(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx)
int gr_mat_permanent_glynn_threaded(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx)

Set `res` to the permanent $\operatorname{perm}(A) = \sum_{\sigma}
\prod_i a_{i,\sigma(i)}$ of a square matrix, the determinant without
signs. Computing it is #P-hard [@Val1979], and all algorithms take
exponential time. See [](gr_mat/det#permanent).
:::

::: function
int gr_mat_trace(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)

Sets `res` to the sum of the diagonal entries of a square matrix.
Returns `GR_DOMAIN` if the matrix is not square.
:::

## Rank, echelon form and nullspace {#rank}

::: function
int gr_mat_rank(slong * rank, const gr_mat_t A, gr_ctx_t ctx)
int gr_mat_rank_lu(slong * rank, const gr_mat_t A, gr_ctx_t ctx)
int gr_mat_rank_fflu(slong * rank, const gr_mat_t A, gr_ctx_t ctx)

Set `rank` to the rank of $A$, by LU decomposition or by fraction-free
LU decomposition. `gr_mat_rank` uses LU over finite fields (and finite
rings pretending to be fields), fraction-free LU over other integral
domains, which keeps the entries small over $\Z$ and $\mathbb{Q}$, and returns
`GR_DOMAIN` for rings that are known not to be integral domains, where
the rank is not defined here. The *lu* and *fflu* versions do not check
the ring; any failure of the decomposition (an impossible inverse, an
undecidable pivot) is reported with `GR_UNABLE`.
:::

::: function
int gr_mat_rref(slong * rank, gr_mat_t R, const gr_mat_t A, gr_ctx_t ctx)
int gr_mat_rref_lu(slong * rank, gr_mat_t R, const gr_mat_t A, gr_ctx_t ctx)
int gr_mat_rref_fflu(slong * rank, gr_mat_t R, const gr_mat_t A, gr_ctx_t ctx)

Set $R$ to the reduced row echelon form of $A$ (pivots equal to 1, zeros
above and below them) and `rank` to the rank, using LU decomposition or
fraction-free LU decomposition followed by back substitution. This
requires a field, or at least that all pivots are invertible.
`gr_mat_rref` uses the LU version.
:::

::: function
int gr_mat_rref_den(slong * rank, gr_mat_t R, gr_ptr den, const gr_mat_t A, gr_ctx_t ctx)
int gr_mat_rref_den_fflu(slong * rank, gr_mat_t R, gr_ptr den, const gr_mat_t A, gr_ctx_t ctx)

Compute $R$ and $\mathit{den}$ such that $R / \mathit{den}$ is the reduced
row echelon form of $A$ over the fraction field of an integral domain,
with entries of $R$ in the ring. The denominator is a common, not
necessarily minimal, denominator: the pivots of $R$ all equal
$\mathit{den}$.
:::

::: function
int gr_mat_nullspace(gr_mat_t X, const gr_mat_t A, gr_ctx_t ctx)
int gr_mat_nullspace_no_resize(slong * nullity, gr_mat_t X, const gr_mat_t A, gr_ctx_t ctx)

Set the columns of $X$ to a basis of the right nullspace $\{x : Ax =
0\}$ of the $m \times n$ matrix $A$, read off from the reduced row
echelon form: one basis vector per non-pivot column. `gr_mat_nullspace`
resizes $X$ to $n$ rows and as many columns as the nullity; the
*no_resize* version requires $X$ with $n$ rows and at least that many
columns, zero-pads the rest and writes the nullity to `nullity`.

::: details
The basis is not canonical. Over an integral domain that is not a field
the basis is computed from `gr_mat_rref_den` and has entries in the
ring; it is a basis over the fraction field. Over other rings the
result may be meaningless.
:::
:::

::: function
int gr_mat_nullspace_from_rref(gr_mat_t X, const gr_mat_t A, gr_srcptr Aden, slong rank, gr_ctx_t ctx)

Sets $X$ to a nullspace basis given a reduced row echelon form $A$ of
rank `rank`, as computed by `gr_mat_rref`, or by `gr_mat_rref_den` with
the denominator `Aden` (`NULL` for an ordinary rref). $X$ must have $n$
rows and $n - \mathit{rank}$ columns. Returns `GR_UNABLE` if a pivot
position cannot be identified because a zero test is undecidable.
:::

## Inverse and adjugate {#inverse}

::: function
int gr_mat_inv(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_inv_generic(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)

Sets `res` to the inverse of a square matrix in the ring. `gr_mat_inv`
calls the ring method `GR_METHOD_MAT_INV`, whose default
`gr_mat_inv_generic` uses `gr_mat_nonsingular_solve` with the identity
matrix. Returns `GR_DOMAIN`
if the matrix is shown to be not invertible over $R$ (this is not
reliable over rings with zero divisors) and `GR_UNABLE` if
invertibility cannot be decided. For the inverse over the fraction
field, use `gr_mat_nonsingular_solve_den` or `gr_mat_adjugate`.
:::

::: function
int gr_mat_adjugate(gr_mat_t adj, gr_ptr det, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_adjugate_cofactor(gr_mat_t adj, gr_ptr det, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_adjugate_charpoly(gr_mat_t adj, gr_ptr det, const gr_mat_t mat, gr_ctx_t ctx)

Set `adj` to the adjugate of $A$ and `det` to its determinant, so that
$\operatorname{adj}(A) A = A \operatorname{adj}(A) = \det(A) I$. The
*cofactor* version computes $n^2$ minors; the *charpoly* version
evaluates $(-1)^{n-1}(A^{n-1} + c_{n-1} A^{n-2} + \cdots + c_1 I)$ from
the characteristic polynomial $x^n + c_{n-1} x^{n-1} + \cdots + c_0$ by
rectangular splitting, with $O(\sqrt{n})$ matrix products. The default
uses cofactors for $n \le 5$ and the characteristic polynomial
otherwise. All versions are division-free when the characteristic
polynomial is.
:::

## Characteristic and minimal polynomials {#charpoly}

::: function
int _gr_mat_charpoly(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_charpoly(gr_poly_t res, const gr_mat_t mat, gr_ctx_t ctx)
int _gr_mat_charpoly_generic(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_charpoly_generic(gr_poly_t res, const gr_mat_t mat, gr_ctx_t ctx)

Set `res` to the characteristic polynomial $\det(xI - A)$ of a square
matrix (`GR_DOMAIN` otherwise). `_gr_mat_charpoly` calls the ring method
`GR_METHOD_MAT_CHARPOLY`; the default `_gr_mat_charpoly_generic` uses
Berkowitz's algorithm, which works over any commutative ring. Rings
with fast arithmetic override it: `nmod` and `fmpz_mod` use Danilevsky's
algorithm above dimension 8 (16 for `fmpz_mod`), falling back to
Berkowitz when it meets a non-invertible pivot; `fq` and its variants
use Berkowitz up to dimension 12 and Danilevsky above; `fmpz` and
`fmpq` call `fmpz_mat_charpoly` and `fmpq_mat_charpoly`; `mpn_mod` has
its own. The underscore versions write $n + 1$ coefficients. See
[](algorithms/linear-algebra/charpoly).
:::

::: function
int _gr_mat_charpoly_berkowitz(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_charpoly_berkowitz(gr_poly_t res, const gr_mat_t mat, gr_ctx_t ctx)

The characteristic polynomial by Berkowitz's division-free algorithm
[@Ber1984], with $O(n^4)$ ring operations (as $O(n)$ matrix-vector
products for each of $n$ leading principal submatrices). Since it
needs no divisions and no zero tests, it succeeds whenever the ring
arithmetic does.
:::

::: function
int _gr_mat_charpoly_danilevsky_inplace(gr_ptr res, gr_mat_t mat, gr_ctx_t ctx)
int _gr_mat_charpoly_danilevsky(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_charpoly_danilevsky(gr_poly_t res, const gr_mat_t mat, gr_ctx_t ctx)
int _gr_mat_charpoly_gauss(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_charpoly_gauss(gr_poly_t res, const gr_mat_t mat, gr_ctx_t ctx)
int _gr_mat_charpoly_householder(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_charpoly_householder(gr_poly_t res, const gr_mat_t mat, gr_ctx_t ctx)

The characteristic polynomial with $O(n^3)$ operations, by Danilevsky's
reduction to a block companion (Frobenius) form [@Dan1937], or by
reduction to Hessenberg form with Gaussian elimination or with
Householder reflections followed by `gr_mat_charpoly_from_hessenberg`.
The *inplace* version overwrites its input.

::: details
These algorithms divide and test for zero, so they need a field; the
Householder version also needs square roots and complex conjugation.
They return `GR_DOMAIN` or `GR_UNABLE` when a division or square root
is impossible or a zero test cannot be decided.
:::
:::

::: function
int _gr_mat_charpoly_faddeev(gr_ptr res, gr_mat_t adj, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_charpoly_faddeev(gr_poly_t res, gr_mat_t adj, const gr_mat_t mat, gr_ctx_t ctx)
int _gr_mat_charpoly_faddeev_bsgs(gr_ptr res, gr_mat_t adj, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_charpoly_faddeev_bsgs(gr_poly_t res, gr_mat_t adj, const gr_mat_t mat, gr_ctx_t ctx)

The characteristic polynomial by the Faddeev-LeVerrier algorithm
[@FadSom1972], which computes $M_k = A M_{k-1} + c_{n-k+1} I$ and $c_{n-k} =
-\operatorname{tr}(A M_k) / k$, and, if `adj` is not `NULL`, the
adjugate matrix at no extra cost. The *bsgs* version is the baby-step
giant-step variant of Preparata and Sarwate [@PreSar1978], with
$O(n^{3.5})$ operations instead of $O(n^4)$ and $\sqrt{n}$ temporary
matrices. Both divide by $1, \ldots, n$ and fail (with `GR_DOMAIN` or
`GR_UNABLE`) in small positive characteristic or without division.
:::

::: function
int _gr_mat_charpoly_from_hessenberg(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_charpoly_from_hessenberg(gr_poly_t res, const gr_mat_t mat, gr_ctx_t ctx)

The characteristic polynomial of a matrix in upper Hessenberg form
(not checked), by the recurrence on its leading principal submatrices,
with $O(n^3)$ operations and no divisions.
:::

::: function
int gr_mat_hessenberg(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_hessenberg_gauss(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
int gr_mat_hessenberg_householder(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)

Set `res` to an upper Hessenberg matrix similar to $A$, computed with
Gaussian elimination (similarity transformations by elementary
matrices) or with Householder reflections; the default uses Gaussian
elimination. They need divisions and zero tests, and the Householder
version also square roots and complex conjugation, and return
`GR_DOMAIN` or `GR_UNABLE` if these fail.
:::

::: function
int gr_mat_minpoly_field(gr_poly_t res, const gr_mat_t mat, gr_ctx_t ctx)

Sets `res` to the minimal polynomial of a square matrix over a field,
by computing the minimal polynomials of Krylov sequences $v, Av, A^2 v,
\ldots$ of unit vectors with incremental Gaussian elimination
(`gr_mat_reduce_row`) and taking their least common multiple, with
$O(n^3)$ field operations. See
[](algorithms/linear-algebra/charpoly#minpoly).
:::

::: function
int _gr_mat_companion(gr_mat_t res, gr_srcptr poly, gr_ctx_t ctx)
int gr_mat_companion(gr_mat_t res, const gr_poly_t poly, gr_ctx_t ctx)

Sets the $n \times n$ matrix `res` to the companion matrix of a
polynomial of degree $n$, whose characteristic polynomial is the monic
associate of the polynomial. The leading coefficient must be invertible
(not checked). The underscore version reads $n + 1$ coefficients.
:::

::: function
int _gr_mat_companion_fraction(gr_mat_t res_num, gr_ptr res_den, gr_srcptr poly, gr_ctx_t ctx)
int gr_mat_companion_fraction(gr_mat_t res_num, gr_ptr res_den, const gr_poly_t poly, gr_ctx_t ctx)

Set `res_num` and the scalar `res_den` such that $\mathit{res\_num} /
\mathit{res\_den}$ is the companion matrix of a polynomial of degree
$n$, without dividing: `res_den` is the leading coefficient. The
underscore version reads $n + 1$ coefficients.
:::

::: function
int gr_mat_apply_row_similarity(gr_mat_t M, slong r, gr_ptr d, gr_ctx_t ctx)

Replaces the $n \times n$ matrix $M$ by $P^{-1} M P$ in place, where $P$
is the identity matrix with the zero entries of row $r$ replaced by $d$
(an elementary similarity transformation). Requires $0 \le r < n$;
returns `GR_DOMAIN` if $M$ is not square or $r$ is out of range. Such
transformations preserve the determinant and the characteristic and
minimal polynomials.
:::

::: function
int gr_mat_reduce_row(slong * column, gr_mat_t A, slong * P, slong * L, slong m, gr_ctx_t ctx)
int gr_mat_reduce_row_generic(slong * column, gr_mat_t A, slong * P, slong * L, slong m, gr_ctx_t ctx)

Reduces row $m$ of $A$ by the previous rows, which are assumed to be in
Gauss form though not necessarily in order, as needed for incremental
elimination (minimal polynomials). Entry $i$ of `P` is the row of $A$
with a pivot in column $i$, or $-1$; entry $i$ of `L` is the number of
possibly nonzero columns of row $i$ (all equal to the number of columns
unless $A$ is zero on the right), with nondecreasing entries. On return,
`column` is the column of the new pivot of row $m$ (the first available
column from the left, also recorded in `P`), or $-1$ if the row reduces
to zero.

::: details
`gr_mat_reduce_row` calls the ring method `GR_METHOD_MAT_REDUCE_ROW`;
`fmpz_mod`, `mpn_mod`, `fq` and `fq_nmod` override the default to
delay modular reductions.
:::
:::

## Eigenvalues and diagonalisation {#eigenvalues}

::: function
int gr_mat_eigenvalues(gr_vec_t lambda, fmpz_vec_t mult, const gr_mat_t mat, int flags, gr_ctx_t ctx)
int gr_mat_eigenvalues_other(gr_vec_t lambda, fmpz_vec_t mult, const gr_mat_t mat, gr_ctx_t mat_ctx, int flags, gr_ctx_t ctx)

Find the eigenvalues of a square matrix in the ring `ctx`, without
repetition, with their multiplicities in `mult`, as the roots of the
characteristic polynomial (`gr_poly_roots`, respectively
`gr_poly_roots_other` for a matrix over another ring `mat_ctx`). As
for `gr_poly_roots`, the multiplicities may add up to less than $n$ if
the ring is not algebraically closed. Over exact real and complex
numbers, see [](algorithms/exact/linear-algebra).
:::

::: function
int gr_mat_diagonalization(gr_vec_t D, gr_mat_t L, gr_mat_t R, const gr_mat_t A, int flags, gr_ctx_t ctx)
int gr_mat_diagonalization_generic(gr_vec_t D, gr_mat_t L, gr_mat_t R, const gr_mat_t A, int flags, gr_ctx_t ctx)
int gr_mat_diagonalization_precomp(gr_vec_t D, gr_mat_t L, gr_mat_t R, const gr_mat_t A, const gr_vec_t eigenvalues, const fmpz_vec_t mult, gr_ctx_t ctx)

Compute a diagonalisation $L A R = D$ of a square matrix: $D$ is the
vector of eigenvalues repeated according to multiplicity, the columns
of $R$ are right eigenvectors and the rows of $L = R^{-1}$ left
eigenvectors, so that $A = R D R^{-1}$. `L` or `R` may be `NULL` to skip
computing it.

::: details
Over a field, `GR_DOMAIN` means that $A$ is not diagonalisable over
the field. The generic version finds the eigenvalues with
`gr_mat_eigenvalues` and eigenvectors as nullspaces of $A - \lambda I$;
the *precomp* version takes eigenvalues and multiplicities computed
beforehand. `gr_mat_diagonalization` calls the ring method
`GR_METHOD_MAT_DIAGONALIZATION`, overridden by `acb` and `acf` with
numerical algorithms; over `acb`, an approximation from the QR algorithm
is certified with `acb_mat_eig_simple`
([](algorithms/ball-arithmetic/eigenvalues)).
:::
:::

## Jordan form and matrix functions {#functions}

::: function
int gr_mat_jordan_blocks(gr_vec_t lambda, slong * num_blocks, slong * block_lambda, slong * block_size, const gr_mat_t A, gr_ctx_t ctx)

Computes the Jordan block structure of a square matrix over a field
containing its eigenvalues: the distinct eigenvalues in `lambda`, the
number of blocks in `num_blocks`, and for each block the index of its
eigenvalue in `block_lambda` and its size in `block_size` (arrays of
length at least $n$). The block sizes are determined from the ranks of
the powers of $A - \lambda I$.
:::

::: function
int gr_mat_set_jordan_blocks(gr_mat_t mat, const gr_vec_t lambda, slong num_blocks, slong * block_lambda, slong * block_size, gr_ctx_t ctx)

Sets `mat` to the Jordan matrix with the given blocks (each block has
its eigenvalue on the diagonal and ones on the superdiagonal).
:::

::: function
int gr_mat_jordan_transformation(gr_mat_t mat, const gr_vec_t lambda, slong num_blocks, slong * block_lambda, slong * block_size, const gr_mat_t A, gr_ctx_t ctx)

Given the block structure computed by `gr_mat_jordan_blocks`, sets
`mat` to an invertible matrix $P$ with $A = P J P^{-1}$, where $J$ is the
Jordan matrix.
:::

::: function
int gr_mat_jordan_form(gr_mat_t J, gr_mat_t P, const gr_mat_t A, gr_ctx_t ctx)

Computes the Jordan form $J$ of a square matrix and, if `P` is not
`NULL`, a transformation with $A = P J P^{-1}$. The eigenvalues must lie
in the ring, which must be a field with decidable equality for the
result to be meaningful (exact fields such as `qqbar` and `ca`). The
algorithm and its use over exact numbers are described in
[linear algebra over exact numbers](algorithms/exact/linear-algebra#jordan).
:::

::: function
int gr_mat_func_jordan(gr_mat_t res, const gr_mat_t A, __gr_method_vec_op jet_func, gr_ctx_t ctx)
int gr_mat_func_param_jordan(gr_mat_t res, const gr_mat_t A, __gr_method_vec_scalar_op jet_func, gr_srcptr c, gr_ctx_t ctx)

Compute the matrix function $f(A) = P f(J) P^{-1}$ from the Jordan
decomposition $A = P J P^{-1}$. The function is given by `jet_func`,
with signature `int jet_func(gr_ptr r, gr_srcptr x, slong n, gr_ctx_t
ctx)`, which writes the Taylor coefficients $f(x), f'(x), \ldots,
f^{(n-1)}(x) / (n-1)!$ to `r`. The *param* version takes a function
`int jet_func(gr_ptr r, gr_srcptr x, slong n, gr_srcptr c, gr_ctx_t ctx)`
of $x$ and an extra parameter $c$, which may in fact point to arbitrary
data understood by `jet_func`.
:::

::: function
int gr_mat_exp(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx)
int gr_mat_exp_jordan(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx)
int gr_mat_log(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx)
int gr_mat_log_jordan(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx)

Set `res` to the exponential or the (principal) logarithm of a square
matrix. The default methods use the Jordan decomposition; `arb` and
`acb` override the exponential with `arb_mat_exp` and `acb_mat_exp`
(scaling and squaring of Taylor series with error bounds). The Jordan
versions need the eigenvalues in the ring and are meant for exact
fields such as `qqbar` and `ca`; elsewhere they typically return
`GR_UNABLE`.
:::

::: function
int gr_mat_pow_scalar(gr_mat_t res, const gr_mat_t A, gr_srcptr c, gr_ctx_t ctx)
int gr_mat_pow_scalar_jordan(gr_mat_t res, const gr_mat_t A, gr_srcptr c, gr_ctx_t ctx)
int gr_mat_pow_fmpq(gr_mat_t res, const gr_mat_t mat, const fmpq_t exp, gr_ctx_t ctx)
int gr_mat_pow_fmpq_jordan(gr_mat_t res, const gr_mat_t mat, const fmpq_t exp, gr_ctx_t ctx)

Set `res` to $A^c$ for a scalar or rational exponent, using the Jordan
decomposition (principal branch). The non-Jordan versions first check
for small integer exponents and use `gr_mat_pow_si` for them.
:::

::: function
int gr_mat_sqrt(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx)
int gr_mat_rsqrt(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx)

Set `res` to a square root $A^{1/2}$ or reciprocal square root
$A^{-1/2}$, computed from the Jordan decomposition (the only algorithm
implemented). This often fails with `GR_UNABLE` when the ring is not an
algebraically closed field, even if a square root exists.
:::

## Orthogonality and QR decomposition {#qr}

::: function
truth_t gr_mat_is_orthogonal(const gr_mat_t A, gr_ctx_t ctx)

Returns whether $A$ is square with $A A^T = A^T A = I$. The ring is
assumed to be commutative.
:::

::: function
truth_t gr_mat_is_row_orthogonal(const gr_mat_t A, gr_ctx_t ctx)
truth_t gr_mat_is_col_orthogonal(const gr_mat_t A, gr_ctx_t ctx)
truth_t gr_mat_is_row_orthonormal(const gr_mat_t A, gr_ctx_t ctx)
truth_t gr_mat_is_col_orthonormal(const gr_mat_t A, gr_ctx_t ctx)

Return whether $A A^T$ (row versions) or $A^T A$ (column versions) is
diagonal (orthogonal) or the identity (orthonormal).
:::

::: function
int gr_mat_randtest_orthogonal(gr_mat_t A, flint_rand_t state, gr_ctx_t ctx)

Sets the square matrix $A$ to a random orthogonal matrix, using
Cayley's construction $(I + S)^{-1}(I - S)$ with $S$ skew-symmetric, or a
permutation matrix if that fails, and then negates random rows. Returns `GR_DOMAIN` if $A$ is
not square. The ring is assumed to be commutative.
:::

::: function
int gr_mat_qr(gr_mat_t Q, gr_mat_t R, const gr_mat_t A, gr_ctx_t ctx)
int gr_mat_lq(gr_mat_t L, gr_mat_t Q, const gr_mat_t A, gr_ctx_t ctx)
int gr_mat_lq_generic(gr_mat_t L, gr_mat_t Q, const gr_mat_t A, gr_ctx_t ctx)
int gr_mat_lq_gso(gr_mat_t L, gr_mat_t Q, const gr_mat_t A, gr_ctx_t ctx)
int gr_mat_lq_recursive(gr_mat_t L, gr_mat_t Q, const gr_mat_t A, gr_ctx_t ctx)

Compute a reduced QR decomposition $A = QR$ of an $m \times n$ matrix of
full column rank ($m \ge n$; $Q$ is $m \times n$ with orthonormal
columns and $R$ is $n \times n$ upper triangular), or the transposed LQ
decomposition $A = LQ$ of a matrix of full row rank ($m \le n$; $L$ is
$m \times m$ lower triangular and $Q$ is $m \times n$ with orthonormal
rows).

::: details
The entries should be real (a subfield of $\mathbb{R}$ or real
floating-point numbers); for other rings the functions may succeed but
the normalisations need not be meaningful. A full QR decomposition with
a square $Q$ is obtained by extending $Q$ with an orthonormal basis of
the complement and padding $R$ (or $L$) with zeros. The *gso* algorithm
is Gram-Schmidt orthogonalisation of the rows and the *recursive*
algorithm a block recursion that benefits from fast matrix
multiplication. `gr_mat_qr` is a wrapper of the LQ decomposition, which
works on rows because the matrices are stored by rows.
`gr_mat_lq` calls the ring method `GR_METHOD_MAT_LQ` (overridden by
`nfloat`). Aliasing of $A$ and $Q$ is handled in place; aliasing of $A$
with $R$ or $L$ uses a temporary copy. The functions return `GR_DOMAIN`
for incompatible shapes, rank-deficient input, or when the divisions or
square roots of the normalisation are not possible in the ring.
:::
:::

## LLL reduction tests {#lll}

Let the rows $a_0, \ldots, a_{n-1}$ of $A$ be linearly independent vectors
over $\mathbb{R}$ with Gram-Schmidt orthogonalisation $b_0, \ldots,
b_{n-1}$ and coefficients $\mu_{i,j} = \langle a_i, b_j \rangle / \| b_j
\|^2$. The basis is LLL-reduced with parameters $(\delta, \eta)$, where
$1/4 < \delta \le 1$ and $1/2 \le \eta < \sqrt{\delta}$, if it is size
reduced, $|\mu_{i,j}| \le \eta$ for $0 \le j < i < n$, and satisfies the
Lovász condition $(\delta - \mu_{i,i-1}^2) \| b_{i-1} \|^2 \le \| b_i
\|^2$ for $1 \le i < n$ [@LLL1982]; see
[](algorithms/linear-algebra/lll).

::: function
truth_t gr_mat_is_row_lll_reduced_naive(const gr_mat_t A, gr_srcptr delta, gr_srcptr eta, gr_ctx_t ctx)
truth_t gr_mat_is_row_lll_reduced_with_removal_naive(const gr_mat_t A, gr_srcptr delta, gr_srcptr eta, gr_srcptr gs_B, slong newd, gr_ctx_t ctx)

Return whether the rows of $A$ are LLL-reduced with parameters
$(\delta, \eta)$, by Gram-Schmidt orthogonalisation in the ring and
checking the conditions row by row. The *with_removal* version checks
only the first `newd` rows, as needed when rows whose Gram-Schmidt
norm exceeds the bound `gs_B` are to be removed (as in `fmpz_lll`).

::: details
With interval arithmetic these functions stop early: `T_UNKNOWN` is
returned as soon as the test of one row is inconclusive, even if a
later row would show that the answer is `T_FALSE`.
:::
:::

## Special matrices {#special}

For the following functions, the output matrix must have the intended
shape.

::: function
int gr_mat_ones(gr_mat_t res, gr_ctx_t ctx)

Sets all entries of `res` to one.
:::

::: function
int gr_mat_pascal(gr_mat_t res, int triangular, gr_ctx_t ctx)

Sets `res` to a Pascal matrix of binomial coefficients: for
`triangular` 0 the symmetric matrix with entries $\binom{i+j}{i}$
(the rows of Pascal's triangle along the antidiagonals), for 1 the
upper triangular matrix with entries $\binom{j}{i}$, and for $-1$ the
lower triangular matrix with entries $\binom{i}{j}$.
:::

::: function
int gr_mat_stirling(gr_mat_t res, int kind, gr_ctx_t ctx)

Sets `res` to the matrix of Stirling numbers: unsigned of the first
kind for `kind` 0, signed of the first kind for 1, and of the second
kind for 2.
:::

::: function
int gr_mat_hilbert(gr_mat_t res, gr_ctx_t ctx)

Sets `res` to the Hilbert matrix with entries $1 / (i + j + 1)$, which
requires these inverses in $R$.
:::

::: function
int gr_mat_hadamard(gr_mat_t res, gr_ctx_t ctx)

Sets the square matrix `res` to a Hadamard matrix of its size if
possible. Returns `GR_DOMAIN` if none exists (the size must be 0, 1, 2
or a multiple of 4), and `GR_UNABLE` if the implementation cannot
construct one.

::: details
The matrix is constructed by `fmpz_mat_hadamard` and converted. It
uses Sylvester's doubling and Paley's constructions [@Pal1933], which
succeed for all sizes $2^e$ and $2^e(q + 1)$ with $q$ an odd prime
power. Sizes for which Hadamard matrices are known to exist but this
construction fails start with 92, 116, 156 (OEIS A046116).
:::
:::

## Random matrices {#random}

::: function
int gr_mat_randtest(gr_mat_t res, flint_rand_t state, gr_ctx_t ctx)

Sets `res` to a random matrix with a nonuniform distribution meant for
testing (sparse matrices and special values are favoured).
:::

::: function
int gr_mat_randops(gr_mat_t mat, flint_rand_t state, slong count, gr_ctx_t ctx)

Randomises `mat` in place by at most `count` random additions or
subtractions of one row (column) to another, which preserves the
determinant and the rank.
:::

::: function
int gr_mat_randpermdiag(int * parity, gr_mat_t mat, flint_rand_t state, gr_ptr diag, slong n, gr_ctx_t ctx)

Sets `mat` to a random row and column permutation of the diagonal
matrix whose first $n$ diagonal entries are given by the vector `diag`,
writing the parity (0 or 1) of the permutation to `parity`. Returns
`GR_DOMAIN` if the diagonal of `mat` is shorter than $n$.
:::

::: function
int gr_mat_randrank(gr_mat_t mat, flint_rand_t state, slong rank, gr_ctx_t ctx)

Sets `mat` to a random sparse matrix with the given rank and exactly
`rank` nonzero entries (meaningful over integral domains). A dense
matrix of the same rank is obtained by calling `gr_mat_randops`
afterwards.
:::

::: function
int gr_mat_randsimilar(gr_mat_t mat, flint_rand_t state, slong opcount, gr_ctx_t ctx)

Randomises the square matrix `mat` in place by at most `opcount`
conjugations by random elementary matrices, which preserve the
characteristic and minimal polynomials.
:::

## Linear differential equations {#ode}

::: function
int gr_mat_gr_poly_solve_lode_newton(gr_mat_t Y, const gr_mat_t A_numerator, const gr_poly_t A_denominator, const gr_mat_t Y0, slong len, gr_ctx_t A_poly_ctx, gr_ctx_t sol_poly_ctx)
int _gr_mat_gr_poly_solve_lode_newton(gr_mat_t Y, gr_mat_t Z, const gr_mat_t A_numerator, const gr_poly_t A_denominator, const gr_mat_t Y0, slong len, gr_ctx_t A_poly_ctx, gr_ctx_t sol_poly_ctx)
int _gr_mat_gr_poly_solve_lode_newton_start(gr_mat_t Y, gr_mat_t Z, gr_poly_t A_denominator_inv, const gr_mat_t A_numerator, const gr_poly_t A_denominator, const gr_mat_t Y0, gr_ctx_t sol_poly_ctx)
int _gr_mat_gr_poly_solve_lode_newton_step(gr_mat_t Y, gr_mat_t Z, gr_poly_t A_denominator_inv, slong len, const gr_mat_t A_numerator, const gr_poly_t A_denominator, int A_is_companion, gr_ctx_t sol_poly_ctx)

Solve the system of linear differential equations $Y'(t) = A(t) Y(t)$
with $Y(0) = Y_0$, where $A(t)$ is a matrix of rational functions given
by a matrix of polynomials `A_numerator` (over the polynomial ring
`A_poly_ctx`) and a common denominator `A_denominator` with nonzero
constant term, and set $Y$ to the power series solution truncated to
length `len`. The entries of $Y$ are polynomials over the ring
`sol_poly_ctx`.

::: details
The algorithm is the Newton iteration of Bostan, Chyzak, Ollivier,
Salvy, Schost and Sedoglavic [@BCOSSS2007], which doubles the precision
of $Y$ and of an approximation $Z$ of $Y^{-1}$ at each step, combined
with a Newton iteration for the power series inverse of the
denominator and a dedicated multiplication by $A(t) \bmod t^m$ when $A$
is a companion matrix. The precisions follow an optimal sequence
ending at `len`, which need not be a power of two. The *start* and
*step* functions perform the initialisation and one step to length
`len` (`A_is_companion` selects the companion matrix multiplication),
for callers that drive the iteration themselves.
:::
:::

## Internal functions {#internal}

::: function
int _gr_mat_check_resize(gr_mat_t mat, slong r, slong c, gr_ctx_t ctx)

Helper for the contexts of matrix rings and spaces (here `ctx` is such a
context, not the coefficient ring): resizes `mat` to $r \times c$ if the
context allows matrices of all sizes, and otherwise returns `GR_DOMAIN`
unless $r \times c$ is the shape of the context.
:::

## Test functions {#testing}

The following functions run `iters` random tests of an implementation
of a method, with matrices of size up to `maxn`. If `ctx` is `NULL`, a
random ring is chosen for each iteration; otherwise the given ring is
used.

::: function
void gr_mat_test_mul(gr_method_mat_binary_op mul_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
void gr_mat_test_lu(gr_method_mat_lu_op lu_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
void gr_mat_test_det(gr_method_mat_unary_op_get_scalar det_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
void gr_mat_test_charpoly(gr_method_mat_unary_op_get_scalar charpoly_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
void gr_mat_test_nonsingular_solve_tril(gr_method_mat_binary_op_with_flag solve_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
void gr_mat_test_nonsingular_solve_triu(gr_method_mat_binary_op_with_flag solve_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)

Test an implementation of `gr_mat_mul`, `gr_mat_lu`, `gr_mat_det`,
`_gr_mat_charpoly`, `gr_mat_nonsingular_solve_tril` or
`gr_mat_nonsingular_solve_triu` by comparing with the generic
algorithms and checking algebraic identities.
:::

::: function
void gr_mat_test_approx_mul_max_norm(gr_method_mat_binary_op mul_impl, gr_srcptr rel_tol, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
void gr_mat_test_approx_mul_pos_entrywise_accurate(gr_method_mat_binary_op mul_impl, gr_srcptr rel_tol, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)

Test the accuracy of a matrix multiplication over an approximate
(floating-point) ring: the first checks $\|C - AB\| \le \|A\| \|B\|
\cdot \mathit{rel\_tol}$ in the max norm against classical
multiplication; the second generates matrices with nonnegative entries
and checks that every entry has relative error at most
$\mathit{rel\_tol}$.
:::
