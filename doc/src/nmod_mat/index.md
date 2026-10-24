---
title: nmod_mat.h
header: nmod_mat.h
description: Dense matrices over $\Z/n\Z$ for word-size $n$, with fast multiplication, LU decomposition, linear solving, echelon forms and characteristic polynomials.
extra-headers: nmod_poly
category: Matrices
order: 4
---

This module implements dense matrices over $\Z/n\Z$ for a word-size modulus
$n \ge 1$. An `nmod_mat_t` stores its entries as `ulong` values in
$[0, n)$ in row-major order, together with its dimensions and the
modulus as an `nmod_t`. Matrix multiplication, the operation on which
everything else rests, chooses among classical multiplication with
delayed reduction, Strassen–Winograd, vectorised kernels for moduli up to
$2^{52}$, bit- and byte-packed kernels for moduli up to 255 and
floating-point `gemm`; see [](nmod_mat/mul). On top of it the module
provides LU decomposition, linear system solving, determinants, inverses,
echelon forms, nullspaces and characteristic and minimal polynomials, most
of which require $n$ to be prime. The multimodular algorithms for
matrices over $\Z$ in [](fmpz_mat) are built on this module.

The guide [](nmod_mat/linear-algebra) explains which functions work for
which moduli and how they fit together.

::: summary algorithms/linear-algebra/matrix-multiplication :::

## Types and conventions {#types}

::: type
nmod_mat_struct
nmod_mat_t

An `nmod_mat_struct` has the fields `entries` (a pointer to the entry in
row 0, column 0), `r` and `c` (the numbers of rows and columns), `stride`
(the distance in words between the starts of consecutive rows) and `mod`
(the modulus as an `nmod_t`). An `nmod_mat_t` is an array of one
`nmod_mat_struct`, so that it is passed by reference.

::: details
Entry $(i, j)$ is `entries[i * stride + j]`. A matrix created by
`nmod_mat_init` owns a block of $rc$ words with `stride` equal to $c$; a
window created by `nmod_mat_window_init` points into another matrix and
has that matrix's stride. All functions except `nmod_mat_clear` and
`nmod_mat_swap` accept windows, so code that accesses entries must always
use the stride to find a row. There is no array of row pointers: rows are
swapped by swapping their entries.
:::
:::

The following conventions apply throughout the module.

- *Reduced entries.* Unless stated otherwise, the entries of input
  matrices must lie in $[0, n)$, and output entries do. Entries may be
  written directly with `nmod_mat_entry`; it is then the caller's
  responsibility to keep them reduced. Nothing is checked.
- *Dimensions.* The output must be initialised with dimensions compatible
  with the operation; most functions do not check this (some check it
  with `FLINT_ASSERT` only). Matrices with zero rows or columns are
  allowed everywhere.
- *Modulus.* All matrices passed to one function must have the same
  modulus. Any $n \ge 1$ is allowed for memory management, arithmetic,
  multiplication, determinants, characteristic polynomials and the Howell
  form. LU decomposition, solving, inversion, rank, echelon forms,
  nullspaces and minimal polynomials divide by pivots and assume that $n$
  is prime. This is not checked: a pivot that is not invertible raises a
  `FLINT_IMPINV` exception, and other results may be silently wrong. See
  [](nmod_mat/linear-algebra#moduli).
- *Aliasing.* Each function states whether the output may be the same
  matrix as an input. Partially overlapping windows are never allowed.

## Memory management {#memory}

::: function
void nmod_mat_init(nmod_mat_t mat, slong rows, slong cols, ulong n)

Initialises `mat` as a `rows` by `cols` matrix modulo $n$ with all
entries zero. Requires $n \ge 1$. Throws an exception if
$\mathit{rows} \cdot \mathit{cols}$ overflows an `slong`.
:::

::: function
void nmod_mat_init_set(nmod_mat_t mat, const nmod_mat_t src)

Initialises `mat` as a copy of `src`, with the same dimensions and
modulus. The copy owns its entries even if `src` is a window.
:::

::: function
void nmod_mat_clear(nmod_mat_t mat)

Frees the memory used by `mat`. Must not be called on a window.
:::

::: function
void nmod_mat_set(nmod_mat_t B, const nmod_mat_t A)

Copies the entries of $A$ to $B$, which must have the same dimensions.
The modulus of $B$ is not changed.
:::

::: function
void nmod_mat_swap(nmod_mat_t mat1, nmod_mat_t mat2)

Swaps `mat1` and `mat2` in constant time by swapping their structs. Must
not be used with windows, whose entries belong to another matrix.
:::

::: function
void nmod_mat_swap_entrywise(nmod_mat_t mat1, nmod_mat_t mat2)

Swaps the entries of two matrices of the same dimensions one by one. This
works for windows.
:::

::: function
nmod_t nmod_mat_mod(const nmod_mat_t mat)

Returns the modulus of `mat`.
:::

::: function
void nmod_mat_set_mod(nmod_mat_t mat, ulong n)

Sets the modulus of `mat` to $n \ge 1$. The entries are not changed, so
they must be reduced afterwards if they may exceed $n - 1$.
:::

## Entries and dimensions {#entries}

::: macro
nmod_mat_entry(mat, i, j)

The entry in row $i$ and column $j$ of `mat` (indexed from zero), as an
lvalue that can be read and assigned. No bounds checking is performed.
:::

::: function
ulong nmod_mat_get_entry(const nmod_mat_t mat, slong i, slong j)
void nmod_mat_set_entry(nmod_mat_t mat, slong i, slong j, ulong x)

Returns, respectively sets, the entry in row $i$ and column $j$. `x` must
be reduced.
:::

::: function
ulong * nmod_mat_entry_ptr(const nmod_mat_t mat, slong i, slong j)
ulong * nmod_mat_row_ptr(const nmod_mat_t mat, slong i)

Returns a pointer to the entry in row $i$ and column $j$, respectively to
the first entry of row $i$. The entries of a row are contiguous; the next
row starts `mat->stride` words later.
:::

::: function
slong nmod_mat_nrows(const nmod_mat_t mat)
slong nmod_mat_ncols(const nmod_mat_t mat)

Returns the number of rows, respectively columns, of `mat`.
:::

::: function
int nmod_mat_is_empty(const nmod_mat_t mat)
int nmod_mat_is_square(const nmod_mat_t mat)

Returns 1 if `mat` has no rows or no columns, respectively if it has as
many rows as columns, and 0 otherwise.
:::

## Windows {#windows}

::: function
void nmod_mat_window_init(nmod_mat_t window, const nmod_mat_t mat, slong r1, slong c1, slong r2, slong c2)

Initialises `window` as the submatrix of `mat` formed by rows $r_1$ to
$r_2 - 1$ and columns $c_1$ to $c_2 - 1$, where
$0 \le r_1 \le r_2 \le$ `mat->r` and $0 \le c_1 \le c_2 \le$ `mat->c`.
The window shares its entries with `mat`, so writing to one changes the
other.

::: details
No memory is allocated. The window has the stride and modulus of `mat`
and can be passed to any function that takes an `nmod_mat_t`, as input or
output, as long as no two arguments overlap partially. The bounds are
only checked by `FLINT_ASSERT`. A window must not be passed to
`nmod_mat_clear` or `nmod_mat_swap`.
:::
:::

::: function
void nmod_mat_window_clear(nmod_mat_t window)

Releases a window. This does nothing at present, but every window should
be cleared for future compatibility.
:::

## Concatenation {#concat}

::: function
void nmod_mat_concat_vertical(nmod_mat_t res, const nmod_mat_t mat1, const nmod_mat_t mat2)
void nmod_mat_concat_horizontal(nmod_mat_t res, const nmod_mat_t mat1, const nmod_mat_t mat2)

Sets `res` to `mat1` placed above `mat2`, respectively to the left of
`mat2`. For an $m \times n$ matrix `mat1`, `mat2` must be $k \times n$
and `res` $(m + k) \times n$, respectively `mat2` $m \times k$ and `res`
$m \times (n + k)$. `res` must not be aliased with the inputs.
:::

## Special matrices and comparison {#comparison}

::: function
void nmod_mat_zero(nmod_mat_t mat)
void nmod_mat_one(nmod_mat_t mat)

Sets `mat` to the zero matrix, respectively to the matrix with ones on
the main diagonal and zeros elsewhere. `mat` need not be square.
:::

::: function
int nmod_mat_equal(const nmod_mat_t mat1, const nmod_mat_t mat2)

Returns 1 if the two matrices have the same dimensions and the same
entries, and 0 otherwise. The moduli are not compared.
:::

::: function
int nmod_mat_is_zero(const nmod_mat_t mat)
int nmod_mat_is_one(const nmod_mat_t mat)
int nmod_mat_is_zero_row(const nmod_mat_t mat, slong i)

Returns 1 if `mat` is zero, if it has ones on the main diagonal and zeros
elsewhere, or if its row $i$ is zero, respectively, and 0 otherwise. Empty
matrices are both zero and one.

::: details
For $n = 1$ every matrix is zero and one: `nmod_mat_is_one` compares the
diagonal entries with $1 \bmod n$, which is 0, as `nmod_mat_one` sets
them.
:::
:::

## Random matrices {#random}

These functions are intended for testing.

::: function
void nmod_mat_rand(nmod_mat_t mat, flint_rand_t state)

Sets the entries to independent integers chosen uniformly at random in
$[0, n)$.
:::

::: function
void nmod_mat_randtest(nmod_mat_t mat, flint_rand_t state)

Sets each row with `_nmod_vec_randtest`: entries from `n_randtest`
reduced modulo $n$, and with probability 1/2 a sparse row.
:::

::: function
void nmod_mat_randfull(nmod_mat_t mat, flint_rand_t state)

Sets the entries to independent integers chosen uniformly at random in
$[1, n)$, so that no entry is zero.

::: details
The entries are $\max(1, x)$ with $x$ uniform in $[0, n)$, so 1 is twice
as likely as any other value. For $n = 1$ the entries are set to 0, the
only residue.
:::
:::

::: function
int nmod_mat_randpermdiag(nmod_mat_t mat, flint_rand_t state, nn_srcptr diag, slong n)

Sets `mat` to $P D Q$ where $D$ is the matrix with $d_0, \ldots, d_{n-1}$
(the entries of `diag`) on the main diagonal and zeros elsewhere, and $P$
and $Q$ are random permutation matrices. Requires $n \le$ the number of
rows and of columns. Returns the parity (0 for even, 1 for odd) of the
two permutations combined.
:::

::: function
void nmod_mat_randrank(nmod_mat_t mat, flint_rand_t state, slong rank)

Sets `mat` to a random matrix of the given rank, with exactly `rank`
nonzero entries, each uniform in $[1, n)$ and no two in the same row or
column. Throws an exception if `rank` is negative or exceeds a dimension.
Follow with `nmod_mat_randops` to obtain a dense matrix of the same rank.

::: details
For composite $n$ the nonzero entries need not be invertible, and the rank
over $\Z/n\Z$ is then not well defined; the matrix has `rank` nonzero
entries in distinct rows and columns.
:::
:::

::: function
void nmod_mat_randrank_profile(nmod_mat_t mat, flint_rand_t state, const slong * mrp_r, const slong * mrp_c, slong rank)

For a prime modulus, sets `mat` to a random matrix of rank $r$ = `rank`
whose rank profile matrix has its $r$ nonzero entries at the positions
(`mrp_r[k]`, `mrp_c[k]`) for $0 \le k < r$, given in any order. Throws an
exception if `rank` is negative or exceeds a dimension, or if two
positions share a row or a column or one lies outside the matrix.

::: details
The *rank profile matrix* of an $m \times n$ matrix $A$ of rank $r$ is
the unique $m \times n$ matrix $\mathcal{R}_A$ with $r$ entries equal to
1, no two in the same row or column, and zeros elsewhere, such that every
leading submatrix of $A$ (first $i$ rows, first $j$ columns) has the same
rank as the corresponding leading submatrix of $\mathcal{R}_A$. It
determines the row and column rank profiles of all leading submatrices of
$A$ [@DPS2017].

The matrix is computed as one product $L U$: $L$ is $m \times r$ with
column $k$ zero above row `mrp_r[k]`, 1 in that row and random below it,
and $U$ is $r \times n$ with row $k$ zero left of column `mrp_c[k]`,
uniformly random and nonzero in that column, and random to the right of
it. The random entries come from `nmod_mat_randtest`, so that the
factors, and therefore the result, are sparse with some probability.
Completing $L$ and $U$ into square unit lower and invertible upper
triangular matrices gives $A = L' \mathcal{R} U'$ where $\mathcal{R}$
has the entries 1 at the given positions, and multiplying by such
triangular matrices does not change the rank profile matrix. If the
modulus is 1, the matrix is zero.
:::
:::

::: function
void nmod_mat_randops(nmod_mat_t mat, flint_rand_t state, slong count)

Performs `count` random elementary operations on `mat`, each adding or
subtracting a row to or from another row, or a column to or from another
column (operations whose two indices coincide are skipped). These
operations have determinant 1, so the rank and, for square matrices, the
determinant are preserved.
:::

::: function
void nmod_mat_randtril(nmod_mat_t mat, flint_rand_t state, int unit)
void nmod_mat_randtriu(nmod_mat_t mat, flint_rand_t state, int unit)

Sets `mat` to a random lower, respectively upper, triangular matrix. The
entries strictly below (above) the diagonal are random, the diagonal
entries are 1 if `unit` is nonzero and random nonzero values otherwise.
`mat` need not be square.

::: details
A diagonal entry that comes out as 0 is replaced by 1. For composite $n$
the diagonal entries need not be invertible.
:::
:::

## Printing {#printing}

::: function
void nmod_mat_print_pretty(const nmod_mat_t mat)
int nmod_mat_fprint_pretty(FILE * file, const nmod_mat_t mat)

Prints a header line with the dimensions and modulus, then one line per
row, with the entries in brackets and right-aligned to the number of
decimal digits of $n$. `nmod_mat_fprint_pretty` returns a positive value
on success and a nonpositive value on a write error. Functions with a
`FILE` argument are declared only if `<stdio.h>` is included before
`nmod_mat.h`.

::: details
For example, a $2 \times 3$ matrix modulo 2903 is printed as

```
<2 x 3 integer matrix mod 2903>
[   0    0 2607]
[ 622    0    0]
```

For a matrix with no rows or no columns only the header line is printed.
:::
:::

::: function
int nmod_mat_print(const nmod_mat_t mat)
int nmod_mat_fprint(FILE * f, const nmod_mat_t mat)

The same as `nmod_mat_fprint_pretty`, writing to `stdout` and `f`
respectively.
:::

## Transposition and permutations {#permutations}

::: function
void nmod_mat_transpose(nmod_mat_t B, const nmod_mat_t A)

Sets $B$ to the transpose of $A$. Throws an exception unless $B$ has as
many rows as $A$ has columns and vice versa. $B$ may be aliased with $A$
if $A$ is square.

::: details
The transposition works on blocks of $32 \times 32$ entries, or, with AVX2
or NEON, on $4 \times 4$ blocks transposed in vector registers.
:::
:::

::: function
void _nmod_mat_transpose(nn_ptr B, slong Bstride, nn_srcptr A, slong Astride, slong m, slong n)

Writes the transpose of the $m \times n$ array $A$ (row stride `Astride`)
to the $n \times m$ array $B$ (row stride `Bstride`). The words are
copied as they are. The two arrays must not overlap.
:::

::: function
void nmod_mat_swap_rows(nmod_mat_t mat, slong * perm, slong r, slong s)
void nmod_mat_swap_cols(nmod_mat_t mat, slong * perm, slong r, slong s)

Swaps rows (columns) $r$ and $s$ of `mat`. If `perm` is not `NULL`,
entries $r$ and $s$ of the array `perm` are swapped as well.
:::

::: function
void nmod_mat_invert_rows(nmod_mat_t mat, slong * perm)
void nmod_mat_invert_cols(nmod_mat_t mat, slong * perm)

Reverses the order of the rows (columns) of `mat`, that is, swaps rows
(columns) $i$ and $r - 1 - i$ for $0 \le i < r/2$, where $r$ is the number
of rows (columns). If `perm` is not `NULL`, the same permutation is
applied to it.
:::

::: function
void nmod_mat_permute_rows(nmod_mat_t mat, const slong * perm_act, slong * perm_store)

Replaces row $i$ of `mat` by the old row `perm_act[i]`, for every $i$,
where `perm_act` is a permutation of the row indices. If `perm_store` is
not `NULL`, it is updated in the same way:
`perm_store[i]` becomes the old `perm_store[perm_act[i]]`.
:::

## Addition and scalar multiplication {#addition}

These functions allow the output to be aliased with any input.

::: function
void nmod_mat_add(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
void nmod_mat_sub(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)

Sets $C = A + B$, respectively $C = A - B$. The three matrices must have
the same dimensions.
:::

::: function
void nmod_mat_neg(nmod_mat_t B, const nmod_mat_t A)

Sets $B = -A$. The matrices must have the same dimensions.
:::

::: function
void nmod_mat_scalar_mul(nmod_mat_t B, const nmod_mat_t A, ulong c)

Sets $B = cA$. Requires $c < n$.

::: details
The cases $c = 0$, $1$ and $n - 1$ are a zeroing, a copy and a negation.
Otherwise, for $n < 2^{\bits - 1}$ and matrices with at least
`FLINT_MULMOD_SHOUP_THRESHOLD` entries, the products use
[Shoup multiplication](algorithms/arithmetic/shoup) with a quotient
precomputed for $c$, and `nmod_mul` otherwise.
:::
:::

::: function
void nmod_mat_scalar_mul_fmpz(nmod_mat_t B, const nmod_mat_t A, const fmpz_t c)

Sets $B = cA$ for any integer $c$, which is first reduced modulo $n$.
:::

::: function
void nmod_mat_scalar_addmul_ui(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, const ulong c)

Sets $C = A + cB$. Unlike `nmod_mat_scalar_mul`, $c$ may be any word; it
is reduced modulo $n$ first. Uses Shoup multiplication under the same
conditions as `nmod_mat_scalar_mul`.
:::

## Matrix multiplication {#mul}

`nmod_mat_mul` is the function to use: it selects the fastest of the
methods below for the dimensions, the modulus, the number of threads and
the instruction set. The other functions are the individual methods,
documented for testing, benchmarking and special needs. The kernels that
return an `int` return 0, without touching $C$, when they do not apply.
All methods work for any modulus within their range, prime or not.

::: summary algorithms/linear-algebra/matrix-multiplication :::

::: function
void nmod_mat_mul(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)

Sets $C = AB$, where $A$ is $m \times k$, $B$ is $k \times n$ and $C$ is
$m \times n$. $C$ may be aliased with $A$ or $B$.
:::

::: function
void nmod_mat_addmul(nmod_mat_t D, const nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
void nmod_mat_submul(nmod_mat_t D, const nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)

Sets $D = C + AB$, respectively $D = C - AB$. $D$ may be aliased with
$C$ but not with $A$ or $B$. The product is either computed by classical
multiplication fused with the addition, or with `nmod_mat_mul` and then
added; see [](nmod_mat/mul#addmul).
:::

::: function
void nmod_mat_mul_classical(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)

Sets $C = AB$ by classical multiplication, computing each entry as a dot
product with a single reduction. $C$ must not be aliased with $A$ or $B$.

::: details
The dot products use the parameters of `_nmod_vec_dot_params` for length
$k$ ([](algorithms/linear-algebra/dot-products)). The code takes one of
three forms.

- If the unreduced sums fit in one word and all dimensions exceed 10, $B$
  is transposed and several of its entries are packed into each word, as
  many as fit with $\lceil \log_2 (k (n-1)^2 + 1) \rceil$ bits each, so
  that one word product computes several dot products at once.
- If $k < 20$, or if $m < 20$ and neither $m \ge 4, k \ge 256$ nor
  $m \ge 2, k \ge 1000$, the dot products read the columns of $B$ in place
  with `NMOD_VEC_DOT`.
- Otherwise $B$ is transposed into a temporary matrix
  (`NMOD_MAT_MUL_TRANSPOSE_CUTOFF` is 20) and the dot products use
  `_nmod_vec_dot` on contiguous rows.
:::
:::

::: function
void nmod_mat_mul_classical_threaded(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)

The same as `nmod_mat_mul_classical`, with the rows of $C$ distributed
over the threads of FLINT's thread pool. Dimensions below 20 are handled
by the serial code. $C$ must not be aliased with $A$ or $B$.
:::

::: function
void nmod_mat_mul_strassen(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)

Sets $C = AB$ by one level of the Strassen–Winograd algorithm: seven
products of half-size blocks, each computed by `nmod_mat_mul` (which may
recurse), and 15 block additions. $C$ must not be aliased with $A$ or
$B$.

::: details
The operation schedule is Bodrato's variant [@Bodrato2010] arranged as
in [@DPZ2009], with two temporary matrices of about a quarter of the size
of $C$. Odd dimensions are handled by virtual padding: the blocks are
treated as if padded with a zero row or column, which costs no copies.
If any dimension is at most 4, `nmod_mat_mul` is called directly.

::: summary algorithms/linear-algebra/strassen-winograd :::
:::
:::

::: function
int nmod_mat_mul_blas(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)

Sets $C = AB$ by floating-point matrix multiplication with `flint_sgemm`
or `flint_dgemm`, which use an external BLAS if FLINT was built with one
and FLINT's own kernels otherwise. Returns 1 on success. Returns 0, and
does nothing, if a dimension is 0 or exceeds `INT_MAX`, if the inner
dimension is so large that no suitable set of at most 12 moduli exists,
or on 32-bit machines. $C$ may be aliased with $A$ or $B$.

::: details
The entries are lifted to the symmetric range $|x| \le \lfloor n/2
\rfloor$. Let $h = k \lfloor n/2 \rfloor^2$. If $h < 2^{24}$, a single
`flint_sgemm` is exact; if $h < 2^{53}$, a single `flint_dgemm`.
Otherwise the product is computed modulo up to 12 pairwise coprime
integers $t_i < \sqrt{2^{55}/k}$ with one `flint_dgemm` each, and
reconstructed by Chinese remaindering in mixed radix form. The
conversions are distributed over the thread pool.

::: summary algorithms/linear-algebra/modular-matrix-kernels :::
:::
:::

::: function
void nmod_mat_mul_u8(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)

Sets $C = AB$ for a modulus $n \le 255$ (checked only by `FLINT_ASSERT`)
using kernels for entries stored in bytes. $C$ may be aliased with $A$
or $B$. Uses several threads when available.

::: details
For $n = 2$ and $n = 3$ the entries are packed into bits (one bit per
entry for $n = 2$, two bit planes for $n = 3$) and multiplied by
Strassen–Winograd over the packed matrices with leaves computed by the
method of four Russians, unless $A$ has fewer than 64 (for $n = 2$) or
192 (for $n = 3$) rows. For $n \le 16$ the entries are converted to bytes
and multiplied by a kernel that replaces each multiplication by a table
lookup within a vector register, under Strassen–Winograd when all
dimensions are at least 256. For $17 \le n \le 255$ the products are
computed by `flint_sgemm` on blocks of at most
$\lfloor (2^{24} - 1)/(n - 1)^2 \rfloor$ columns of $A$, under
Strassen–Winograd on bytes when all dimensions are at least 2048, or by a
single `flint_sgemm` when several threads are available, $m, n \ge 256$
and $k (n - 1)^2 < 2^{24}$.

::: summary algorithms/linear-algebra/small-modulus-matrix-multiplication :::
:::
:::

::: function
void _nmod_mat_mul_u8(uint8_t * C, slong Cstride, const uint8_t * A, slong Astride, const uint8_t * B, slong Bstride, slong m, slong k, slong n, nmod_t mod)

The multiplication of `nmod_mat_mul_u8` on byte arrays: sets the
$m \times n$ array $C$ to $AB$, where $A$ is $m \times k$ and $B$ is
$k \times n$, each stored by rows with the given row strides (in bytes),
with entries reduced modulo `mod.n`, which must be at most 255. $C$ may
be equal to $A$ or $B$ (same pointer) but must not otherwise overlap
them.
:::

::: function
int nmod_mat_mul_u32(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)

Sets $C = AB$ for a modulus $n < 2^{32}$ with an integer SIMD kernel
(AVX-512, AVX2, NEON, or portable C), and returns 1. Returns 0 if
$n \ge 2^{32}$ or on 32-bit machines. $C$ may be aliased with $A$ or $B$.
Uses several threads when available.

::: details
The entries are lifted to signed 32-bit integers $|x| \le n/2$, and their
products are accumulated exactly in signed 64-bit lanes with widening
$32 \times 32 \to 64$-bit multiplications. When an accumulator could
overflow, it is folded: writing it as $h 2^{32} + l$, it is replaced by
$h c + l$ with $c = 2^{32} \bmod n$. For $n$ below about $2^{28.5}$ no
fold is needed within a block of 256 products, and a single pass serves
any $n < 2^{32}$.

::: summary algorithms/linear-algebra/modular-matrix-kernels :::
:::
:::

::: function
int _nmod_mat_mul_u32(uint32_t * C, slong Cstride, const uint32_t * A, slong Astride, const uint32_t * B, slong Bstride, slong m, slong k, slong n, nmod_t mod)

The multiplication of `nmod_mat_mul_u32` on arrays of 32-bit entries: $A$
is $m \times k$, $B$ is $k \times n$ and $C$ is $m \times n$, each stored
by rows with the given row strides, with entries reduced modulo `mod.n`.
Returns 0 if `mod.n` $\ge 2^{32}$ or on 32-bit machines and 1 otherwise.
$C$ may be equal to $A$ or $B$ (same pointer and stride) but must not
otherwise overlap them.
:::

::: function
int nmod_mat_mul_u52(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)

Sets $C = AB$ for a modulus $n \le 2^{52}$ with the AVX512-IFMA
instructions, and returns 1. Returns 0 if $n > 2^{52}$, or if FLINT was
not compiled for a processor with AVX512-IFMA (or on 32-bit machines).
$C$ may be aliased with $A$ or $B$. Uses several threads when available.

::: details
Each product of two entries below $2^{52}$ is accumulated as its low and
high 52-bit halves in two 64-bit lanes, which leaves room for 4095
products without any reduction. For $n \le 2^{26}$ the high halves vanish
and a single instruction per product suffices.
:::
:::

::: function
int nmod_mat_mul_k52(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)

Sets $C = AB$ for a modulus $n \le 2^{52}$ with integer SIMD instructions
(AVX-512, AVX2, NEON, or portable C), and returns 1. Returns 0 if
$n > 2^{52}$, on 32-bit machines, and for the portable C code compiled by
MSVC. $C$ may be aliased with $A$ or $B$. Uses several threads when
available.

::: details
Each entry, lifted to $|x| \le n/2$, is split into two limbs of 27 and 25
bits, and each product is formed as a two-limb Karatsuba product with
three widening $32 \times 32 \to 64$-bit multiplications, accumulated in
three 64-bit lanes without reduction for a block of 256 products. The
three sums are combined modulo $n$ in double precision.
:::
:::

::: function
int nmod_mat_mul_fp50(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)

Sets $C = AB$ for a modulus $n < 2^{50}$ in double precision, and returns
1. Returns 0 if $n \ge 2^{50}$, on 32-bit machines, and for the portable
C code compiled by MSVC. $C$ may be aliased with $A$ or $B$. Uses several
threads when available.

::: details
Each product is reduced as it is formed with the floating-point modular
multiplication of `fft_small` (one rounding and two fused multiply-adds),
and the results, of absolute value below $\frac98 n$, are summed in
double precision. The sums are reduced every
$\lfloor \frac89 (2^{53}/n - 1) \rfloor$ steps, which is every 6 steps
for $n$ close to $2^{50}$ and never within a block of 256 steps for
$n < 2^{53}/289 \approx 2^{44.8}$.
:::
:::

## Matrix-vector products {#matvec}

::: function
void nmod_mat_mul_nmod_vec(ulong * c, const nmod_mat_t A, const ulong * b, slong blen)
void nmod_mat_mul_nmod_vec_ptr(ulong * const * c, const nmod_mat_t A, const ulong * const * b, slong blen)

Sets $c = Ab$, where $A$ is $m \times k$ and $b$ is the vector of length
`blen`, truncated or extended by zeros to length $k$. The $m$ entries of
$c$ are written. In the `_ptr` version, entry $i$ of $b$ and $c$ is
`*b[i]` and `*c[i]`. Each entry of $c$ is a dot product with one
reduction. $c$ must not overlap $b$, except in the `_ptr` version, which
copies.
:::

::: function
void nmod_mat_nmod_vec_mul(ulong * c, const ulong * a, slong alen, const nmod_mat_t B)
void nmod_mat_nmod_vec_mul_ptr(ulong * const * c, const ulong * const * a, slong alen, const nmod_mat_t B)

Sets $c = aB$, where $B$ is $k \times n$ and $a$ is the vector of length
`alen`, truncated or extended by zeros to length $k$. The $n$ entries of
$c$ are written. In the `_ptr` version, entry $i$ of $a$ and $c$ is
`*a[i]` and `*c[i]`. $c$ must not overlap $a$, except in the `_ptr`
version, which copies.

::: details
The rows of $B$ are accumulated with delayed reductions where a vector
code path exists: with AVX512-IFMA for $n \le 2^{52}$ (as in
`nmod_mat_mul_u52`), with AVX2 or AVX-512 without IFMA for $n \le 2^{32}$
(splitting each product at bit 56 as in
[](algorithms/linear-algebra/dot-products#split); the next path is used
instead for short vectors: $k \le 6$, or $k \le 32$ and $B$ has at
least $\max(k, 16)$ columns), with a vector FMA
for $n < 2^{50}$ (as in `nmod_mat_mul_fp50`), and with AVX2 or AVX-512
for the other moduli, by accumulating the products of limbs of at most
32 bits of the entries as in
[](algorithms/linear-algebra/dot-products#split-limbs), with one
reduction of each entry of $c$ at the end
([](algorithms/linear-algebra/modular-matrix-kernels#thin)). These code
paths are not used when $B$ is too narrow for them, with thresholds on
the length $k$ of $a$ and the number of columns that depend on the path
(for instance, with AVX512-IFMA, from 16 columns for $k = 2$ down to 3 or
4 columns from $k = 6$ on; 5 columns for the moduli of the last path).
Otherwise, with one column, or up to 8 columns and $4 \le k \le 32$
($8 \le k \le 32$ for $n > 2^{62}$), each entry of $c$ is a dot product;
else each row is added with `_nmod_vec_scalar_addmul_nmod` when $k$ is
less than 6 (for $n \le 2^{32}$), 8 (up to $2^{62}$) or 10, or when $B$
has at least 2048 columns (4096 for $n \ge 2^{63}$), and the
product is done with `nmod_mat_mul_classical` on blocks of rows of $B$
for longer vectors.
:::
:::

## Powers and trace {#pow}

::: function
void nmod_mat_pow(nmod_mat_t dest, const nmod_mat_t mat, ulong pow)
void _nmod_mat_pow(nmod_mat_t dest, const nmod_mat_t mat, ulong pow)

Sets `dest` to `mat` raised to the power `pow`, by binary
exponentiation. The matrix must be square; the power 0 gives the identity
matrix. `nmod_mat_pow` allows `dest` to be aliased with `mat`;
`_nmod_mat_pow` does not.
:::

::: function
ulong nmod_mat_trace(const nmod_mat_t mat)

Returns the sum of the diagonal entries of the square matrix `mat`.
:::

## Triangular solving {#triangular}

These functions solve $LX = B$ or $UX = B$ for a square triangular matrix
with $r$ rows and an $r \times c$ matrix $B$, or, for the `_right`
versions, $XL = B$ or $XU = B$ for a $k \times r$ matrix $B$. Only the
triangle and, if `unit` is 0, the diagonal of $L$ or $U$ are read, so they
may be stored in a matrix holding other data, such as the output of
`nmod_mat_lu`. If `unit` is nonzero the diagonal is taken to be 1.
Otherwise every diagonal entry must be invertible modulo $n$, or a
`FLINT_IMPINV` exception is raised. $X$ may be aliased with $B$, but no
other aliasing is allowed.

::: summary algorithms/linear-algebra/gaussian-elimination :::

::: function
void nmod_mat_solve_tril(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit)
void nmod_mat_solve_triu(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit)

Sets $X = L^{-1} B$ for lower triangular $L$, respectively
$X = U^{-1} B$ for upper triangular $U$. Uses substitution if $B$ has
fewer than 64 rows or fewer than 64 columns, and the recursive algorithm
otherwise.
:::

::: function
void nmod_mat_solve_tril_classical(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit)
void nmod_mat_solve_triu_classical(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit)

Forward, respectively back, substitution: each column of $X$ is computed
entry by entry, each entry with one dot product and one reduction. The
inverses of the diagonal entries are computed once.
:::

::: function
void nmod_mat_solve_tril_recursive(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit)
void nmod_mat_solve_triu_recursive(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit)

Splits $L$ (or $U$) into blocks of $\lfloor r/2 \rfloor$ and
$\lceil r/2 \rceil$ rows and solves recursively, reducing the work to
matrix multiplications.

::: details
For lower triangular $L$, the block formula is

$$
\begin{pmatrix} A & 0 \\ C & D \end{pmatrix}^{-1}
\begin{pmatrix} X \\ Y \end{pmatrix} =
\begin{pmatrix} A^{-1} X \\ D^{-1} ( Y - C A^{-1} X ) \end{pmatrix},
$$

and for upper triangular $U$

$$
\begin{pmatrix} A & B \\ 0 & D \end{pmatrix}^{-1}
\begin{pmatrix} X \\ Y \end{pmatrix} =
\begin{pmatrix} A^{-1} (X - B D^{-1} Y) \\ D^{-1} Y \end{pmatrix}.
$$

The update of the right-hand side uses `nmod_mat_submul` and the two
half-size systems use `nmod_mat_solve_tril` (`nmod_mat_solve_triu`),
which recurses or switches to substitution.
:::
:::

::: function
void nmod_mat_solve_tril_right(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit)
void nmod_mat_solve_triu_right(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit)

Sets $X = B L^{-1}$ for lower triangular $L$, respectively $X = B U^{-1}$
for upper triangular $U$. Uses substitution if $B$ has fewer than 64 rows
or fewer than 64 columns, and the recursive algorithm otherwise.
:::

::: function
void nmod_mat_solve_tril_right_classical(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit)
void nmod_mat_solve_triu_right_classical(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit)

Substitution: each row of $X$ is computed entry by entry, from the last
to the first for $L$ and from the first to the last for $U$, each entry
with one dot product and one reduction. The columns of the triangle are
copied to rows 64 at a time, so that the dot products read contiguous
memory, and the inverses of the diagonal entries are computed once.
:::

::: function
void nmod_mat_solve_tril_right_recursive(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit)
void nmod_mat_solve_triu_right_recursive(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit)

Splits $L$ (or $U$) into blocks of $\lfloor r/2 \rfloor$ and
$\lceil r/2 \rceil$ rows and columns and solves recursively, reducing the
work to matrix multiplications.

::: details
For lower triangular $L$, the block formula is

$$
\begin{pmatrix} X & Y \end{pmatrix}
\begin{pmatrix} A & 0 \\ C & D \end{pmatrix}^{-1} =
\begin{pmatrix} (X - Y D^{-1} C) A^{-1} & Y D^{-1} \end{pmatrix},
$$

and for upper triangular $U$

$$
\begin{pmatrix} X & Y \end{pmatrix}
\begin{pmatrix} A & B \\ 0 & D \end{pmatrix}^{-1} =
\begin{pmatrix} X A^{-1} & (Y - X A^{-1} B) D^{-1} \end{pmatrix}.
$$

The update of the right-hand side uses `nmod_mat_submul` and the two
half-size systems use `nmod_mat_solve_tril_right`
(`nmod_mat_solve_triu_right`), which recurses or switches to
substitution.
:::
:::

## LU decomposition {#lu}

::: summary algorithms/linear-algebra/gaussian-elimination :::

::: function
slong nmod_mat_lu(slong * P, nmod_mat_t A, int rank_check)
slong nmod_mat_lu_classical(slong * P, nmod_mat_t A, int rank_check)
slong nmod_mat_lu_classical_delayed(slong * P, nmod_mat_t A, int rank_check)
slong nmod_mat_lu_recursive(slong * P, nmod_mat_t A, int rank_check)

Computes an LU decomposition with row pivoting of the $m \times n$ matrix
$A$ modulo a prime, overwriting $A$ with $L$ and $U$, and returns the rank
$r$. On return, row $i$ of $LU$ is row `P[i]` of the original $A$, where
`P` is an array of length $m$ (it need not be initialised). $U$ is in row
echelon form with $r$ nonzero rows, and $L$ is $m \times r$ unit lower
triangular. If `rank_check` is nonzero, the function may stop as soon as
it finds that the rank is less than $n$, and then returns 0 and leaves
$A$ undefined.

::: details
All four functions compute the same output. `nmod_mat_lu` chooses the
algorithm: `nmod_mat_lu_classical` for $\min(m, n) \le 3$ (or $\le 7$ when
$n \ge 2^{28}$), and otherwise `nmod_mat_lu_classical_delayed` below a
crossover that depends on the bit size of $n$ and
`nmod_mat_lu_recursive` above it. `_classical` is Gaussian elimination
with a reduction after each update, `_classical_delayed` accumulates the
updates of each entry in one to three words and reduces an entry only
when it is about to be used as a pivot or multiplier, and `_recursive`
splits the columns in half and reduces the work to triangular solving and
matrix multiplication. See [](nmod_mat/lu) for the storage format.
:::
:::

::: function
slong nmod_mat_lu_with_pivots(slong * P, slong * pivots_nonpivots, nmod_mat_t A)

The same as `nmod_mat_lu` with `rank_check` 0, and in addition writes to
`pivots_nonpivots`, an array of length $n$, the indices of the $r$ pivot
columns of $U$ in increasing order followed by the indices of the
$n - r$ other columns in increasing order.
:::

## PLE decomposition and rank profiles {#ple}

::: summary algorithms/linear-algebra/gaussian-elimination :::

::: function
slong nmod_mat_ple(slong * P, slong * pivots_nonpivots, nmod_mat_t A, int flags)
slong nmod_mat_ple_classical(slong * P, slong * pivots_nonpivots, nmod_mat_t A, int flags)
slong nmod_mat_ple_classical_delayed(slong * P, slong * pivots_nonpivots, nmod_mat_t A, int flags)
slong nmod_mat_ple_crout(slong * P, slong * pivots_nonpivots, nmod_mat_t A, int flags)
slong nmod_mat_ple_recursive(slong * P, slong * pivots_nonpivots, nmod_mat_t A, int flags)

Computes the PLE decomposition of the $m \times n$ matrix $A$ modulo a
prime that reveals its rank profile matrix, overwrites $A$ with $L$ and
$E$ in the storage of `nmod_mat_lu`, and returns the rank $r$. Row $i$
of $LE$ is row `P[i]` of the original $A$ (`P` has length $m$). The
array `pivots_nonpivots`, of length $n$, receives the pivot columns of
$E$ followed by the other columns, both in increasing order; it may be
`NULL`. The rows `P[r]`, ..., `P[m-1]` are also in increasing order, and
the rank profile matrix of $A$ has its nonzero entries at
(`P[k]`, `pivots_nonpivots[k]`) for $k < r$. If `flags` contains
`NMOD_MAT_PLE_RANK_CHECK` and the rank is less than $\min(m, n)$, the
function returns 0, possibly before the end, and leaves `A`, `P` and
`pivots_nonpivots` undefined.

::: details
These conditions determine the output, so all the algorithms return the
same; see [](nmod_mat/ple). `nmod_mat_ple_classical` is Gaussian
elimination with a reduction after each update, and
`nmod_mat_ple_classical_delayed` accumulates the updates of each entry
unreduced in one to three words, as `nmod_mat_lu_classical_delayed`.
`nmod_mat_ple_crout` computes each entry when it is needed, with dot
products and vector-matrix products that delay their reductions.
`nmod_mat_ple_recursive` splits the columns in half and reduces the work
to triangular solving and matrix multiplication, as
`nmod_mat_lu_recursive`. `nmod_mat_ple` chooses the algorithm.
:::
:::

::: function
void nmod_mat_ple_nullspace(nmod_mat_t X, const nmod_mat_t LU, const slong * pivots_nonpivots, slong rank)
void nmod_mat_ple_left_nullspace(nmod_mat_t X, const nmod_mat_t LU, const slong * P, slong rank)

Given the output `LU`, `P`, `pivots_nonpivots` and `rank` $= r$ of
`nmod_mat_ple` for an $m \times n$ matrix $A$, sets the
$n \times (n - r)$ matrix $X$ to the basis of the right nullspace
$\{x : Ax = 0\}$ computed by `nmod_mat_nullspace`, respectively the
$(m - r) \times m$ matrix $X$ to the basis of the left nullspace
$\{x : xA = 0\}$ computed by `nmod_mat_left_nullspace`. $X$ must have
exactly these dimensions.

::: details
Let $c_0 < \cdots < c_{r-1}$ be the pivot columns and
$f_0 < \cdots < f_{n-r-1}$ the other columns, as in `pivots_nonpivots`.
Vector $k$ of the right nullspace has a 1 in position $f_k$, zeros at
the other $f_j$, and entry $(j, k)$ of $C = -E_1^{-1} E_2$ in position
$c_j$, where $E_1$ and $E_2$ are the columns $c_j$ and $f_j$ of $E$; this
costs a triangular solve, about $r^2 (n - r)$ operations. Vector $k$ of
the left nullspace has a 1 in position `P[r+k]`, zeros at the other
`P[r+j]`, and entry $(k, j)$ of $D = -L_2 L_1^{-1}$ in position `P[j]`,
where $L_1$ and $L_2$ are the first $r$ and the last $m - r$ rows of $L$;
this costs a triangular solve on the right, about $(m - r) r^2$
operations. These are the only bases with an identity matrix in these
positions: the columns $f_j$ for the right nullspace, and for the left
nullspace the rows that are not in the row rank profile, in increasing
order.
:::
:::

::: function
void nmod_mat_ple_nullspace_compact(nmod_mat_t X, const nmod_mat_t LU, const slong * pivots_nonpivots, slong rank)
void nmod_mat_ple_left_nullspace_compact(nmod_mat_t X, slong * rows, const nmod_mat_t LU, const slong * P, slong rank)

The same bases without their identity matrix. For the right nullspace,
$X$ is $r \times (n - r)$ and receives $C$: its row $j$ gives the entries
in position `pivots_nonpivots[j]` and its column $k$ corresponds to the
vector with a 1 in position `pivots_nonpivots[r+k]`. For the left
nullspace, $X$ is $(m - r) \times r$, `rows` (length $m$) receives the row
rank profile of $A$ in increasing order followed by the other rows in
increasing order, and $X$ receives the columns of $D$ sorted accordingly:
its column $j$ gives the entries in position `rows[j]` and its row $k$
corresponds to the vector with a 1 in position `rows[r+k]`.
:::

::: macro
NMOD_MAT_PLE_RANK_CHECK

Flag for `nmod_mat_ple`: return 0 as soon as the rank is known to be
less than $\min(m, n)$.
:::

## Solving linear systems {#solve}

These functions require a prime modulus.

::: function
int nmod_mat_solve(nmod_mat_t X, const nmod_mat_t A, const nmod_mat_t B)

Solves $AX = B$ for a square matrix $A$. If $A$ is nonsingular, sets $X$
to the unique solution and returns 1. Otherwise returns 0 and leaves $X$
unchanged. $X$ may be aliased with $A$ or $B$.
:::

::: function
int nmod_mat_solve_vec(nn_ptr x, const nmod_mat_t A, nn_srcptr b)

Solves $Ax = b$ for a square $m \times m$ matrix $A$ and vectors $x$, $b$
of length $m$, as `nmod_mat_solve`. $x$ may be equal to $b$.
:::

::: function
int nmod_mat_can_solve(nmod_mat_t X, const nmod_mat_t A, const nmod_mat_t B)

Decides whether $AX = B$ has a solution for an arbitrary $m \times n$
matrix $A$ and $m \times c$ matrix $B$. If so, sets the $n \times c$
matrix $X$ to a solution and returns 1; otherwise sets $X$ to zero and
returns 0. Returns 0 if the dimensions are incompatible. $X$ may be
aliased with $B$ but not with $A$.
:::

::: function
int nmod_mat_can_solve_inner(slong * rank, slong * perm, slong * pivots, nmod_mat_t X, const nmod_mat_t A, const nmod_mat_t B)

The same as `nmod_mat_can_solve`, and in addition, for each pointer that
is not `NULL`, stores the rank of $A$ in `*rank`, the row permutation of
the LU decomposition of $A$ in `perm` (length $m$), and the indices of
the pivot columns in the first `*rank` entries of `pivots` (length at
least $\min(m, n)$).
:::

## Determinant, rank and inverse {#det}

::: function
ulong nmod_mat_det(const nmod_mat_t A)

Returns the determinant of the square matrix $A$, for any modulus. Throws
an exception if $A$ is not square.
:::

::: function
ulong nmod_mat_det_howell(const nmod_mat_t A)

Returns the determinant of the square matrix $A$, for any modulus, by
triangularising $A$ with Euclidean row operations. Throws an exception if
$A$ is not square.
:::

::: function
slong nmod_mat_rank(const nmod_mat_t A)

Returns the rank of $A$, which may have any shape. The modulus must be
prime.
:::

::: function
int nmod_mat_inv(nmod_mat_t B, const nmod_mat_t A)

If the square matrix $A$ is invertible, sets $B = A^{-1}$ and returns 1.
Otherwise returns 0, and the entries of $B$ are unspecified. The modulus
must be prime. $B$ may be aliased with $A$.
:::

## Row echelon form and nullspace {#rref}

These functions require a prime modulus.

::: function
slong nmod_mat_rref(nmod_mat_t A)

Puts $A$ in reduced row echelon form and returns its rank.
:::

::: function
slong nmod_mat_nullspace(nmod_mat_t X, const nmod_mat_t A)

Sets the first $d$ columns of $X$ to a basis of the right nullspace
$\{x : Ax = 0\}$ of the $m \times n$ matrix $A$, sets the other entries of
$X$ to zero, and returns the nullity $d = n - \operatorname{rank}(A)$. $X$
must have $n$ rows and at least $d$ columns; $n \times n$ always
suffices. $X$ may be aliased with $A$. The basis is the one with an
identity matrix in the rows of the columns of $A$ without pivot, computed
with `nmod_mat_ple` and `nmod_mat_ple_nullspace`.
:::

::: function
slong nmod_mat_left_nullspace(nmod_mat_t X, const nmod_mat_t A)

Sets the first $d$ rows of $X$ to a basis of the left nullspace
$\{v : vA = 0\}$ of the $m \times n$ matrix $A$, sets the other entries of
$X$ to zero, and returns the nullity $d = m - \operatorname{rank}(A)$. $X$
must have at least $d$ rows and $m$ columns; $m \times m$ always
suffices. The basis is the one with an identity matrix in the columns of
the rows of $A$ that are not in its row rank profile, computed with
`nmod_mat_ple` and `nmod_mat_ple_left_nullspace`.
:::

::: function
slong nmod_mat_reduce_row(nmod_mat_t M, slong * P, slong * L, slong m)

Reduces row $m$ of $M$ against rows that are already reduced, for
incremental Gaussian elimination, and returns the column of its new
pivot, or $-1$ if it reduces to zero.

::: details
`P` has one entry per column: `P[i]` is the index of the row whose pivot
(leading entry, which must be 1) is in column $i$, or $-1$ if there is
none. Row $m$ is reduced from left to right by subtracting multiples of
these rows. At the first column $i$ where it has a nonzero entry and
`P[i]` is $-1$, the row is divided by that entry, `P[i]` is set to $m$ and
$i$ is returned. If the row becomes zero, $-1$ is returned. `L[j]` is the
number of leading columns of row $j$ that may be nonzero; entries beyond
it are taken to be zero, which saves work when the rows get longer
towards the bottom. Setting every entry of `L` to the number of columns
is always valid. The pivot rows may be any rows of $M$, but at most $m$ of
them; the updates are accumulated without reduction in up to three words.
This is the elimination step of `nmod_mat_minpoly`.
:::
:::

## Howell form {#howell}

These functions work for any modulus. For a prime modulus they compute
the reduced row echelon form, up to the position of zero rows.

The algorithms are described in
[Howell form](algorithms/linear-algebra/gaussian-elimination#howell).

::: function
void nmod_mat_strong_echelon_form(nmod_mat_t A)

Puts $A$ in strong echelon form, using the algorithm of [@FieHof2014].
$A$ must have at least as many rows as columns; add zero rows otherwise.

::: details
The strong echelon form is an upper triangular form in which the
diagonal entry of row $i$ is the pivot of row $i$, a divisor of $n$, or
zero, the entries above each pivot are reduced modulo it, and every row
of the module spanned by $A$ whose first $i$ entries vanish is a
combination of the rows below row $i$. It is the Howell form
[@How1986; @StoMul1998] with zero rows inserted so that pivots lie on
the diagonal. [@FieHof2014] defines a lower left variant; FLINT computes
the upper right form.
:::
:::

::: function
slong nmod_mat_howell_form(nmod_mat_t A)

Puts $A$ in Howell form [@StoMul1998] and returns the number $k$ of
nonzero rows, which come first. $A$ must have at least as many rows as
columns.

::: details
The Howell form is computed as the strong echelon form, followed by moving
the zero rows to the bottom. Two matrices span the same submodule of
$(\Z/n\Z)^c$ if and only if their Howell forms are equal.
:::
:::

## Characteristic and minimal polynomials {#charpoly}

These functions are declared in `nmod_poly.h`. The polynomial must be
initialised with the same modulus as the matrix.

::: summary algorithms/linear-algebra/charpoly :::

::: function
void nmod_mat_charpoly(nmod_poly_t p, const nmod_mat_t M)
void nmod_mat_charpoly_berkowitz(nmod_poly_t p, const nmod_mat_t M)
int nmod_mat_charpoly_danilevsky(nmod_poly_t p, const nmod_mat_t M)

Sets $p$ to the characteristic polynomial $\det(xI - M)$ of the square
matrix $M$. An exception is thrown if $M$ is not square.
`nmod_mat_charpoly` and the division-free `_berkowitz` work for any
modulus. `_danilevsky` returns 1 on success and 0 if it needs an
inverse that does not exist, which cannot happen for a prime modulus.
For $n = 1$ the result is the zero polynomial.
:::

::: function
void nmod_mat_minpoly(nmod_poly_t p, const nmod_mat_t M)
void nmod_mat_minpoly_with_gens(nmod_poly_t p, const nmod_mat_t M, ulong * P)

Sets $p$ to the minimal polynomial of the square matrix $M$, the monic
polynomial of least degree with $p(M) = 0$, and throws an exception if
$M$ is not square. The modulus must be prime. If `P` is not `NULL`,
`nmod_mat_minpoly_with_gens` sets `P[j]` to 1 for each index $j$ such that
the unit vector $e_j$ was used as a Krylov generator, and leaves the other
entries of the array (of length $r$) unchanged.
:::

::: function
void nmod_mat_similarity(nmod_mat_t M, slong r, ulong d)

Replaces the square matrix $M$ by $P^{-1} M P$, where $P$ is the identity
matrix with the entries in row $r$ and every column $j \notin \{r-1, r\}$
replaced by $d$. Requires $d < n$.

::: details
Equivalently, $d$ times column $r$ is added to every column
$j \notin \{r - 1, r\}$, and then $d$ times the sum of the rows
$j \notin \{r - 1, r\}$ is subtracted from row $r$. Such transforms
preserve the determinant and the characteristic and minimal polynomials.
:::
:::

## Tuning parameters {#tuning}

These macros are defined in `nmod_mat.h`. The parameters of
`nmod_mat_mul` are in `flint-mparam.h`; see [](nmod_mat/mul#dispatch).

::: macro
NMOD_MAT_MUL_TRANSPOSE_CUTOFF

The dimension (20) from which classical multiplication transposes $B$
before forming dot products, and below which the threaded classical
multiplication runs serially.
:::

::: macro
NMOD_MAT_SOLVE_TRI_ROWS_CUTOFF
NMOD_MAT_SOLVE_TRI_COLS_CUTOFF

The triangular solving functions use substitution when the right-hand
side has fewer than `NMOD_MAT_SOLVE_TRI_ROWS_CUTOFF` (64) rows or fewer
than `NMOD_MAT_SOLVE_TRI_COLS_CUTOFF` (64) columns, and the recursive
algorithm otherwise.
:::

::: macro
NMOD_MAT_OPTIMAL_MODULUS_BITS

A suggested bit size, `FLINT_BITS - 5`, for the primes of multimodular
algorithms that rely on matrix multiplication, used for example by
`fmpz_mat_mul_multi_mod` and the Dixon solver of `fmpz_mat`. It was
chosen when classical multiplication with two-word dot products was the
fastest method; with the SIMD kernels, primes of at most 52 bits can be
cheaper per bit of modulus.
:::

## Internal functions {#internal}

These functions are used by the functions above and by other modules.

::: function
void _nmod_mat_mul_classical_op(nmod_mat_t D, const nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, int op)

Sets $D = AB$ if `op` is 0, $D = C + AB$ if `op` is 1 and $D = C - AB$ if
`op` is $-1$, by the classical algorithm of `nmod_mat_mul_classical`.
$D$ may be aliased with $C$ but not with $A$ or $B$. $C$ is not read if
`op` is 0.
:::

::: function
void _nmod_mat_mul_classical_threaded_pool_op(nmod_mat_t D, const nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, int op, thread_pool_handle * threads, slong num_threads)

The threaded version of `_nmod_mat_mul_classical_op`, using the
`num_threads` threads in `threads`, obtained from
`flint_request_threads`, in addition to the calling thread. If $k = 0$ or
the modulus is 1, $D$ is not written.
:::

::: function
void _nmod_mat_scalar_mul_generic(nmod_mat_t B, const nmod_mat_t A, ulong c)
void _nmod_mat_scalar_mul_precomp(nmod_mat_t B, const nmod_mat_t A, ulong c, ulong c_pr)

Sets $B = cA$ with `nmod_mul`, respectively with `n_mulmod_shoup` and
`c_pr = n_mulmod_precomp_shoup(c, n)`, which requires
$n < 2^{\bits - 1}$. Both require $c < n$.
:::

::: function
void _nmod_mat_scalar_addmul_ui_generic(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, const ulong c)
void _nmod_mat_scalar_addmul_ui_precomp(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, const ulong c, const ulong c_pr)

Sets $C = A + cB$, as the two functions above. Both require $c < n$.
:::

::: function
ulong _nmod_mat_det(nmod_mat_t A)

Returns the determinant of the square matrix $A$ modulo a prime, computed
by `nmod_mat_lu`, which destroys $A$.
:::

::: function
ulong _nmod_mat_det_howell(nmod_mat_t A)

Returns the determinant of the square matrix $A$ for any modulus, as
`nmod_mat_det_howell`, destroying $A$.
:::

::: function
slong _nmod_mat_rref(nmod_mat_t A, slong * pivots_nonpivots, slong * P)

Puts $A$ in reduced row echelon form and returns the rank $r$, as
`nmod_mat_rref`, and in addition stores in `pivots_nonpivots` (length the
number of columns) the pivot columns followed by the other columns, and in
`P` (length the number of rows) the row permutation of the LU
decomposition used.
:::

::: function
void _nmod_mat_charpoly_berkowitz(nn_ptr p, const nmod_mat_t M, nmod_t mod)

Writes the $r + 1$ coefficients of the characteristic polynomial of the
$r \times r$ matrix $M$ to `p`, by Berkowitz's algorithm. Works for any
modulus.
:::
