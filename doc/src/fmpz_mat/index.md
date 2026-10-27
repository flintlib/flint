---
title: fmpz_mat.h
header: fmpz_mat.h
description: Dense matrices over the integers, with multiplication, determinants, exact solving over $\Q$, echelon forms, Hermite and Smith normal forms, characteristic polynomials and lattice reduction.
extra-headers: d_mat
category: Matrices
order: 1
---

This module implements dense matrices with entries in $\Z$, stored as
arrays of `fmpz`. Besides arithmetic it provides the exact linear algebra
that is specific to the integers: matrix multiplication through several
multimodular and floating-point methods, determinants and linear system
solving by $p$-adic lifting and Chinese remaindering with output-sensitive
termination, fraction-free elimination, rank, echelon forms and
nullspaces, Hermite and Smith normal forms, characteristic and minimal
polynomials, and the classical LLL algorithms; the fast floating-point
LLL is in [](fmpz_lll).

Results that are rational, such as solutions of $AX = B$, inverses and
reduced row echelon forms, are returned as an integer matrix together with
a common denominator. The guide [](fmpz_mat/linear-algebra) explains this
representation, how the algorithms are chosen and how to organise
computations.

::: summary algorithms/linear-algebra/integer-matrix-multiplication :::

::: summary algorithms/linear-algebra/multimodular :::

## Types and conventions {#types}

::: type
fmpz_mat_struct
fmpz_mat_t

An `fmpz_mat_struct` has the fields `entries` (a pointer to the entry in
row 0, column 0), `r` and `c` (the numbers of rows and columns) and
`stride` (the distance, in entries, between the starts of consecutive
rows). An `fmpz_mat_t` is an array of one `fmpz_mat_struct`, so that it
is passed by reference.

::: details
Entry $(i, j)$ is `entries[i * stride + j]`. A matrix created by
`fmpz_mat_init` owns $rc$ entries with `stride` equal to $c$; a window
created by `fmpz_mat_window_init` points into another matrix and has that
matrix's stride. Code that accesses entries directly must therefore use
`fmpz_mat_entry` or `fmpz_mat_row` rather than assume contiguous storage.
:::
:::

The following conventions apply throughout the module.

- *Indexing.* Rows and columns are indexed from zero. Matrices with zero
  rows or zero columns are allowed everywhere.
- *Dimensions.* Outputs must be initialised with the dimensions of the
  result; no function resizes its output. Most functions check dimensions
  only with `FLINT_ASSERT`.
- *Aliasing.* Unless a function states otherwise, an output may be the
  same matrix as an input only for entrywise operations (addition, scalar
  operations and the like). Partially overlapping windows are never
  allowed.
- *Rational results.* A rational matrix $X$ is returned as a pair
  (`X`, `den`) with $X/\mathit{den}$ the true result. Unless stated
  otherwise the denominator is not guaranteed to be minimal or positive;
  each solver states which denominator it returns.
- *Windows.* All functions accept windows as inputs. As outputs, windows
  are accepted except by functions that replace or free the storage of
  their argument: `fmpz_mat_clear`, `fmpz_mat_swap`, `fmpz_mat_fread`
  and `fmpz_mat_next_col_van_hoeij`.

## Memory management {#memory}

::: function
void fmpz_mat_init(fmpz_mat_t mat, slong rows, slong cols)

Initialises `mat` as a `rows` by `cols` matrix with all entries zero.
Throws an exception if $\mathit{rows} \cdot \mathit{cols}$ overflows an
`slong`.
:::

::: function
void fmpz_mat_init_set(fmpz_mat_t mat, const fmpz_mat_t src)

Initialises `mat` as a copy of `src`. The copy owns its entries even if
`src` is a window.
:::

::: function
void fmpz_mat_clear(fmpz_mat_t mat)

Frees the memory used by `mat`. Must not be called on a window.
:::

::: function
void fmpz_mat_set(fmpz_mat_t mat1, const fmpz_mat_t mat2)

Copies the entries of `mat2` to `mat1`, which must have the same
dimensions. Nothing happens if the two are the same matrix.
:::

::: function
void fmpz_mat_swap(fmpz_mat_t mat1, fmpz_mat_t mat2)

Swaps `mat1` and `mat2` in constant time by swapping their structs. Must
not be used with windows.
:::

::: function
void fmpz_mat_swap_entrywise(fmpz_mat_t mat1, fmpz_mat_t mat2)

Swaps the entries of two matrices of the same dimensions one by one, in
time linear in the number of entries. This works for windows.
:::

## Entries and dimensions {#entries}

::: function
fmpz * fmpz_mat_entry(const fmpz_mat_t mat, slong i, slong j)

Returns a pointer to the entry in row $i$ and column $j$. The pointer can
be passed to any function of [](fmpz) as input or output. No bounds
checking is performed.
:::

::: function
fmpz * fmpz_mat_row(const fmpz_mat_t mat, slong i)

Returns a pointer to the first entry of row $i$. The entries of a row are
contiguous, so the result can be passed to the vector functions of
[](fmpz_vec).
:::

::: function
slong fmpz_mat_nrows(const fmpz_mat_t mat)
slong fmpz_mat_ncols(const fmpz_mat_t mat)

Returns the number of rows, respectively columns, of `mat`.
:::

::: function
int fmpz_mat_is_empty(const fmpz_mat_t mat)
int fmpz_mat_is_square(const fmpz_mat_t mat)

Returns 1 if `mat` has no rows or no columns, respectively if it has as
many rows as columns, and 0 otherwise.
:::

::: function
slong fmpz_mat_max_bits(const fmpz_mat_t mat)

Returns $\pm b$ where $b$ is the largest bit length of the absolute
value of an entry of `mat`, with a negative sign if some entry is
negative. Returns 0 for an empty or zero matrix.
:::

## Windows and concatenation {#windows}

::: function
void fmpz_mat_window_init(fmpz_mat_t window, const fmpz_mat_t mat, slong r1, slong c1, slong r2, slong c2)

Initialises `window` as the submatrix of `mat` formed by rows $r_1$ to
$r_2 - 1$ and columns $c_1$ to $c_2 - 1$, where
$0 \le r_1 \le r_2 \le$ `mat->r` and $0 \le c_1 \le c_2 \le$ `mat->c`.
The window shares its entries with `mat`; no memory is allocated.

::: details
The bounds are checked only with `FLINT_ASSERT`.
:::
:::

::: function
void fmpz_mat_window_clear(fmpz_mat_t window)

Releases a window. This does nothing at present, but every window should
be cleared for future compatibility.
:::

::: function
void fmpz_mat_concat_horizontal(fmpz_mat_t res, const fmpz_mat_t mat1, const fmpz_mat_t mat2)
void fmpz_mat_concat_vertical(fmpz_mat_t res, const fmpz_mat_t mat1, const fmpz_mat_t mat2)

Sets `res` to `mat1` placed to the left of, respectively above, `mat2`.
For an $m \times n$ matrix `mat1`, `mat2` must be $m \times k$ and `res`
$m \times (n + k)$, respectively `mat2` $k \times n$ and `res`
$(m + k) \times n$. Wrappers of `gr_mat_concat_horizontal` and
`gr_mat_concat_vertical`.
:::

## Input and output {#io}

::: function
int fmpz_mat_print(const fmpz_mat_t mat)
int fmpz_mat_fprint(FILE * file, const fmpz_mat_t mat)

Prints `mat` to `stdout`, respectively `file`, in the format read by
`fmpz_mat_fread`: the number of rows and of columns, followed by two
spaces, followed by the entries in row-major order separated by single
spaces. For example, the $2 \times 2$ matrix with rows $(1, 2)$ and
$(3, 4)$ is printed as `2 2  1 2 3 4`. Returns a positive value on
success and a nonpositive value if a write failed.
:::

::: function
int fmpz_mat_print_pretty(const fmpz_mat_t mat)
int fmpz_mat_fprint_pretty(FILE * file, const fmpz_mat_t mat)

Prints `mat` to `stdout`, respectively `file`, with each row in square
brackets followed by a newline, and the whole matrix enclosed in square
brackets: the matrix above is printed as `[[1 2]`, `[3 4]`, `]` on three
lines, with no final newline. Returns a positive value on success and a
nonpositive value if a write failed.
:::

::: function
int fmpz_mat_read(fmpz_mat_t mat)
int fmpz_mat_fread(FILE * file, fmpz_mat_t mat)

Reads a matrix in the format written by `fmpz_mat_fprint` from `stdin`,
respectively `file`. Returns 1 on success and 0 if reading failed.

::: details
If `mat` is $0 \times 0$, it is cleared and reinitialised with the
dimensions read. Otherwise its dimensions must equal the ones read, and an
exception is thrown if they do not, or if a dimension does not fit in an
`slong`. Because `mat` may be reinitialised, it must not be a window.
:::
:::

## Random matrices {#random}

These functions are intended for testing. Several of them produce the
structured matrices that occur in applications of lattice reduction.

::: function
void fmpz_mat_randbits(fmpz_mat_t mat, flint_rand_t state, flint_bitcnt_t bits)

Sets each entry to a random signed integer whose absolute value has
exactly `bits` bits (see `fmpz_randbits`).
:::

::: function
void fmpz_mat_randtest(fmpz_mat_t mat, flint_rand_t state, flint_bitcnt_t bits)
void fmpz_mat_randtest_unsigned(fmpz_mat_t mat, flint_rand_t state, flint_bitcnt_t bits)

Sets the entries to random signed, respectively nonnegative, integers of
at most `bits` bits drawn with `fmpz_randtest`, respectively
`fmpz_randtest_unsigned`, which favour special values. With probability
1/2 `fmpz_mat_randtest` produces a sparse matrix, with a random density.
:::

::: function
void fmpz_mat_randintrel(fmpz_mat_t mat, flint_rand_t state, flint_bitcnt_t bits)

Sets `mat` to a random *integer relations* matrix: column 0 holds random
signed integers of `bits` bits and the remaining columns form an identity
matrix. Requires the number of columns to be one more than the number of
rows, and throws otherwise.
:::

::: function
void fmpz_mat_randsimdioph(fmpz_mat_t mat, flint_rand_t state, flint_bitcnt_t bits, flint_bitcnt_t bits2)

Sets the square matrix `mat` to a random *simultaneous diophantine*
matrix: the first row is $2^{\mathit{bits2}}$ followed by random signed
integers of `bits` bits, and the other rows are $2^{\mathit{bits}}$ times
the corresponding rows of the identity matrix. Throws if `mat` is not
square.
:::

::: function
void fmpz_mat_randntrulike(fmpz_mat_t mat, flint_rand_t state, flint_bitcnt_t bits, ulong q)
void fmpz_mat_randntrulike2(fmpz_mat_t mat, flint_rand_t state, flint_bitcnt_t bits, ulong q)

Sets the square matrix `mat` of even dimension $2d$ to a random
*NTRU-like* matrix $\begin{pmatrix} I & H \\ 0 & qI \end{pmatrix}$,
respectively $\begin{pmatrix} qI & 0 \\ H & I \end{pmatrix}$, where the
$d \times d$ block $H$ is circulant, $H_{ij} = h_{(i + j) \bmod d}$, for
a random vector $h$ of signed integers of `bits` bits. Throws if the
matrix is not square of even dimension.
:::

::: function
void fmpz_mat_randajtai(fmpz_mat_t mat, flint_rand_t state, double alpha)

Sets the square matrix `mat` to a random lower triangular *Ajtai*
matrix. With $b_i = \lfloor (2n - i)^\alpha \rfloor$, the diagonal entry
$(i, i)$ is random in $[1, 2^{b_i - 1}]$, the entries below it in column
$i$ are random in $(-2^{b_i} + 1, 2^{b_i} - 1)$, and the entries above
the diagonal are zero. Throws if `mat` is not square.
:::

::: function
int fmpz_mat_randpermdiag(fmpz_mat_t mat, flint_rand_t state, const fmpz * diag, slong n)

Sets `mat` to a random permutation of the rows and of the columns of the
diagonal matrix whose first $n$ diagonal entries are `diag` (and whose
other entries are zero). Returns the parity of the permutation, 0 if it
is even and 1 if it is odd. Requires $n \le \min(r, c)$. Wrapper of
`gr_mat_randpermdiag`.
:::

::: function
void fmpz_mat_randrank(fmpz_mat_t mat, flint_rand_t state, slong rank, flint_bitcnt_t bits)

Sets `mat` to a random sparse matrix of the given rank, with exactly
`rank` nonzero entries, random integers of at most `bits` bits, in
distinct rows and columns. Throws if `rank` exceeds the number of rows
or columns. `fmpz_mat_randops` turns it into a dense matrix of the same
rank.
:::

::: function
void fmpz_mat_randdet(fmpz_mat_t mat, flint_rand_t state, const fmpz_t det)

Sets the square matrix `mat` to a random sparse matrix with determinant
`det` and as few nonzero entries as possible: `det` is factored and its
prime factors are distributed over $n$ entries in distinct rows and
columns. The matrix is zero if `det` is zero. `fmpz_mat_randops` turns it
into a dense matrix with the same determinant.
:::

::: function
void fmpz_mat_randops(fmpz_mat_t mat, flint_rand_t state, slong count)

Applies `count` random elementary operations to `mat`, each adding or
subtracting a random row to another row, or a random column to another
column. These operations are unimodular, so the rank, the determinant and
the Hermite and Smith normal forms are preserved.
:::

## Comparison and special values {#comparison}

::: function
int fmpz_mat_equal(const fmpz_mat_t mat1, const fmpz_mat_t mat2)

Returns 1 if the matrices have the same dimensions and the same entries,
and 0 otherwise.
:::

::: function
int fmpz_mat_is_zero(const fmpz_mat_t mat)
int fmpz_mat_is_one(const fmpz_mat_t mat)

Returns 1 if `mat` is zero, respectively if it has ones on the main
diagonal and zeros elsewhere (it need not be square), and 0 otherwise.
Empty matrices are both zero and one.
:::

::: function
int fmpz_mat_is_zero_row(const fmpz_mat_t mat, slong i)

Returns 1 if row $i$ of `mat` is zero and 0 otherwise.
:::

::: function
int fmpz_mat_equal_row(fmpz_mat_t M, slong m, slong n)
int fmpz_mat_equal_col(fmpz_mat_t M, slong m, slong n)

Returns 1 if rows, respectively columns, $m$ and $n$ of $M$ are equal,
and 0 otherwise.
:::

::: function
int fmpz_mat_is_diagonal(const fmpz_mat_t A)

Returns 1 if all entries of $A$ outside the main diagonal are zero, and 0
otherwise. $A$ need not be square.
:::

::: function
void fmpz_mat_zero(fmpz_mat_t mat)
void fmpz_mat_one(fmpz_mat_t mat)

Sets `mat` to the zero matrix, respectively to the matrix with ones on
the main diagonal and zeros elsewhere. `mat` need not be square.
:::

## Row and column permutations {#permutations}

In these functions `perm` may be `NULL`. Otherwise it is an array of
length equal to the number of rows (respectively columns), to which the
same permutation is applied, so that it records the operations.

::: function
void fmpz_mat_swap_rows(fmpz_mat_t mat, slong * perm, slong r, slong s)
void fmpz_mat_swap_cols(fmpz_mat_t mat, slong * perm, slong r, slong s)

Swaps rows, respectively columns, $r$ and $s$ of `mat`, and entries $r$
and $s$ of `perm`. Rows are swapped by swapping their entries, so this
costs $O(c)$ even for multiprecision entries.
:::

::: function
void fmpz_mat_invert_rows(fmpz_mat_t mat, slong * perm)
void fmpz_mat_invert_cols(fmpz_mat_t mat, slong * perm)

Reverses the order of the rows, respectively columns, of `mat`, and of the
entries of `perm`.
:::

## Arithmetic {#arithmetic}

::: function
void fmpz_mat_transpose(fmpz_mat_t B, const fmpz_mat_t A)

Sets $B$ to the transpose of $A$; $B$ must have the transposed
dimensions. $B$ may be $A$ if it is square. Wrapper of
`gr_mat_transpose`.
:::

::: function
void fmpz_mat_add(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
void fmpz_mat_sub(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
void fmpz_mat_neg(fmpz_mat_t B, const fmpz_mat_t A)

Sets $C = A + B$, $C = A - B$, respectively $B = -A$, for matrices of the
same dimensions. Aliasing is allowed.
:::

::: function
void fmpz_mat_scalar_mul_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t c)
void fmpz_mat_scalar_mul_si(fmpz_mat_t B, const fmpz_mat_t A, slong c)
void fmpz_mat_scalar_mul_ui(fmpz_mat_t B, const fmpz_mat_t A, ulong c)

Sets $B = cA$. Aliasing is allowed.
:::

::: function
void fmpz_mat_scalar_addmul_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t c)
void fmpz_mat_scalar_addmul_si(fmpz_mat_t B, const fmpz_mat_t A, slong c)
void fmpz_mat_scalar_addmul_ui(fmpz_mat_t B, const fmpz_mat_t A, ulong c)
void fmpz_mat_scalar_submul_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t c)
void fmpz_mat_scalar_submul_si(fmpz_mat_t B, const fmpz_mat_t A, slong c)
void fmpz_mat_scalar_submul_ui(fmpz_mat_t B, const fmpz_mat_t A, ulong c)

Sets $B = B + cA$, respectively $B = B - cA$. Aliasing is allowed.
:::

::: function
void fmpz_mat_scalar_addmul_nmod_mat_fmpz(fmpz_mat_t B, const nmod_mat_t A, const fmpz_t c)
void fmpz_mat_scalar_addmul_nmod_mat_ui(fmpz_mat_t B, const nmod_mat_t A, ulong c)

Sets $B = B + cA$, where the entries of the `nmod_mat_t` $A$ are taken
as integers in $[0, n)$. Used to accumulate $p$-adic digits.
:::

::: function
void fmpz_mat_scalar_divexact_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t c)
void fmpz_mat_scalar_divexact_si(fmpz_mat_t B, const fmpz_mat_t A, slong c)
void fmpz_mat_scalar_divexact_ui(fmpz_mat_t B, const fmpz_mat_t A, ulong c)

Sets $B = A / c$, assuming that $c \ne 0$ divides every entry of $A$;
the result is undefined otherwise. Aliasing is allowed.

::: details
`fmpz_mat_scalar_divexact_fmpz` precomputes the division by $c$ once for
all entries (see `_fmpz_vec_scalar_divexact_fmpz_strided`), which makes it
much faster than entrywise division when $c$ is large.
:::
:::

::: function
void fmpz_mat_scalar_mul_2exp(fmpz_mat_t B, const fmpz_mat_t A, ulong exp)
void fmpz_mat_scalar_tdiv_q_2exp(fmpz_mat_t B, const fmpz_mat_t A, ulong exp)

Sets $B = 2^{\mathit{exp}} A$, respectively $B$ to the entries of
$A / 2^{\mathit{exp}}$ rounded towards zero. Aliasing is allowed.
:::

::: function
void fmpz_mat_scalar_mod_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t m)
void fmpz_mat_scalar_smod(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t P)

Reduces the entries of $A$ modulo $m$ into $[0, |m|)$ (see `fmpz_mod`),
respectively modulo $P > 0$ into the symmetric range $(-P/2, P/2]$, and
writes them to $B$. Aliasing is allowed.
:::

::: function
void fmpz_mat_kronecker_product(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)

Sets $C$ to the Kronecker product $A \otimes B$. For $A$ of size
$m \times n$ and $B$ of size $k \times l$, $C$ must be
$mk \times nl$. $C$ must not be aliased with $A$ or $B$.
:::

::: function
void fmpz_mat_content(fmpz_t ret, const fmpz_mat_t A)

Sets `ret` to the nonnegative greatest common divisor of the entries of
$A$, which is 0 if $A$ is zero or empty. The rows are passed in turn to
`_fmpz_vec_content_chained`, stopping when the GCD becomes 1.
:::

::: function
void fmpz_mat_trace(fmpz_t trace, const fmpz_mat_t mat)

Sets `trace` to the sum of the diagonal entries of the square matrix
`mat` (0 if it is empty). Throws if `mat` is not square.
:::

## Matrix multiplication {#mul}

::: function
void fmpz_mat_mul(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
void fmpz_mat_mul_classical(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
void fmpz_mat_mul_waksman(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
void fmpz_mat_mul_strassen(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
void fmpz_mat_mul_multi_mod(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
void fmpz_mat_mul_fft(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)

Sets $C = AB$ for an $m \times k$ matrix $A$ and a $k \times n$ matrix
$B$; $C$ must be $m \times n$. `fmpz_mat_mul` chooses a method from the
dimensions and the sizes of the entries and allows $C$ to be $A$ or $B$;
the other functions force one method and do not allow aliasing.

::: details
The methods are classical dot products, Waksman's algorithm with about
half the multiplications, Strassen–Winograd recursion, multimodular
multiplication through `nmod_mat_mul`, and the Schönhage–Strassen
transform of the [](fft) module. See [](fmpz_mat/mul) and
[](algorithms/linear-algebra/integer-matrix-multiplication).
:::
:::

::: function
int fmpz_mat_mul_blas(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)

Sets $C = AB$ using floating-point matrix multiplication (`flint_dgemm`),
directly or modulo several primes of about 20 bits, and returns 1, or
returns 0 without computing anything if the method is not applicable.
The result is exact.

::: details
This works on 64-bit machines whether or not FLINT is linked with a BLAS
library, since `flint_dgemm` has its own kernels. It returns 0 on 32-bit
machines, if a dimension exceeds `INT_MAX`, or if $k$ is so large that
no primes are left. See [](fmpz_mat/mul#blas).
:::
:::

::: function
int fmpz_mat_mul_fft_small(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
int fmpz_mat_mul_fft_small_trunc(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B, slong lo_limbs)

Sets $C = AB$ by transforming every entry once with the small-prime FFT
of [](fft_small) and accumulating pointwise products, and returns 1; or
returns 0 if the method is unavailable or not applicable, in which case
$C$ may have been partially overwritten but $A$ and $B$ are unchanged.
The `trunc` variant returns for each entry only the limbs above limb
`lo_limbs` of its absolute value, with its sign, with an error of at most
one unit in the last returned limb.

::: details
Both return 0 if `fft_small` is not available, if a dimension is zero, if
$A$ or $B$ is zero, if $C$ is $A$ or $B$, or if the transform storage
would exceed the budget `flint_fft_small_max_transformed_ring_size` even
with blocking; the `trunc` variant also returns 0 for negative
`lo_limbs` and when the inner dimension would have to be split. The
method pays off only for entries of thousands of bits; see
[](fmpz_mat/mul#fft-small).
:::
:::

::: function
void fmpz_mat_sqr(fmpz_mat_t B, const fmpz_mat_t A)
void fmpz_mat_sqr_bodrato(fmpz_mat_t B, const fmpz_mat_t A)

Sets $B = A^2$ for a square matrix $A$. `fmpz_mat_sqr` allows $B$ to be
$A$; `fmpz_mat_sqr_bodrato` does not.

::: details
`fmpz_mat_sqr_bodrato` uses explicit formulas for dimensions up to 3 and
otherwise Bodrato's squaring variant of Strassen–Winograd
[@Bodrato2010], which needs four squarings and three products of half
size per level. `fmpz_mat_sqr` uses it for
dimension at most 3, and for dimension 4 with entries of at least 1024
bits, as long as the entries have at most 12000 bits; otherwise it calls
`fmpz_mat_mul`.
:::
:::

::: function
void fmpz_mat_pow(fmpz_mat_t B, const fmpz_mat_t A, ulong exp)

Sets $B = A^{\mathit{exp}}$ for a square matrix $A$ by binary
exponentiation, using `fmpz_mat_sqr` and `fmpz_mat_mul`. $A^0$ is the
identity. $B$ may be $A$, but must not be a window.
:::

## Matrix-vector products {#matvec}

::: function
void fmpz_mat_mul_fmpz_vec(fmpz * c, const fmpz_mat_t A, const fmpz * b, slong blen)
void fmpz_mat_mul_fmpz_vec_ptr(fmpz * const * c, const fmpz_mat_t A, const fmpz * const * b, slong blen)

Sets the column vector $c$ of length `A->r` to $Ab$, where $b$ is a
column vector of length `blen`. If `blen` differs from the number of
columns of $A$, the vector is truncated or zero-extended to that length.
The `ptr` variant takes arrays of pointers to the entries. $c$ must not
overlap $b$.
:::

::: function
void fmpz_mat_fmpz_vec_mul(fmpz * c, const fmpz * a, slong alen, const fmpz_mat_t B)
void fmpz_mat_fmpz_vec_mul_ptr(fmpz * const * c, const fmpz * const * a, slong alen, const fmpz_mat_t B)

Sets the row vector $c$ of length `B->c` to $aB$, where $a$ is a row
vector of length `alen`, truncated or zero-extended to the number of rows
of $B$. The `ptr` variant takes arrays of pointers to the entries. $c$
must not overlap $a$.
:::

## Modular reduction and reconstruction {#modular}

::: function
void fmpz_mat_get_nmod_mat(nmod_mat_t Amod, const fmpz_mat_t A)

Sets `Amod` to $A$ reduced modulo the modulus of `Amod`, which must have
the same dimensions.
:::

::: function
void fmpz_mat_set_nmod_mat(fmpz_mat_t A, const nmod_mat_t Amod)
void fmpz_mat_set_nmod_mat_unsigned(fmpz_mat_t A, const nmod_mat_t Amod)

Sets $A$ to the integer matrix congruent to `Amod` modulo its modulus $n$
with entries in $(-n/2, n/2]$, respectively in $[0, n)$.
:::

::: function
void fmpz_mat_CRT_ui(fmpz_mat_t res, const fmpz_mat_t mat1, const fmpz_t m1, const nmod_mat_t mat2, int sign)

Given `mat1` with entries in $[0, m_1)$ (or $(-m_1/2, m_1/2]$ if `sign`
is 1) and `mat2` modulo $m_2$ coprime to $m_1$, sets `res` to the matrix
congruent to `mat1` modulo $m_1$ and to `mat2` modulo $m_2$ with entries
in $[0, m_1 m_2)$ if `sign` is 0 and in $(-m_1 m_2/2, m_1 m_2/2]$ if
`sign` is 1. `res` may be `mat1`. Throws if $m_1$ is not invertible
modulo $m_2$.
:::

::: function
void fmpz_mat_multi_mod_ui(nmod_mat_t * residues, slong nres, const fmpz_mat_t mat)
void fmpz_mat_multi_mod_ui_precomp(nmod_mat_t * residues, slong nres, const fmpz_mat_t mat, const fmpz_comb_t comb)

Reduces `mat` modulo the moduli of the `nres` matrices `residues`, which
must be initialised with the dimensions of `mat` and pairwise coprime
moduli, writing the residues into them. The `precomp` variant uses a
`comb` precomputed with `fmpz_comb_init` for the same moduli in the same
order; `comb` may be `NULL`.

::: details
The reduction uses the subproduct tree of `fmpz_comb_t`, so that reducing
an entry of $b$ bits modulo $k$ word-size moduli costs about
$M(b) \log k$ instead of $kb$ word operations. Rows are distributed over
threads when the number of entries times the number of moduli is large
enough.
:::
:::

::: function
void fmpz_mat_multi_CRT_ui(fmpz_mat_t mat, nmod_mat_t * const residues, slong nres, int sign)
void fmpz_mat_multi_CRT_ui_precomp(fmpz_mat_t mat, nmod_mat_t * const residues, slong nres, const fmpz_comb_t comb, int sign)

Sets `mat` to the matrix congruent to each of the `nres` matrices
`residues` modulo their pairwise coprime moduli, with entries in $[0, M)$
if `sign` is 0 and in $(-M/2, M/2]$ if `sign` is 1, where $M$ is the
product of the moduli. The `precomp` variant uses a precomputed `comb`,
which may be `NULL`. Uses several threads for large matrices.
:::

## Row reduction and echelon forms {#rref}

::: function
slong fmpz_mat_find_pivot_any(const fmpz_mat_t mat, slong start_row, slong end_row, slong c)
slong fmpz_mat_find_pivot_smallest(const fmpz_mat_t mat, slong start_row, slong end_row, slong c)

Returns the index of a row $r$ with `start_row` $\le r <$ `end_row`
whose entry in column $c$ is nonzero, or $-1$ if there is none.
`fmpz_mat_find_pivot_any` returns the first such row,
`fmpz_mat_find_pivot_smallest` the first one whose entry has the smallest
absolute value.
:::

::: function
slong fmpz_mat_fflu(fmpz_mat_t B, fmpz_t den, slong * perm, const fmpz_mat_t A, int rank_check)

Computes a fraction-free LU decomposition of $A$ by Bareiss elimination,
stored in $B$, and returns the rank $r$ of $A$. `den` is set to the last
pivot, a nonzero $r \times r$ minor of $A$ (1 if $r = 0$). Row exchanges
are applied to `perm` if it is not `NULL`. $B$ may be $A$.

::: details
If `rank_check` is nonzero, the function stops as soon as a column without
pivot is found, and then returns 0 with `den` set to 0. See
[](fmpz_mat/rref) for the layout of $B$ and
[](algorithms/linear-algebra/fraction-free).
:::
:::

::: function
slong fmpz_mat_rref(fmpz_mat_t B, fmpz_t den, const fmpz_mat_t A)
slong fmpz_mat_rref_fflu(fmpz_mat_t B, fmpz_t den, const fmpz_mat_t A)
slong fmpz_mat_rref_mul(fmpz_mat_t B, fmpz_t den, const fmpz_mat_t A)

Sets ($B$, `den`) to the reduced row echelon form of $A$, so that
$B/\mathit{den}$ is in reduced row echelon form, and returns the rank of
$A$. The pivot entries of $B$ equal `den`. $B$ may be $A$.

::: details
`fmpz_mat_rref_fflu` uses fraction-free elimination followed by
fraction-free back substitution, and `den` is a nonzero $r \times r$
minor of $A$. `fmpz_mat_rref_mul` computes the echelon structure modulo a
prime, solves for the echelon form with the nonsingular pivot block and
certifies the result exactly; it returns `den` = 1 if $A$ has full column
rank. `fmpz_mat_rref` chooses between them from the dimensions. See
[](fmpz_mat/rref).
:::
:::

::: function
int fmpz_mat_is_in_rref_with_rank(const fmpz_mat_t A, const fmpz_t den, slong rank)

Returns 1 if $A/\mathit{den}$ is in reduced row echelon form with exactly
`rank` nonzero rows, the first `rank` rows, and 0 otherwise. The pivot
entries of $A$ must equal `den`.
:::

::: function
int fmpz_mat_rref_upper_certify_lu_mod_p(fmpz_mat_t E, fmpz_t den, const fmpz_mat_t A, slong rank, const slong * P, const slong * pivs)

Given the echelon structure of the $m \times n$ matrix $A$ modulo a
prime, as computed by `nmod_mat_lu_with_pivots` (the rank $r$ modulo $p$,
the row permutation `P` and the array `pivs` of the $r$ pivot columns
followed by the other columns), tries to compute and certify the first
$r$ rows of the reduced row echelon form of $A$ over $\Q$. Returns 1 on
success, with ($E$, `den`) holding these rows, and 0 if the prime was
unlucky.

::: details
$E$ must be an $r \times n$ zero matrix not aliased with $A$. The
function solves $B E_2 = \mathit{den} \cdot C$ exactly, where $B$ and $C$
are the pivot rows of $A$ restricted to the pivot and the other columns,
checks that $(E, \mathit{den})$ is in reduced row echelon form, and checks
that the other rows of $A$ are annihilated by the nullspace basis that
$E$ defines. Both checks are exact, so a result of 1 is a proof. With
probability close to 1 for a random prime, the rank modulo $p$ equals the
rank of $A$ and the function succeeds; otherwise it returns 0 and $E$ and
`den` are unspecified.
:::
:::

::: function
int fmpz_mat_rank_certify_lu_mod_p(const fmpz_mat_t A, slong rank, const slong * P, const slong * pivs)

Returns 1 if the rank of $A$ over $\Q$ is certified to equal `rank`,
given the echelon structure modulo a prime as for
`fmpz_mat_rref_upper_certify_lu_mod_p`, and 0 if the certification
failed.
:::

::: function
slong fmpz_mat_rank(const fmpz_mat_t A)

Returns the rank of $A$ over $\Q$.

::: details
For entries of at most $(\beta - 4)/2$ bits ($\beta$ = `FLINT_BITS`) the
rank is computed by `fmpz_mat_rank_small_inplace`, integer elimination in
single words that falls back to the general code on overflow; otherwise
by `fmpz_mat_fflu` if $\min(m, n) < 25$ and by `fmpz_mat_rref` (which
then uses the multimodular `fmpz_mat_rref_mul`) for larger matrices. The
result is always exact.
:::
:::

::: function
slong fmpz_mat_rank_small_inplace(fmpz_mat_t B)

Returns the rank of $B$ and destroys $B$. Intended for matrices with
small entries: $B$ is triangularised with integer row operations (a
Euclidean algorithm on the entries of each column, with the smallest
nonzero entry as pivot) in single-word arithmetic, and the remaining
submatrix is passed to `fmpz_mat_fflu` or `fmpz_mat_rref` as soon as an
entry reaches $2^{(\beta - 4)/2}$.
:::

## Strong echelon form and Howell form {#howell}

::: function
void fmpz_mat_strong_echelon_form_mod(fmpz_mat_t A, const fmpz_t mod)

Transforms $A$ so that $A$ modulo `mod` is the strong echelon form of the
input modulo `mod` [@FieHof2014]. $A$ must have at least as many rows as
columns. The entries of the result are reduced into $[0, \mathit{mod})$.
:::

::: function
slong fmpz_mat_howell_form_mod(fmpz_mat_t A, const fmpz_t mod)

Transforms $A$ so that $A$ modulo `mod` is the Howell form of the input
modulo `mod` [@How1986; @StoMul1998], and returns the number of nonzero
rows. $A$ must have at least as many rows as columns.

::: details
The Howell form is computed as the strong echelon form with the zero rows
moved to the bottom; see
[](algorithms/linear-algebra/gaussian-elimination#howell) for both forms.
Two matrices span the same submodule of $(\Z/N\Z)^c$ if and only if their
Howell forms are equal.
:::
:::

## Determinants {#det}

::: function
void fmpz_mat_det(fmpz_t det, const fmpz_mat_t A)
void fmpz_mat_det_cofactor(fmpz_t det, const fmpz_mat_t A)
void fmpz_mat_det_bareiss(fmpz_t det, const fmpz_mat_t A)

Sets `det` to the determinant of the square matrix $A$; the determinant
of an empty matrix is 1. `fmpz_mat_det` throws if $A$ is not square and
chooses a method from the dimension and the entry size;
`fmpz_mat_det_cofactor` uses cofactor expansion and throws for
dimensions above 4; `fmpz_mat_det_bareiss` uses fraction-free
elimination. See [](fmpz_mat/det).
:::

::: function
void fmpz_mat_det_modular(fmpz_t det, const fmpz_mat_t A, int proved)
void fmpz_mat_det_modular_accelerated(fmpz_t det, const fmpz_mat_t A, int proved)
void fmpz_mat_det_modular_given_divisor(fmpz_t det, const fmpz_mat_t A, const fmpz_t d, int proved)

Sets `det` to the determinant of the square matrix $A$, computed modulo
word-size primes and reconstructed by the Chinese remainder theorem.
`fmpz_mat_det_modular_given_divisor` takes a known divisor $d$ of the
determinant and only reconstructs $\det(A)/d$;
`fmpz_mat_det_modular_accelerated` first computes such a divisor with
`fmpz_mat_det_divisor`. If `proved` is 0, the result is only correct with
high probability.

::: details
With `proved` set, primes are used until their product exceeds twice the
Hadamard bound divided by $d$. With `proved` 0, the iteration also stops
once the reconstructed value has not changed while more than 100 bits of
primes were added, which is correct with high probability but not proven.
A zero determinant is detected early and certified with an exact kernel
vector. See [](fmpz_mat/det) and
[](algorithms/linear-algebra/multimodular#det).
:::
:::

::: function
void fmpz_mat_det_bound(fmpz_t bound, const fmpz_mat_t A)
void fmpz_mat_det_bound_submatrix(fmpz_t bound, const fmpz_mat_t A)

Sets `bound` to an upper bound for $|\det(A)|$: the smaller of the
products of the Euclidean norms of the rows and of the columns
(Hadamard's inequality [@Had1893]), rounded up. The `submatrix` variant
omits zero rows and columns from the products, which yields a bound for
the absolute value of every minor of $A$ (of any size), also for
nonsquare $A$. The bound of an empty matrix is 1.
:::

::: function
void fmpz_mat_det_divisor(fmpz_t d, const fmpz_mat_t A)

Sets $d$ to a divisor of the determinant of the square matrix $A$, which
is likely to be its largest invariant factor and hence close to the
determinant, or to 0 if $A$ is singular.

::: details
The function solves $Ax = b$ for $b = (-1, 1, -1, 1, \ldots)^T$ by
$p$-adic lifting (`fmpz_mat_solve_dixon`), reconstructs the entries of
$x$ as fractions ([](algorithms/number-theory/rational-reconstruction))
and returns the least common multiple of their
denominators, which divides the largest invariant factor of $A$ and hence
$\det(A)$ [@AbbottBronsteinMulders1999; @Pan1988].
:::
:::

::: function
int fmpz_mat_is_singular(const fmpz_mat_t A)

Returns 1 if the square matrix $A$ is singular and 0 otherwise. The
result is always proven. Throws if $A$ is not square. The empty matrix is
nonsingular.

::: details
Dimensions up to 4 are decided from the determinant modulo $2^\beta$ (and
the exact determinant if that is zero). Larger matrices are reduced
modulo primes of 17 bits: full rank modulo a prime proves
nonsingularity. A rank drop is confirmed either by computing the
determinant (small matrices) or by an exact kernel vector, obtained by
solving the nonsingular pivot block against a non-pivot column; if that
fails the prime was unlucky and another one is tried, until the product
of the unlucky primes exceeds the Hadamard bound.
:::
:::

## Permanent {#permanent}

::: function
int fmpz_mat_permanent(fmpz_t res, const fmpz_mat_t A)

Sets `res` to the permanent of the square matrix $A$ and returns 1.
Squareness is not checked. Uses Glynn's formula [@Gly2010], with
fixed-size multiword arithmetic when the entries are small and with
threads for dimensions above 10, and otherwise
`gr_mat_permanent_cofactor` or `gr_mat_permanent_glynn`. The cost is
$O(2^n n)$ operations.
:::

## Characteristic and minimal polynomials {#charpoly}

::: function
void fmpz_mat_charpoly(fmpz_poly_t cp, const fmpz_mat_t mat)
void _fmpz_mat_charpoly(fmpz * cp, const fmpz_mat_t mat)
void fmpz_mat_charpoly_berkowitz(fmpz_poly_t cp, const fmpz_mat_t mat)
void _fmpz_mat_charpoly_berkowitz(fmpz * cp, const fmpz_mat_t mat)
void fmpz_mat_charpoly_modular(fmpz_poly_t cp, const fmpz_mat_t mat)
void _fmpz_mat_charpoly_modular(fmpz * cp, const fmpz_mat_t mat)

Sets `cp` to the characteristic polynomial $\det(xI - A)$ of the square
$n \times n$ matrix $A$ = `mat`. The underscore versions write the
$n + 1$ coefficients to `cp`. Only `fmpz_mat_charpoly` checks that the
matrix is square, and throws otherwise.

::: details
`fmpz_mat_charpoly_berkowitz` is a wrapper of `_gr_mat_charpoly_berkowitz`,
division-free with $O(n^4)$ operations on integers.
`fmpz_mat_charpoly_modular` computes the characteristic polynomial modulo
enough 63-bit primes to exceed twice the bound of
`fmpz_mat_charpoly_bound`, in parallel, and combines the results by the
Chinese remainder theorem; the result is proven. `fmpz_mat_charpoly` uses
explicit formulas for $n \le 3$, Berkowitz's algorithm for $n \le 8$ and
when $n$ is smaller than the bit length of the bit size of the entries,
and the modular algorithm otherwise. See [](fmpz_mat/charpoly).
:::
:::

::: function
void fmpz_mat_charpoly_bound(fmpz_t bound, const fmpz_mat_t A)

Sets `bound` to an upper bound for the absolute values of the
coefficients of the characteristic polynomial of the square matrix $A$.

::: details
The coefficient of $x^{n-k}$ is $\pm$ the sum of the $\binom{n}{k}$
principal $k \times k$ minors. By Hadamard's inequality each minor is at
most the product of the $k$ largest row norms of $A$, and also of the $k$
largest column norms. The function returns the smaller of the maxima over
$k$ of $\binom{n}{k}$ times these products.
:::
:::

::: function
void fmpz_mat_minpoly(fmpz_poly_t cp, const fmpz_mat_t mat)
slong _fmpz_mat_minpoly(fmpz * cp, const fmpz_mat_t mat)
void fmpz_mat_minpoly_modular(fmpz_poly_t cp, const fmpz_mat_t mat)
slong _fmpz_mat_minpoly_modular(fmpz * cp, const fmpz_mat_t mat)

Sets `cp` to the minimal polynomial of the square matrix `mat`. The
underscore versions write the coefficients to `cp`, which needs space for
$n + 1$ of them, and return the length. `fmpz_mat_minpoly` throws if
`mat` is not square; the others assume it.

::: details
All four use the multimodular algorithm: minimal polynomials modulo
63-bit primes (`nmod_mat_minpoly_with_gens`) are combined by the Chinese
remainder theorem until the result stabilises, primes giving a smaller
degree are discarded, and the candidate $f$ is accepted when
$f(A) v = 0$ for every Krylov generator $v$ recorded modulo the primes,
which proves that $f$ is the minimal polynomial. A coefficient bound
from [@Dum2007] caps the number of primes. See [](fmpz_mat/charpoly).
:::
:::

::: function
void fmpz_mat_similarity(fmpz_mat_t A, slong r, fmpz_t d)

Replaces the square matrix $A$ by $P^{-1} A P$, where $P$ is the
identity matrix with the entries in row $r$ and every column
$j \notin \{r - 1, r\}$ replaced by $d$. Requires $0 \le r < n$.

::: details
Equivalently, $d$ times column $r$ is added to every column
$j \notin \{r - 1, r\}$, and then $d$ times the sum of the rows
$j \notin \{r - 1, r\}$ is subtracted from row $r$. Such transforms
preserve the determinant and the characteristic and minimal polynomials.
:::
:::

## Solving linear systems {#solve}

The solvers compute a solution of $AX = B$ over $\Q$ as an integer matrix
$X$ and a denominator `den` with $AX = \mathit{den} \cdot B$. The
nonsingular solvers return 1 if $A$ is nonsingular and 0 if it is
singular, which they prove. See [](fmpz_mat/solve) and
[](algorithms/linear-algebra/multimodular#solve).

::: function
int fmpz_mat_solve(fmpz_mat_t X, fmpz_t den, const fmpz_mat_t A, const fmpz_mat_t B)
int fmpz_mat_solve_cramer(fmpz_mat_t X, fmpz_t den, const fmpz_mat_t A, const fmpz_mat_t B)
int fmpz_mat_solve_fflu(fmpz_mat_t X, fmpz_t den, const fmpz_mat_t A, const fmpz_mat_t B)
int fmpz_mat_solve_dixon_den(fmpz_mat_t X, fmpz_t den, const fmpz_mat_t A, const fmpz_mat_t B)
int fmpz_mat_solve_multi_mod_den(fmpz_mat_t X, fmpz_t den, const fmpz_mat_t A, const fmpz_mat_t B)

For a square matrix $A$ of dimension $n$ and an $n \times k$ matrix $B$,
sets ($X$, `den`) to the solution of $AX = B$ and returns 1 if $A$ is
nonsingular; returns 0 if $A$ is singular. `fmpz_mat_solve` chooses the
algorithm; the others force Cramer's rule ($n \le 3$ only), fraction-free
elimination, $p$-adic lifting, respectively multimodular solving.

::: details
`fmpz_mat_solve_cramer` and `fmpz_mat_solve_fflu` return
$\mathit{den} = \det(A)$; `fmpz_mat_solve_dixon_den` and
`fmpz_mat_solve_multi_mod_den` return the smallest positive common
denominator of the solution. $X$ may be $A$ or $B$. On a return value
of 0, $X$ is unspecified and `den` is 0 for the fraction-free solvers.
$X$ may be a window matrix.
:::
:::

::: function
int fmpz_mat_solve_fflu_precomp(fmpz_mat_t X, const slong * perm, const fmpz_mat_t FFLU, const fmpz_mat_t B)

Solves $AX = \mathit{den} \cdot B$ given the output `FFLU` of
`fmpz_mat_fflu` for $A$, with `den` its denominator and `perm` its row
permutation (which must not be `NULL`). Returns 1 on success. For a
nonsingular $A$ this always succeeds; for a rank-deficient $A$ it
computes the solution whose non-pivot entries are zero, if exact
divisions allow it, and returns 0 if a division is inexact.

::: details
The sign of the permutation is not taken into account: if `perm` is odd,
the result satisfies $AX = -\mathit{den} \cdot B$. $X$ may be $B$.
:::
:::

::: function
void fmpz_mat_solve_bound(fmpz_t N, fmpz_t D, const fmpz_mat_t A, const fmpz_mat_t B)

Sets $N$ and $D$ to bounds for the solution of $AX = B$ written over the
denominator $\det(A)$: $D \ge |\det(A)|$ and $N \ge$ the absolute value of
every entry of $\det(A) A^{-1} B$.

::: details
$D$ is the smaller of the row and column Hadamard bounds. $N$ is the
product of the column norms of $A$ with one factor replaced by the largest
column norm of $B$, by Cramer's rule and Hadamard's inequality applied to
columns; the row bound cannot be used for $N$. All norms are accumulated
with upward rounding in `mag_t` arithmetic.
:::
:::

::: function
int fmpz_mat_solve_dixon(fmpz_mat_t X, fmpz_t mod, const fmpz_mat_t A, const fmpz_mat_t B)

Solves $AX = B$ for a nonsingular square $A$ by $p$-adic lifting and
returns 1, or returns 0 if $A$ is singular. The solution is not
reconstructed: $X$ is set to $A^{-1} B$ modulo `mod`, a power of a prime
$p$ exceeding $2 \max(N, D)^2$ with $N$ and $D$ from
`fmpz_mat_solve_bound`, with entries in $[0, \mathit{mod})$, from which
the rational entries can be reconstructed.

::: details
If $A$ or $B$ is empty, the function returns 1 without setting `mod`.
Most applications should call `fmpz_mat_solve_dixon_den`, which
reconstructs the solution
([](algorithms/number-theory/rational-reconstruction#vectors)) and stops
early when it is small.
:::
:::

::: function
int fmpz_mat_can_solve(fmpz_mat_t X, fmpz_t den, const fmpz_mat_t A, const fmpz_mat_t B)
int fmpz_mat_can_solve_fflu(fmpz_mat_t X, fmpz_t den, const fmpz_mat_t A, const fmpz_mat_t B)
int fmpz_mat_can_solve_dixon_den(fmpz_mat_t X, fmpz_t den, const fmpz_mat_t A, const fmpz_mat_t B)
int fmpz_mat_can_solve_multi_mod_den(fmpz_mat_t X, fmpz_t den, const fmpz_mat_t A, const fmpz_mat_t B)

For an $m \times n$ matrix $A$ of any rank and an $m \times k$ matrix
$B$, returns 1 if $AX = B$ has a solution over $\Q$, with ($X$, `den`)
set to one such solution ($X$ is $n \times k$), and 0 if the system is
inconsistent. The answer is proven.

::: details
The solution returned has zeros in the rows corresponding to the
non-pivot columns of $A$. `fmpz_mat_can_solve_fflu` uses fraction-free
elimination; the other two find the rank structure modulo primes, solve
the nonsingular pivot block with `fmpz_mat_solve_dixon_den`,
respectively `fmpz_mat_solve_multi_mod_den`, and verify the other
equations exactly. `fmpz_mat_can_solve` uses `fmpz_mat_can_solve_fflu`
for $m \le 15$ and otherwise the modular method, solving the pivot block
with `fmpz_mat_solve`. On a return value of 0, $X$ and `den` are
unspecified. $X$ must not be aliased with $A$ or $B$. See
[](fmpz_mat/solve#can-solve).
:::
:::

## Nullspace and inverse {#nullspace}

::: function
slong fmpz_mat_nullspace(fmpz_mat_t res, const fmpz_mat_t mat)

Sets the first $n - r$ columns of `res` to a basis of the right nullspace
of the $m \times n$ matrix `mat` over $\Q$, where $r$ is the rank, and
returns $n - r$. The other columns of `res` are set to zero. `res` must
be $n \times n$.

::: details
The basis is read off from `fmpz_mat_rref`: for each non-pivot column $f$
the basis vector has $-\mathit{den}$ in position $f$ and the entries of
column $f$ of the scaled echelon form in the pivot positions. The vectors
are integral but in general neither primitive nor a basis of the integer
kernel lattice $\{x \in \Z^n : Ax = 0\}$; for that, compute the Hermite
normal form with transformation of $A^T$ or use LLL.
:::
:::

::: function
int fmpz_mat_inv(fmpz_mat_t B, fmpz_t den, const fmpz_mat_t A)

Sets ($B$, `den`) to the inverse of the square matrix $A$, so that
$AB = BA = \mathit{den} \cdot I$, and returns 1 if $A$ is nonsingular;
returns 0 if $A$ is singular. $B$ may be $A$.

::: details
Dimensions up to 2 use explicit formulas with $\mathit{den} = \det(A)$;
larger matrices are solved against the identity matrix with
`fmpz_mat_solve`, so that `den` is $\det(A)$ for dimensions up to 15 and
the minimal positive denominator (the largest invariant factor of $A$)
above.
:::
:::

## Hermite normal form {#hnf}

FLINT computes the *row-style* Hermite normal form: $H = UA$ with $U$
unimodular, $H$ in row echelon form with its zero rows last, positive
pivots, and the entries above each pivot reduced into $[0, \text{pivot})$.
The rows of $H$ are a basis of the lattice spanned by the rows of $A$.
See [](fmpz_mat/hnf) and [](algorithms/linear-algebra/hermite-smith).

::: function
void fmpz_mat_hnf(fmpz_mat_t H, const fmpz_mat_t A)
void fmpz_mat_hnf_classical(fmpz_mat_t H, const fmpz_mat_t A)
void fmpz_mat_hnf_xgcd(fmpz_mat_t H, const fmpz_mat_t A)
void fmpz_mat_hnf_minors(fmpz_mat_t H, const fmpz_mat_t A)
void fmpz_mat_hnf_pernet_stein(fmpz_mat_t H, const fmpz_mat_t A, flint_rand_t state)

Sets $H$, of the same dimensions as $A$, to the Hermite normal form of
$A$. $H$ may be $A$. `fmpz_mat_hnf` chooses the algorithm; the others
force the classical algorithm, elimination with extended gcds, the
algorithm of Kannan and Bachem, respectively that of Pernet and Stein.

::: details
`fmpz_mat_hnf_minors` requires $A$ to have full column rank; the others
accept any matrix. `fmpz_mat_hnf_pernet_stein` uses randomness from
`state`, but its result does not depend on it since the Hermite normal
form is unique.
:::
:::

::: function
void fmpz_mat_hnf_modular(fmpz_mat_t H, const fmpz_mat_t A, const fmpz_t D)
void fmpz_mat_hnf_modular_eldiv(fmpz_mat_t A, const fmpz_t D)

Computes the Hermite normal form of an $m \times n$ matrix $A$ with
$m \ge n$ and full column rank by computing modulo $D$.
`fmpz_mat_hnf_modular` requires $D$ to be a positive multiple of the
determinant of the lattice spanned by the rows of $A$;
`fmpz_mat_hnf_modular_eldiv` works in place and only requires $D$ to be
a positive multiple of the largest elementary divisor of $A$. Neither
condition is checked.

::: details
`fmpz_mat_hnf_modular` is the algorithm of Domich, Kannan and Trotter
[@DomKanTro1987]; `fmpz_mat_hnf_modular_eldiv` computes the strong
echelon form modulo $D$ (with `nmod_mat_strong_echelon_form` when $D$
fits in a word) and replaces zero diagonal entries by $D$. See
[](algorithms/linear-algebra/hermite-smith#modular).
:::
:::

::: function
void fmpz_mat_hnf_transform(fmpz_mat_t H, fmpz_mat_t U, const fmpz_mat_t A)
void fmpz_mat_hnf_minors_transform(fmpz_mat_t H, fmpz_mat_t U, const fmpz_mat_t A)

Sets $H$ to the Hermite normal form of the $m \times n$ matrix $A$ and
$U$ to an $m \times m$ unimodular matrix with $UA = H$.
`fmpz_mat_hnf_minors_transform` requires full column rank.

::: details
`fmpz_mat_hnf_transform` uses `fmpz_mat_hnf_minors_transform` if
$m \ge n$ and $A$ has full column rank modulo a random 59-bit prime, and
otherwise computes the Hermite normal form of $(A \mid I_m)$, whose right
block is then $U$. The transformation matrix is not unique when $A$ does
not have full row rank.
:::
:::

::: function
int fmpz_mat_is_in_hnf(const fmpz_mat_t A)

Returns 1 if $A$ is in row-style Hermite normal form and 0 otherwise.
:::

## Smith normal form {#snf}

The Smith normal form of an $m \times n$ matrix $A$ is the diagonal
matrix $S = UAV$, with $U$ and $V$ unimodular, whose diagonal entries
$d_1 \mid d_2 \mid \cdots \mid d_r$, followed by zeros, are positive;
$d_1, \ldots, d_r$ are the elementary divisors (invariant factors) of
$A$. See [](fmpz_mat/snf) and
[](algorithms/linear-algebra/hermite-smith#smith).

::: function
void fmpz_mat_snf(fmpz_mat_t S, const fmpz_mat_t A)
void fmpz_mat_snf_kannan_bachem(fmpz_mat_t S, const fmpz_mat_t A)
void fmpz_mat_snf_iliopoulos(fmpz_mat_t S, const fmpz_mat_t A, const fmpz_t mod)

Sets $S$, of the same dimensions as $A$, to the Smith normal form of $A$.
$S$ may be $A$. `fmpz_mat_snf` and `fmpz_mat_snf_kannan_bachem` accept
any matrix; `fmpz_mat_snf_iliopoulos` works modulo `mod` and requires $A$
to be square and nonsingular with `mod` a positive multiple of its
largest elementary divisor, for example of $|\det(A)|$.

::: details
`fmpz_mat_snf` returns at once for a diagonal matrix
(`fmpz_mat_snf_diagonal`) and otherwise alternates row and column Hermite
normal forms (`fmpz_mat_hnf` of $X$ and of $X^T$) until the matrix is
diagonal, then fixes the divisibility of the diagonal with gcd and lcm
operations.
:::
:::

::: function
void fmpz_mat_snf_diagonal(fmpz_mat_t S, const fmpz_mat_t A)

Sets $S$ to the Smith normal form of the diagonal matrix $A$, by
replacing pairs of diagonal entries by their gcd and lcm. $A$ need not be
square; that it is diagonal is not checked, and only its diagonal is
used. $S$ may be $A$.
:::

::: function
void fmpz_mat_snf_transform(fmpz_mat_t S, fmpz_mat_t U, fmpz_mat_t V, const fmpz_mat_t A)

Sets $S$ to the Smith normal form of the $m \times n$ matrix $A$ and $U$
($m \times m$) and $V$ ($n \times n$) to unimodular matrices with
$S = UAV$. $U$ or $V$ may be `NULL` if it is not wanted. $S$ may be $A$;
$U$ and $V$ must be distinct from $A$, $S$ and each other.
:::

::: function
slong fmpz_mat_elementary_divisors(fmpz * ed, const fmpz_mat_t A)

Sets `ed` to the nonzero elementary divisors $d_1 \mid \cdots \mid d_r$
of $A$, where $r$ is the rank, and returns $r$. `ed` must have room for
$\min(m, n)$ entries.

::: details
The function computes the Hermite normal form of $A$, factors its pivots
and, if all their prime factors fit in a word, determines the $p$-adic
valuations of the elementary divisors for each such prime $p$ with
Lübeck's algorithm [@Lüb2002], which needs only ranks and nullspaces
modulo $p$. If a pivot has more than $2\beta$ bits or cannot be factored
cheaply, it falls back to `fmpz_mat_snf`.
:::
:::

::: function
int fmpz_mat_is_in_snf(const fmpz_mat_t A)

Returns 1 if $A$ is in Smith normal form and 0 otherwise.
:::

## Gram matrices and floating-point conversion {#gram}

::: function
void fmpz_mat_gram(fmpz_mat_t B, const fmpz_mat_t A)

Sets $B$ to the Gram matrix $AA^T$ of the rows of the $m \times n$
matrix $A$; $B$ must be $m \times m$, and may be $A$ if $A$ is square.
Throws if the dimensions are wrong.
:::

::: function
int fmpz_mat_is_spd(const fmpz_mat_t A)

Returns 1 if $A$ is symmetric positive definite and 0 otherwise. The
answer is proven.

::: details
After checking symmetry, the function tries an $LDL^T$ decomposition in
ball arithmetic (`arb_mat_ldl`) at precisions doubling from 32 bits up to
four times the size of the entries; if that does not certify positive
definiteness, it decides exactly from the signs of the coefficients of
the characteristic polynomial (Descartes' rule of signs).
:::
:::

::: function
int fmpz_mat_get_d_mat(d_mat_t B, const fmpz_mat_t A)
int fmpz_mat_get_d_mat_transpose(d_mat_t B, const fmpz_mat_t A)

Sets $B$ to $A$, respectively to its transpose, with each entry rounded
to a `double`, and returns 0. Returns $-1$ if some entry exceeds
`DBL_MAX` in absolute value, in which case $B$ is partially written.
:::

::: function
void fmpz_mat_chol_d(d_mat_t R, const fmpz_mat_t A)

Sets $R$ to the lower triangular Cholesky factor of the symmetric
positive definite matrix $A$, with $A = RR^T$, computed in `double`
arithmetic. Only the lower triangle of $R$ is written. Throws if the
dimensions are incompatible; positive definiteness is not checked.
:::

## LLL reduction {#lll}

These functions check and compute LLL-reduced bases with exact rational
arithmetic. They serve as references; for actual computations use
`fmpz_lll`. See [](algorithms/linear-algebra/lll).

::: function
int fmpz_mat_is_reduced(const fmpz_mat_t A, double delta, double eta)
int fmpz_mat_is_reduced_gram(const fmpz_mat_t A, double delta, double eta)

Returns 1 if the rows of $A$ are $(\delta, \eta)$-LLL-reduced, and 0
otherwise. `fmpz_mat_is_reduced_gram` takes the Gram matrix of the basis
instead of the basis. The answer is proven.

::: details
`fmpz_mat_is_reduced` calls `gr_mat_is_row_lll_reduced_naive` in real
ball arithmetic at doubling precisions from 64 bits and, once the
precision exceeds about $(64 + b) \max(m, n)$ for $b$-bit entries, in
exact rational arithmetic. `fmpz_mat_is_reduced_gram` computes the
Gram–Schmidt data over $\Q$. A basis with linearly dependent rows is not
reduced.
:::
:::

::: function
int fmpz_mat_is_reduced_with_removal(const fmpz_mat_t A, double delta, double eta, const fmpz_t gs_B, int newd)
int fmpz_mat_is_reduced_gram_with_removal(const fmpz_mat_t A, double delta, double eta, const fmpz_t gs_B, int newd)

Returns 1 if the first `newd` rows of $A$ are $(\delta, \eta)$-LLL-reduced
and the squared Gram–Schmidt norms of the other rows are at least
`gs_B`, and 0 otherwise. The `gram` variant takes the Gram matrix. These
define the output condition of LLL with removals.
:::

::: function
void fmpz_mat_lll_original(fmpz_mat_t A, const fmpq_t delta, const fmpq_t eta)

Reduces the rows of $A$ in place with the original LLL algorithm
[@LLL1982] in exact rational arithmetic, with parameters $\delta$ and
$\eta$. The rows must be linearly independent.
:::

::: function
void fmpz_mat_lll_storjohann(fmpz_mat_t A, const fmpq_t delta, const fmpq_t eta)

Reduces the rows of $A$ in place with Storjohann's modified LLL
[@Sto1996], which works with the integral Gram–Schmidt data obtained by
fraction-free elimination of the Gram matrix and has a better complexity
in the dimension. The rows must be linearly independent.
:::

## Special matrices {#special}

::: function
int fmpz_mat_is_hadamard(const fmpz_mat_t A)

Returns 1 if $A$ is a Hadamard matrix (a square matrix with entries
$\pm 1$ and $AA^T = nI$) and 0 otherwise.
:::

::: function
int fmpz_mat_hadamard(fmpz_mat_t A)

Sets the square matrix $A$ of dimension $n$ to a Hadamard matrix and
returns 1 if FLINT can construct one, and returns 0 otherwise. The
constructions are Sylvester's for $n = 2^v$ and Paley's for
$n = 2^v(q + 1)$ and $n = 2^v \cdot 2(q + 1)$ with $q$ a prime power
[@Syl1867; @Pal1933].
:::

## Helpers for polynomial factoring {#factoring}

::: function
int fmpz_mat_col_partition(slong * part, fmpz_mat_t M, int short_circuit)

Returns the number $p$ of distinct columns of $M$ and sets `part`, of
length the number of columns, so that `part[i]` and `part[j]` are equal
elements of $[0, p)$ exactly when columns $i$ and $j$ are equal. If
`short_circuit` is nonzero, returns 0 as soon as $p$ is known to exceed
the number of rows. Used in van Hoeij's factoring algorithm.
:::

::: function
int fmpz_mat_next_col_van_hoeij(fmpz_mat_t M, fmpz_t P, fmpz_mat_t col, slong exp, slong U_exp)

One step of the knapsack lattice construction of van Hoeij's factoring
algorithm [@vHN2010]. Given the lattice $M$ whose leading block $U$ is
the combinatorial part scaled by $2^{\mathit{U\_exp}}$, a new column
`col` of data modulo $P$ and its bound $2^{\mathit{exp}}$, decides
whether the data carries enough information to be worth an LLL call. If
not, returns 0; otherwise appends the truncated column
$U \cdot \mathit{col}$ and a row with the scaled modulus to $M$ and
returns the scaling of the new column. $M$ is reallocated, so it must
not be a window.
:::

## Internal functions {#internal}

::: variable
int flint_fmpz_mat_force_small_primes

A testing hook, not part of the public API. When nonzero, the multimodular
and $p$-adic algorithms for solving, determinants, rank and echelon forms
start their prime search at 2 instead of at word-size primes, so that bad
primes actually occur. Matrix multiplication is not affected.
:::

::: macro
FMPZ_MAT_SOLVE_MIN_WORK_PER_THREAD
FMPZ_MAT_CRT_MIN_WORK_PER_THREAD

The minimum amount of work (word operations, respectively entries times
primes) for which the solvers, respectively the multimodular reduction
and reconstruction, spawn a thread.
:::

::: function
void _fmpz_mat_window_readonly_init_strip_initial_zero_rows(fmpz_mat_t A, const fmpz_mat_t B)

Initialises $A$ as a window of $B$ without its leading zero rows. Release
with `_fmpz_mat_window_readonly_clear`, a macro that does nothing.
:::

::: function
void _fmpz_mat_mul_multi_mod(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B, int sign, flint_bitcnt_t Cbits)
int _fmpz_mat_mul_blas(fmpz_mat_t C, const fmpz_mat_t A, flint_bitcnt_t Abits, const fmpz_mat_t B, flint_bitcnt_t Bbits, int sign, flint_bitcnt_t Cbits)
void _fmpz_mat_mul_fft(fmpz_mat_t C, const fmpz_mat_t A, slong abits, const fmpz_mat_t B, slong bbits, int sign)

Versions of `fmpz_mat_mul_multi_mod`, `fmpz_mat_mul_blas` and
`fmpz_mat_mul_fft` that take precomputed size information:
$|A_{ij}| < 2^{\mathit{Abits}}$, $|B_{ij}| < 2^{\mathit{Bbits}}$,
$|C_{ij}| < 2^{\mathit{Cbits}}$, and `sign` = 1 if an entry of $A$ or $B$
may be negative, 0 otherwise. The bounds must be correct.
:::

::: function
void _fmpz_mat_mul_small_1(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
void _fmpz_mat_mul_small_2a(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
void _fmpz_mat_mul_small_2b(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)

Classical multiplication for matrices whose entries are small `fmpz`
(at most `SMALL_FMPZ_BITCOUNT_MAX` bits), accumulating each dot product
in one signed word (the result entries must be small `fmpz`), in two
words (the sums must fit in a signed two-word integer), or in three
words.
:::

::: function
void _fmpz_mat_mul_small(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
void _fmpz_mat_mul_small_internal(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B, flint_bitcnt_t Cbits)

Multiplies matrices whose entries are small `fmpz` with blocked,
multithreaded dot products accumulated in one, two or three words,
according to the bound $|C_{ij}| < 2^{\mathit{Cbits}}$. The first
function computes the bound itself.
:::

::: function
void _fmpz_mat_mul_double_word(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
void _fmpz_mat_mul_double_word_internal(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B, int sign, flint_bitcnt_t bits)

Multiplies matrices whose entries fit in two words (with one more bit
when signed) by packing them into two-word integers and accumulating the
products in four or five words, using several threads. The internal
version takes the sign flag and the bound $|C_{ij}| < 2^{\mathit{bits}}$.
:::

::: function
void fmpz_mat_multi_mod_2_ui_precomp(nmod_mat_t * residues_A, nmod_mat_t * residues_B, slong nres, const fmpz_mat_t A, const fmpz_mat_t B, const fmpz_comb_t comb)

Reduces both $A$ and $B$ modulo the same `nres` moduli in one pass, as
`fmpz_mat_multi_mod_ui_precomp`. $B$ and `residues_B` may be `NULL`.
:::

::: function
ulong fmpz_mat_find_good_prime_and_invert(nmod_mat_t Ainv, const fmpz_mat_t A, const fmpz_t det_bound)

Searches primes $p$ above $2^{59}$ for one modulo which the square
matrix $A$ is nonsingular, sets `Ainv` to $A^{-1} \bmod p$ and returns
$p$. Returns 0 if $A$ is proven singular, by an exact kernel vector or by
the product of the primes tested exceeding `det_bound`, an upper bound
for $|\det(A)|$.
:::

::: function
ulong _fmpz_mat_find_good_prime_and_lu(nmod_mat_t LU, slong * P, const fmpz_mat_t A, const fmpz_t det_bound)
ulong _fmpz_mat_find_good_prime_and_lu2(nmod_mat_t LU, slong * P, const fmpz_mat_t A, const fmpz_t det_bound, slong bits, int A_is_known_nonsingular)

As `fmpz_mat_find_good_prime_and_invert`, but returns the LU
decomposition of $A$ modulo $p$ in `LU` and `P` as computed by
`nmod_mat_lu_with_pivots`. The second version starts the search at
$2^{\mathit{bits}}$; if `A_is_known_nonsingular` is set it attempts no
certification and `det_bound` may be `NULL`.
:::

::: function
int _fmpz_mat_certify_singular_lu_mod_p(const fmpz_mat_t A, slong rank, const slong * P, const slong * pivs, ulong p)

Given the rank $r < n$ of the $m \times n$ matrix $A$ modulo the prime
$p$ and its echelon structure (`P` and `pivs` as from
`nmod_mat_lu_with_pivots`), tries to prove that $A$ has rank less than
$n$ by exhibiting a nonzero integer kernel vector, obtained by solving
the pivot block against the first non-pivot column. Returns 1 if this
succeeds and 0 otherwise; a return value of 1 is a proof whatever $p$.
:::

::: function
ulong * fmpz_mat_dixon_get_crt_primes(slong * num_primes, const fmpz_mat_t A, ulong p)

Returns an array of `num_primes` primes, $p$ and the primes following
it, whose product exceeds $2n \lVert A \rVert_\infty (p - 1)$, so that
products $Ay$ with $y$ reduced modulo $p$ can be computed
multimodularly. The caller frees the array with `flint_free`.
:::

::: function
void _fmpz_mat_solve_dixon(fmpz_mat_t X, fmpz_t mod, const fmpz_mat_t A, const fmpz_mat_t B, const nmod_mat_t Ainv, ulong p, const fmpz_t N, const fmpz_t D)

The lifting loop of `fmpz_mat_solve_dixon`, given $A^{-1} \bmod p$ and the
bounds $N$, $D$ of `fmpz_mat_solve_bound`.
:::

::: function
int _fmpz_mat_solve_dixon_den_given_lu(fmpz_mat_t X, fmpz_t den, const fmpz_mat_t A, const fmpz_mat_t B, const nmod_mat_t LU, const slong * P, ulong p)

Solves $AX = \mathit{den} \cdot B$ like `fmpz_mat_solve_dixon_den`, given
a prime $p$ modulo which $A$ is nonsingular and the LU decomposition
(`LU`, `P`) of $A$ modulo $p$ as produced by `nmod_mat_lu`. Performs no
prime search and no singularity certification. Returns 1 on success.
:::

::: function
int _fmpz_mat_reconstruct_matwise(fmpz_mat_t Z, fmpz_t den, const fmpz_mat_t x, const fmpz_t m, const fmpz_t N, const fmpz_t D, const fmpz_t den_init)

Tries to reconstruct a matrix $Z/\mathit{den}$ congruent to $x$ modulo
$m$, with entries of $Z$ bounded by $N$ in absolute value and a common
denominator at most $D$ that is a multiple of `den_init`. Each entry is
first tried with the denominator found so far, and otherwise rationally
reconstructed with numerator bound $N$ and denominator bound
$D/\mathit{den}$, the new factor being multiplied into `den` and the
earlier entries. Returns 1 on success and 0 if an entry could not be
reconstructed. The result is unique when $2ND < m$.
:::

::: function
int _fmpz_mat_solve_probe_reconstruct(fmpz_t num, fmpz_t den, const fmpz_t s, const fmpz_t m, const fmpz_t den_prev, int have_prev)

Tries to identify the rational number $\mathit{num}/\mathit{den}$
congruent to $s$ modulo $m$, where $s$ is a random linear combination of
the entries of a solution known modulo $m$: first with the denominators 1
and `den_prev` (if `have_prev` is set), then from the continued fraction
of $s/m$, accepting only results with a margin of 30 bits. Returns 1 on
success. A wrong result is possible with probability about $2^{-30}$; the
callers never rely on it for correctness.
:::

::: function
int _fmpz_mat_solve_attempt_worthwhile(const fmpz_t m, const fmpz_t snum, const fmpz_t sden, const fmpz_t Amax, const fmpz_t Bmax, slong n, slong cols, slong step_bits, double step_cost)

Decides from a probe result (`snum`/`sden`) whether reconstructing the
full solution modulo $m$ can succeed: either the size certificate can
hold at this modulus, or an explicit verification is estimated to be
cheaper than the lifting steps still needed to reach it.
:::

::: function
int _fmpz_mat_solve_reconstruct_attempt(fmpz_mat_t Z, fmpz_t den, const fmpz_mat_t x, const fmpz_t m, const fmpz_t N, const fmpz_t D, const fmpz_t sden, const fmpz_t Amax, const fmpz_t Bmax, const fmpz_mat_t A, const fmpz_mat_t B, slong step_bits, double step_cost, slong den_margin_bits)

Tries to obtain a certified solution $Z/\mathit{den}$ of
$AZ = \mathit{den} \cdot B$ from the solution $x$ modulo $m$, using the
candidate denominator `sden` from a probe with a margin of
`den_margin_bits` bits, `Amax` $= n \lVert A \rVert_\infty$ and
`Bmax` $= \lVert B \rVert_\infty$. Returns 1 if the result is certified,
either by $n \lVert A \rVert \lVert Z \rVert + \mathit{den} \lVert B
\rVert < m$ or by an explicit multiplication, and 0 otherwise. See
[](algorithms/linear-algebra/multimodular#termination).
:::

::: function
int _fmpz_mat_can_solve_auto_den(fmpz_mat_t X, fmpz_t den, const fmpz_mat_t A, const fmpz_mat_t B)

The modular method of `fmpz_mat_can_solve`, with the pivot block solved
by `fmpz_mat_solve`.
:::

::: function
slong _fmpz_mat_snf_iter_bound(const fmpz_mat_t A)

Returns the number of alternating row and column Hermite normal forms
after which `fmpz_mat_snf` and `fmpz_mat_snf_transform` throw, as a guard
against an infinite loop: $16 + 4 M^2 (b + \lfloor \log_2 M \rfloor + 2)$
for $M = \max(m, n)$ and $b$-bit entries, after the worst case bound of
Kannan and Bachem [@KanBac1979].
:::
