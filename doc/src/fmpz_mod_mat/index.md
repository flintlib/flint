---
title: fmpz_mod_mat.h
header: fmpz_mod_mat.h
description: Dense matrices over $\Z/n\Z$ for a modulus of any size, with entries stored as `fmpz` residues, multiplication by exact integer products, and Gaussian elimination through the generic matrix algorithms.
extra-headers: fmpz_mod_poly
category: Matrices
order: 5
---

This module implements dense matrices over $\Z/n\Z$ for an arbitrary
modulus $n \ge 1$ given by an `fmpz_mod_ctx_t`. A matrix is an
`fmpz_mat_t` whose entries are canonical residues in $[0, n)$; it does
not store the modulus, so every function takes the context. Many
functions are thin wrappers: memory management, windows, comparison and
printing are those of [`fmpz_mat`](fmpz_mat); multiplication computes the exact
integer product with `fmpz_mat_mul` (multimodular and BLAS-based
methods) and reduces each entry once; and the Gaussian elimination
functions (LU decomposition, solving, rank, reduced row echelon form,
nullspace, determinant, characteristic and minimal polynomials) call the
generic [`gr_mat`](gr_mat) algorithms over the `fmpz_mod` ring.

Most elimination functions require $n$ to be prime, which is not
checked; the guide [](fmpz_mod_mat/composite) explains what happens for
composite moduli and which functions work for every modulus (the
determinant, the characteristic polynomial, Howell and strong echelon
forms). For moduli of 2 to 16 limbs, generic `gr_mat` code over the
`mpn_mod` ring is usually faster ([](mpn_mod/mat-mul)).

::: summary algorithms/linear-algebra/gaussian-elimination :::

## Types and conventions {#types}

::: type
fmpz_mod_mat_struct
fmpz_mod_mat_t

`fmpz_mod_mat_struct` is the same type as `fmpz_mat_struct`, and
`fmpz_mod_mat_t` an array of one of them, passed by reference. The
fields `entries`, `r`, `c` and `stride` are those of `fmpz_mat`, so an
`fmpz_mod_mat_t` can be passed to `fmpz_mat` functions (and to `gr_mat`
functions over a `gr_ctx_init_fmpz_mod` ring) by a cast.
:::

- Entries must be canonical, in $[0, n)$, unless stated otherwise, and
  results are canonical. Entries can be modified directly through
  `fmpz_mod_mat_entry`, but the caller must keep them canonical.
- Dimensions must be compatible; this is not checked unless stated.
- Unless stated otherwise, the output may alias the inputs.

## Element access {#entries}

::: function
fmpz * fmpz_mod_mat_entry(const fmpz_mod_mat_t mat, slong i, slong j)
fmpz * fmpz_mod_mat_row(const fmpz_mod_mat_t mat, slong i)

Return a pointer to the entry in row $i$ and column $j$, respectively to
the first entry of row $i$ (the entries of a row are contiguous).
Indices are 0-based and not checked.
:::

::: function
void fmpz_mod_mat_set_entry(fmpz_mod_mat_t mat, slong i, slong j, const fmpz_t val, const fmpz_mod_ctx_t ctx)
void fmpz_mod_mat_get_entry(fmpz_t x, const fmpz_mod_mat_t mat, slong i, slong j, const fmpz_mod_ctx_t ctx)

Copy `val` into the entry in row $i$ and column $j$ (without reduction:
`val` must be canonical), respectively copy that entry to $x$.
:::

## Memory management {#memory}

::: function
void fmpz_mod_mat_init(fmpz_mod_mat_t mat, slong rows, slong cols, const fmpz_mod_ctx_t ctx)

Initialises `mat` as a `rows` by `cols` zero matrix.
:::

::: function
void fmpz_mod_mat_init_set(fmpz_mod_mat_t mat, const fmpz_mod_mat_t src, const fmpz_mod_ctx_t ctx)

Initialises `mat` as a copy of `src`, with the same dimensions.
:::

::: function
void fmpz_mod_mat_clear(fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)

Frees the memory used by `mat`.
:::

## Basic properties and manipulation {#basic}

::: function
slong fmpz_mod_mat_nrows(const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)
slong fmpz_mod_mat_ncols(const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)

Return the number of rows and of columns of `mat`.
:::

::: function
int fmpz_mod_mat_is_empty(const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)
int fmpz_mod_mat_is_square(const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)

Return 1 if `mat` has zero rows or zero columns, respectively as many
rows as columns, and 0 otherwise.
:::

::: function
void fmpz_mod_mat_zero(fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)
void fmpz_mod_mat_one(fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)

Set `mat` to the zero matrix, respectively to the matrix with ones on the
main diagonal and zeros elsewhere (the zero matrix if $n = 1$). `mat`
need not be square.
:::

::: function
void fmpz_mod_mat_set(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, const fmpz_mod_ctx_t ctx)

Sets $B$ to a copy of $A$, which must have the same dimensions.
:::

::: function
void fmpz_mod_mat_swap(fmpz_mod_mat_t mat1, fmpz_mod_mat_t mat2, const fmpz_mod_ctx_t ctx)
void fmpz_mod_mat_swap_entrywise(fmpz_mod_mat_t mat1, fmpz_mod_mat_t mat2, const fmpz_mod_ctx_t ctx)

Swap two matrices. `fmpz_mod_mat_swap` swaps the structures in constant
time; `fmpz_mod_mat_swap_entrywise` swaps the entries one by one, which
also works for windows, and requires equal dimensions.
:::

::: function
void fmpz_mod_mat_transpose(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, const fmpz_mod_ctx_t ctx)

Sets $B$ to the transpose of $A$. Aliasing is allowed for square
matrices.
:::

## Conversions {#conversions}

::: function
void fmpz_mod_mat_set_fmpz_mat(fmpz_mod_mat_t A, const fmpz_mat_t B, const fmpz_mod_ctx_t ctx)

Sets $A$ to the integer matrix $B$ with every entry reduced modulo $n$.
$B$ may have any entries and the same dimensions as $A$. Large matrices
are reduced in parallel.
:::

::: function
void fmpz_mod_mat_get_fmpz_mat(fmpz_mat_t A, const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx)

Sets $A$ to the integer matrix with the entries of $B$ (in $[0, n)$).
:::

::: function
void fmpz_mod_mat_set_nmod_mat(fmpz_mod_mat_t A, const nmod_mat_t B, const fmpz_mod_ctx_t ctx)

Sets $A$ to the matrix with the entries of $B$, read as unsigned
integers. No reduction is done, so the modulus of $B$ should be at most
$n$ (normally they are equal).
:::

::: function
void _fmpz_mod_mat_reduce(fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)

Reduces every entry of `mat` modulo $n$, in place; the entries may be
any integers.
:::

## Random matrices {#random}

::: function
void fmpz_mod_mat_rand(fmpz_mod_mat_t mat, flint_rand_t state, const fmpz_mod_ctx_t ctx)

Sets the entries of `mat` to uniformly random residues in $[0, n)$.
:::

::: function
void fmpz_mod_mat_randfull(fmpz_mod_mat_t mat, flint_rand_t state, const fmpz_mod_ctx_t ctx)

Sets the entries to uniformly random nonzero residues in $[1, n)$ (to
zero if $n = 1$).
:::

::: function
void fmpz_mod_mat_randtest(fmpz_mod_mat_t mat, flint_rand_t state, const fmpz_mod_ctx_t ctx)

Sets `mat` to a random matrix for testing: `fmpz_mat_randtest` with the
bit size of $n$ (which favours special values and produces sparse
matrices with probability 1/2), reduced modulo $n$.
:::

::: function
void fmpz_mod_mat_randrank(fmpz_mod_mat_t mat, flint_rand_t state, slong rank, const fmpz_mod_ctx_t ctx)

Sets `mat` to a random sparse matrix of the given rank: `rank` random
nonzero residues placed in distinct rows and columns, the other entries
zero. For prime $n$ the rank is exactly `rank`. Throws an exception if
`rank` is negative or exceeds a dimension. `fmpz_mod_mat_randops` turns
the result into a dense matrix of the same rank.
:::

::: function
void fmpz_mod_mat_randtril(fmpz_mod_mat_t mat, flint_rand_t state, int unit, const fmpz_mod_ctx_t ctx)
void fmpz_mod_mat_randtriu(fmpz_mod_mat_t mat, flint_rand_t state, int unit, const fmpz_mod_ctx_t ctx)

Set `mat` to a random lower, respectively upper, triangular matrix with
uniformly random entries below (above) the diagonal and nonzero diagonal
entries, which are 1 if `unit` is nonzero.
:::

::: function
void fmpz_mod_mat_randops(fmpz_mod_mat_t mat, flint_rand_t state, slong count, const fmpz_mod_ctx_t ctx)

Applies up to `count` random elementary operations (additions and
subtractions of distinct rows or columns, `fmpz_mat_randops`) to `mat`
and reduces the result. The rank and the determinant are unchanged.
:::

## Windows and concatenation {#windows}

::: function
void fmpz_mod_mat_window_init(fmpz_mod_mat_t window, const fmpz_mod_mat_t mat, slong r1, slong c1, slong r2, slong c2, const fmpz_mod_ctx_t ctx)
void fmpz_mod_mat_window_clear(fmpz_mod_mat_t window, const fmpz_mod_ctx_t ctx)

`fmpz_mod_mat_window_init` initialises `window` as the
$(r_2 - r_1) \times (c_2 - c_1)$ submatrix of `mat` with top left entry
in row $r_1$ and column $c_1$, sharing its entries with `mat`.
`fmpz_mod_mat_window_clear` releases the window without freeing the
entries.
:::

::: function
void fmpz_mod_mat_concat_horizontal(fmpz_mod_mat_t res, const fmpz_mod_mat_t mat1, const fmpz_mod_mat_t mat2, const fmpz_mod_ctx_t ctx)
void fmpz_mod_mat_concat_vertical(fmpz_mod_mat_t res, const fmpz_mod_mat_t mat1, const fmpz_mod_mat_t mat2, const fmpz_mod_ctx_t ctx)

Set `res` to $(\mathit{mat1} \; \mathit{mat2})$, the $m \times (k + l)$
matrix with `mat1` ($m \times k$) on the left and `mat2` ($m \times l$)
on the right, respectively to the $(m + k) \times l$ matrix with `mat1`
($m \times l$) above `mat2` ($k \times l$).
:::

## Input and output {#io}

::: function
int fmpz_mod_mat_fprint(FILE * file, const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)
int fmpz_mod_mat_print(const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)

Print `mat` to `file` or to `stdout` in the format of `fmpz_mat_fprint`:
the numbers of rows and columns, two spaces, and the entries in
row-major order separated by spaces. Return a positive value on success
and a nonpositive value if a write failed.
:::

::: function
int fmpz_mod_mat_fprint_pretty(FILE * file, const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)
void fmpz_mod_mat_print_pretty(const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)

Print `mat` to `file` or to `stdout` in the format of
`fmpz_mat_fprint_pretty`: rows in square brackets, one per line, inside
an outer pair of brackets.
:::

## Comparison {#comparison}

::: function
int fmpz_mod_mat_equal(const fmpz_mod_mat_t mat1, const fmpz_mod_mat_t mat2, const fmpz_mod_ctx_t ctx)

Returns 1 if the matrices have the same dimensions and entries, and 0
otherwise.
:::

::: function
int fmpz_mod_mat_is_zero(const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)
int fmpz_mod_mat_is_one(const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)

Return 1 if `mat` is the zero matrix, respectively has ones on the main
diagonal and zeros elsewhere (always true for $n = 1$), and 0 otherwise.
:::

## Addition and scalar multiplication {#addition}

::: function
void fmpz_mod_mat_add(fmpz_mod_mat_t C, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx)
void fmpz_mod_mat_sub(fmpz_mod_mat_t C, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx)
void fmpz_mod_mat_neg(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, const fmpz_mod_ctx_t ctx)

Set $C = A + B$, $C = A - B$ and $B = -A$, row by row with the
`fmpz_mod_vec` functions.
:::

::: function
void fmpz_mod_mat_scalar_mul_si(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, slong c, const fmpz_mod_ctx_t ctx)
void fmpz_mod_mat_scalar_mul_ui(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, ulong c, const fmpz_mod_ctx_t ctx)
void fmpz_mod_mat_scalar_mul_fmpz(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, fmpz_t c, const fmpz_mod_ctx_t ctx)

Set $B = cA$ for any integer $c$. The products are computed exactly and
reduced.
:::

## Matrix multiplication {#mul}

::: function
void fmpz_mod_mat_mul(fmpz_mod_mat_t C, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx)
void fmpz_mod_mat_sqr(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, const fmpz_mod_ctx_t ctx)

Set $C = AB$, respectively $B = A^2$ for square $A$. The number of
columns of $A$ must equal the number of rows of $B$.

::: details
The exact integer product is computed with `fmpz_mat_mul`
(respectively `fmpz_mat_sqr`), which chooses between classical,
Strassen, multimodular and BLAS-based algorithms and is multithreaded,
and each entry is then reduced modulo $n$
([](algorithms/linear-algebra/matrix-multiplication)).

:::
:::

::: function
void fmpz_mod_mat_mul_classical_threaded(fmpz_mod_mat_t C, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx)

Sets $C = AB$ with the classical algorithm: each entry is an exact dot
product (`_fmpz_vec_dot` on a transposed copy of $B$) reduced once.
Blocks of entries are distributed over the available threads when all
dimensions are at least 10. `fmpz_mod_mat_mul` does not call this
function.
:::

::: function
void _fmpz_mod_mat_mul_classical_threaded_pool_op(fmpz_mod_mat_t D, const fmpz_mod_mat_t C, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B, int op, thread_pool_handle * threads, slong num_threads, const fmpz_mod_ctx_t ctx)

Sets $D = AB + \mathit{op} \cdot C$ for $\mathit{op} \in \{-1, 0, 1\}$
($C$ is not read for $\mathit{op} = 0$), with the classical algorithm on
the calling thread and the `num_threads` threads of the pool given by
`threads`. All dimensions must be positive.
:::

::: function
void fmpz_mod_mat_submul(fmpz_mod_mat_t D, const fmpz_mod_mat_t C, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx)

Sets $D = C - AB$, using `fmpz_mod_mat_mul` and a temporary matrix.
:::

::: function
void fmpz_mod_mat_pow_ui(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, ulong e, const fmpz_mod_ctx_t ctx)

Sets $B = A^e$ for a square matrix $A$, with $A^0$ the identity. A
wrapper of `gr_mat_pow_ui` (binary exponentiation). Aliasing is
allowed.
:::

::: function
void fmpz_mod_mat_mul_fmpz_vec(fmpz * c, const fmpz_mod_mat_t A, const fmpz * b, slong blen, const fmpz_mod_ctx_t ctx)
void fmpz_mod_mat_mul_fmpz_vec_ptr(fmpz * const * c, const fmpz_mod_mat_t A, const fmpz * const * b, slong blen, const fmpz_mod_ctx_t ctx)

Set the column vector $c$, of length the number of rows of $A$, to $Ab$
modulo $n$, where $b$ is truncated or zero-extended to the number of
columns of $A$. The `ptr` variant takes arrays of pointers to the
entries. $b$ may contain any integers; $c$ must not overlap $b$. The
products are computed exactly with `fmpz_mat_mul_fmpz_vec` and reduced.
:::

::: function
void fmpz_mod_mat_fmpz_vec_mul(fmpz * c, const fmpz * a, slong alen, const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx)
void fmpz_mod_mat_fmpz_vec_mul_ptr(fmpz * const * c, const fmpz * const * a, slong alen, const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx)

Set the row vector $c$, of length the number of columns of $B$, to $aB$
modulo $n$, where $a$ is truncated or zero-extended to the number of
rows of $B$. As above, $a$ may contain any integers and $c$ must not
overlap $a$.
:::

## Trace {#trace}

::: function
void fmpz_mod_mat_trace(fmpz_t trace, const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)

Sets `trace` to the sum of the diagonal entries of the square matrix
`mat`, modulo $n$.
:::

## Permutations {#permutations}

::: function
void fmpz_mod_mat_swap_rows(fmpz_mod_mat_t mat, slong * perm, slong r, slong s, const fmpz_mod_ctx_t ctx)
void fmpz_mod_mat_swap_cols(fmpz_mod_mat_t mat, slong * perm, slong r, slong s, const fmpz_mod_ctx_t ctx)

Swap rows (columns) $r$ and $s$ of `mat` and, if `perm` is not `NULL`,
entries $r$ and $s$ of the permutation array `perm`.
:::

::: function
void fmpz_mod_mat_invert_rows(fmpz_mod_mat_t mat, slong * perm, const fmpz_mod_ctx_t ctx)
void fmpz_mod_mat_invert_cols(fmpz_mod_mat_t mat, slong * perm, const fmpz_mod_ctx_t ctx)

Reverse the order of the rows (columns) of `mat` and, if `perm` is not
`NULL`, of the entries of `perm`.
:::

## Gaussian elimination {#elimination}

These functions are wrappers of the generic `gr_mat` algorithms over the
`fmpz_mod` ring, which in turn use `fmpz_mod_mat_mul` for the products
in block recursive algorithms. Unless stated otherwise they require $n$
to be prime, which is not checked. For a composite modulus they work as
long as every pivot they meet is invertible; otherwise they throw an
exception or, for the functions returning a success flag,
return 0 even if the result exists. See [](fmpz_mod_mat/composite).

::: function
slong fmpz_mod_mat_lu(slong * P, fmpz_mod_mat_t A, int rank_check, const fmpz_mod_ctx_t ctx)

Computes a generalised LU decomposition $PA = LU$ of $A$ in place and
returns the rank of $A$. *P* must have room for as many entries as $A$
has rows; it receives the row permutation. A wrapper of `gr_mat_lu`.

::: details
If $A$ is a nonsingular square matrix, it is overwritten with a unit
lower triangular matrix $L$ (whose diagonal is not stored) and an upper
triangular matrix $U$. For a matrix of rank $r$, $U$ is in row echelon
form with $r$ nonzero rows, and $L$ is lower triangular, truncated to
$r$ columns, with implicit ones on its first $r$ diagonal entries.

If `rank_check` is nonzero, the function gives up as soon as $A$ is
found to be rank-deficient, returns 0 and leaves $A$ in an undefined
state.

:::
:::

::: function
void fmpz_mod_mat_solve_tril(fmpz_mod_mat_t X, const fmpz_mod_mat_t L, const fmpz_mod_mat_t B, int unit, const fmpz_mod_ctx_t ctx)
void fmpz_mod_mat_solve_triu(fmpz_mod_mat_t X, const fmpz_mod_mat_t U, const fmpz_mod_mat_t B, int unit, const fmpz_mod_ctx_t ctx)

Set $X = L^{-1} B$, respectively $X = U^{-1} B$, for a square lower
(upper) triangular matrix whose diagonal entries are invertible; if
`unit` is nonzero, the diagonal is taken to be 1 and not read. Only the
lower (upper) triangle is read. $X$ may alias $B$. Wrappers of
`gr_mat_nonsingular_solve_tril` and `gr_mat_nonsingular_solve_triu`,
which choose between substitution and a block recursive algorithm; they
throw an exception if a diagonal entry is not invertible.
:::

::: function
int fmpz_mod_mat_solve(fmpz_mod_mat_t X, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx)

Solves $AX = B$ for a square matrix $A$. Returns 1 if $A$ is
nonsingular; otherwise returns 0 and leaves $X$ undefined. A wrapper of
`gr_mat_nonsingular_solve_lu`.
:::

::: function
int fmpz_mod_mat_can_solve(fmpz_mod_mat_t X, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx)

Solves $AX = B$ for a matrix $A$ of any shape and rank. Returns 1 and
sets $X$ to one solution if a solution exists; otherwise returns 0 and
leaves $X$ undefined. A wrapper of `gr_mat_solve_field`.
:::

::: function
int fmpz_mod_mat_inv(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, const fmpz_mod_ctx_t ctx)

Sets $B = A^{-1}$ and returns 1 if the square matrix $A$ is invertible;
otherwise returns 0 and leaves $B$ undefined. $A$ and $B$ must have the
same dimensions. The inverse is computed by solving $AB = I$ with
`fmpz_mod_mat_solve`; a $1 \times 1$ matrix is inverted with
`fmpz_mod_inv`, which throws an exception for a nonzero non-unit.
:::

::: function
slong fmpz_mod_mat_rank(const fmpz_mod_mat_t A, const fmpz_mod_ctx_t ctx)

Returns the rank of $A$, computed by LU decomposition (`gr_mat_rank_lu`).
:::

::: function
slong fmpz_mod_mat_rref(fmpz_mod_mat_t res, const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)

Sets `res` to the reduced row echelon form of `mat`, which has the same
dimensions, and returns the rank. A wrapper of `gr_mat_rref_lu`.
:::

::: function
slong fmpz_mod_mat_nullspace(fmpz_mod_mat_t X, const fmpz_mod_mat_t A, const fmpz_mod_ctx_t ctx)

Computes a basis of the right nullspace $\{x : Ax = 0\}$ of $A$ and
returns its dimension (the nullity) $k$. The basis vectors are written
to the first $k$ columns of $X$, and all other entries of $X$ are set to
zero. $X$ is not resized: it must have as many rows as $A$ has columns
and at least $k$ columns (allocate it square, with the number of
columns of $A$, when $k$ is not known).
:::

::: function
void fmpz_mod_mat_det(fmpz_t res, const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)

Sets `res` to the determinant of the square matrix `mat`, and throws an
exception if `mat` is not square. A wrapper of `gr_mat_det`; this works
for every modulus, since a failed LU decomposition falls back to a
division-free algorithm.
:::

::: function
int fmpz_mod_mat_reduce_row(slong * column, fmpz_mod_mat_t A, slong * P, slong * L, slong m, const fmpz_mod_ctx_t ctx)

Reduces row $m$ of $A$ against the earlier rows, with the conventions of
`gr_mat_reduce_row`: `P[i]` is the row with its pivot in column $i$, or
$-1$; `L[i]` bounds the nonzero columns of row $i$. Sets `*column` to the
pivot column of the reduced row (made monic and recorded in *P*), or to
$-1$ if the row becomes zero. Returns `GR_SUCCESS`, or `GR_DOMAIN` if
the new pivot is not invertible. The sums are accumulated exactly and
each entry is reduced once, when it is reached.
:::

## Howell form and strong echelon form {#howell}

These functions work for every modulus. They are wrappers of
`fmpz_mat_strong_echelon_form_mod` and `fmpz_mat_howell_form_mod`; see
[](algorithms/linear-algebra/hermite-smith) and
[](algorithms/linear-algebra/gaussian-elimination) for the algorithms.

::: function
void fmpz_mod_mat_strong_echelon_form(fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)

Transforms `mat` into its strong echelon form [@FieHof2014], which
equals the Howell form up to the order of the rows. `mat` must have at
least as many rows as columns.
:::

::: function
slong fmpz_mod_mat_howell_form(fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)

Transforms `mat` into its Howell form [@How1986; @StoMul1998] and
returns the number of nonzero rows. `mat` must have at least as many
rows as columns. Two matrices generate the same submodule of
$(\Z/n\Z)^c$ if and only if their Howell forms are equal.
:::

## Similarity transforms and polynomials {#charpoly}

::: function
void fmpz_mod_mat_similarity(fmpz_mod_mat_t A, slong r, fmpz_t d, const fmpz_mod_ctx_t ctx)

Applies a similarity transform to the square matrix $A$ in place: with
$P$ the identity matrix whose row $r$ (0-indexed) has its zero entries
replaced by $d$, sets $A = P^{-1} A P$. $d$ must be canonical. A wrapper
of `gr_mat_apply_row_similarity`; it uses no divisions and works for
every modulus. Similarity transforms preserve the determinant and the
characteristic and minimal polynomials.
:::

::: function
void _fmpz_mod_mat_charpoly(fmpz * p, const fmpz_mod_mat_t M, const fmpz_mod_ctx_t ctx)
void fmpz_mod_mat_charpoly(fmpz_mod_poly_t p, const fmpz_mod_mat_t M, const fmpz_mod_ctx_t ctx)

Set $p$ to the characteristic polynomial $\det(xI - M)$ of the square
matrix $M$; the underscore version writes its $m + 1$ coefficients for an
$m \times m$ matrix. `fmpz_mod_mat_charpoly` throws an exception if $M$
is not square. These work for every modulus: the generic `fmpz_mod`
ring uses Danilevsky's algorithm for $m \ge 16$ and falls back to the
division-free Berkowitz algorithm when it fails, and uses Berkowitz for
smaller $m$.
:::

::: function
void fmpz_mod_mat_minpoly(fmpz_mod_poly_t p, const fmpz_mod_mat_t M, const fmpz_mod_ctx_t ctx)

Sets $p$ to the minimal polynomial of the square matrix $M$, using
`gr_mat_minpoly_field`. Throws an exception if $M$ is not square. The
modulus must be prime; the ring is declared a field without checking.
:::
