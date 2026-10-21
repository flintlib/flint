---
title: ca_mat.h
header: ca_mat.h
description: Dense matrices with exact real and complex entries: arithmetic, LU decomposition, solving, rank, determinants, characteristic polynomials, eigenvalues, Jordan forms and matrix functions.
category: Exact real and complex numbers
order: 5
---

This module implements dense matrices whose entries are exact real or
complex numbers of type `ca_t` (see [](ca)). It offers arithmetic,
Gaussian elimination, solving, rank, reduced row echelon forms and
nullspaces, determinants and characteristic polynomials, and, building
on [](ca_poly) root finding, eigenvalues with multiplicities,
diagonalisation, the Jordan canonical form and the matrix exponential
and logarithm, all in exact arithmetic.

Most functions are wrappers around the generic linear algebra of
[](gr_mat) with a `ca` context. The module adds a dispatch layer that
recognises matrices with rational entries, which are handed to
[](fmpz_mat) or [](fmpq_mat), and matrices whose entries lie in one
number field, for which multiplication is done with integer polynomial
matrices and elimination with divisions is preferred. The algorithms,
and the difficulties caused by zero tests that may fail, are discussed
in [](algorithms/exact/linear-algebra).

::: summary algorithms/exact/linear-algebra :::

## Types and conventions {#conventions}

::: type
ca_mat_struct
ca_mat_t

A matrix of `ca_t` entries. The structure has the fields `entries` (a
pointer to the first entry), `r` and `c` (the numbers of rows and
columns) and `stride` (the distance in entries between the starts of
consecutive rows). Entry $(i, j)$ is `entries[i * stride + j]`. A
`ca_mat_t` is an array of one structure, so it is passed by reference.
:::

- The dimensions are fixed at initialisation; either may be zero. The
  caller must pass outputs of compatible dimensions.
- Functions allow aliasing between inputs and outputs unless stated
  otherwise.
- Entries should be numbers; special values (*Unknown* in particular)
  make most algorithms fail or return *Unknown* entries.
- Operations that must decide whether an entry is zero (pivoting, rank,
  inversion, eigenvalues) use `ca_check_is_zero` and fail or return
  `T_UNKNOWN` when it cannot decide. Such functions return an `int`
  success flag or a `truth_t`. Entries that are proved to be zero during
  pivot search are overwritten with an exact zero.

::: macro
ca_mat_entry(mat, i, j)
ca_mat_nrows(mat)
ca_mat_ncols(mat)

`ca_mat_entry` expands to a pointer (`ca_ptr`) to the entry in row `i`
and column `j` (indices not checked); the other macros give the numbers
of rows and columns.
:::

::: function
ca_ptr ca_mat_entry_ptr(ca_mat_t mat, slong i, slong j)

Returns a pointer to the entry in row `i` and column `j`, like
`ca_mat_entry`.
:::

::: function
int ca_mat_is_empty(const ca_mat_t mat)
int ca_mat_is_square(const ca_mat_t mat)

Returns whether `mat` has no entries (zero rows or zero columns), or
whether it is square.
:::

## Memory management {#memory}

::: function
void ca_mat_init(ca_mat_t mat, slong r, slong c, ca_ctx_t ctx)

Initialises `mat` as the zero matrix with `r` rows and `c` columns.
:::

::: function
void ca_mat_clear(ca_mat_t mat, ca_ctx_t ctx)

Clears `mat`, freeing its entries.
:::

::: function
void ca_mat_swap(ca_mat_t mat1, ca_mat_t mat2, ca_ctx_t ctx)

Swaps `mat1` and `mat2` in constant time.
:::

::: function
void ca_mat_window_init(ca_mat_t window, const ca_mat_t mat, slong r1, slong c1, slong r2, slong c2, ca_ctx_t ctx)

Initialises `window` as a view of the submatrix of `mat` formed by rows
$r_1 \le i < r_2$ and columns $c_1 \le j < c_2$, sharing the entries of
`mat` (nothing is copied). The window must not outlive `mat`.
:::

::: function
void ca_mat_window_clear(ca_mat_t window, ca_ctx_t ctx)

Clears a window. This does nothing, since a window owns no memory.
:::

## Assignment and conversions {#assignment}

::: function
void ca_mat_set(ca_mat_t dest, const ca_mat_t src, ca_ctx_t ctx)
void ca_mat_set_fmpz_mat(ca_mat_t dest, const fmpz_mat_t src, ca_ctx_t ctx)
void ca_mat_set_fmpq_mat(ca_mat_t dest, const fmpq_mat_t src, ca_ctx_t ctx)

Sets `dest` to `src`, which must have the same dimensions.
:::

::: function
void ca_mat_set_ca(ca_mat_t y, const ca_t x, ca_ctx_t ctx)

Sets `y` to the matrix with `x` on the main diagonal and zeros
elsewhere (it need not be square).
:::

::: function
void ca_mat_transfer(ca_mat_t res, ca_ctx_t res_ctx, const ca_mat_t src, ca_ctx_t src_ctx)

Sets `res`, belonging to the context `res_ctx`, to the matrix `src`,
belonging to `src_ctx`, transferring each entry with `ca_transfer`.
:::

::: function
void ca_mat_get_fexpr(fexpr_t res, const ca_mat_t A, ulong flags, ca_ctx_t ctx)

Sets `res` to a symbolic expression `Matrix(Row(...), ..., Row(...))`
for `A`, wrapped in `Where(..., Def(a_1, ...), ...)` defining the extension
numbers that occur, as for `ca_get_fexpr`, which also explains `flags`.
:::

## Random generation {#random}

::: function
void ca_mat_randtest(ca_mat_t mat, flint_rand_t state, slong len, slong bits, ca_ctx_t ctx)

Sets `mat` to a random matrix whose entries are generated by
`ca_randtest` with depth `len` and size parameter `bits`.
:::

::: function
void ca_mat_randtest_rational(ca_mat_t mat, flint_rand_t state, slong bits, ca_ctx_t ctx)

Sets `mat` to a random matrix with rational entries of up to `bits`
bits.
:::

::: function
void ca_mat_randops(ca_mat_t mat, flint_rand_t state, slong count, ca_ctx_t ctx)

Modifies `mat` in place by at most `count` random additions or
subtractions of one row (or column) to another. This preserves the rank
and, for a square matrix, the determinant. Wrapper of `gr_mat_randops`.
:::

## Input and output {#io}

::: function
void ca_mat_print(const ca_mat_t mat, ca_ctx_t ctx)

Prints the line `ca_mat of size r x c:` followed by the entries in row
order, one per line, printed with `ca_print` and indented by four
spaces, and a blank line.
:::

::: function
void ca_mat_printn(const ca_mat_t mat, slong digits, ca_ctx_t ctx)

Prints decimal approximations of the entries with `digits` significant
digits as a nested list `[[a, b],\n[c, d]]`, one row per line, followed
by a newline.
:::

## Special matrices {#special}

::: function
void ca_mat_zero(ca_mat_t mat, ca_ctx_t ctx)
void ca_mat_one(ca_mat_t mat, ca_ctx_t ctx)
void ca_mat_ones(ca_mat_t mat, ca_ctx_t ctx)

Sets `mat` to the zero matrix, to the matrix with ones on the main
diagonal and zeros elsewhere, or to the matrix with all entries 1.
:::

::: function
void ca_mat_pascal(ca_mat_t mat, int triangular, ca_ctx_t ctx)

Sets `mat` to a Pascal matrix of binomial coefficients: for
`triangular` = 0 the symmetric matrix with entries $\binom{i+j}{i}$,
for 1 the upper triangular matrix with entries $\binom{j}{i}$, and for
$-1$ the lower triangular matrix with entries $\binom{i}{j}$. Wrapper of
`gr_mat_pascal`.
:::

::: function
void ca_mat_stirling(ca_mat_t mat, int kind, ca_ctx_t ctx)

Sets `mat` to the matrix of Stirling numbers $S(i, j)$: unsigned
Stirling numbers of the first kind for `kind` = 0, signed ones for
`kind` = 1, Stirling numbers of the second kind for `kind` = 2. Wrapper
of `gr_mat_stirling`.
:::

::: function
void ca_mat_hilbert(ca_mat_t mat, ca_ctx_t ctx)

Sets `mat` to the Hilbert matrix with entries $1/(i + j + 1)$, $i, j \ge
0$. Wrapper of `gr_mat_hilbert`.
:::

::: function
void ca_mat_dft(ca_mat_t res, int type, ca_ctx_t ctx)

Sets `res` to a matrix of the discrete Fourier transform of length $n$,
the smaller of the two dimensions (along the larger dimension the
entries repeat periodically). With $\omega = e^{2 \pi i / n}$, the
entries $A_{j,k}$ are $\omega^{-jk}$ for `type` = 0, $\omega^{jk}/n$ for
type 1, $\omega^{-jk}/\sqrt{n}$ for type 2 and $\omega^{jk}/\sqrt{n}$
for type 3. Types 0 and 1 are inverse to each other, and so are types 2
and 3. The entries are elements of the cyclotomic field
$\mathbb{Q}(\omega)$ (and $\sqrt{n}$).
:::

## Comparisons {#comparisons}

::: function
truth_t ca_mat_check_equal(const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)
truth_t ca_mat_check_is_zero(const ca_mat_t A, ca_ctx_t ctx)
truth_t ca_mat_check_is_one(const ca_mat_t A, ca_ctx_t ctx)

Checks whether $A = B$ entrywise (`T_FALSE` if the dimensions differ),
whether $A$ is zero, or whether $A$ has ones on the main diagonal and
zeros elsewhere. Wrappers of `gr_mat_equal`, `gr_mat_is_zero` and
`gr_mat_is_one`.
:::

## Conjugate and transpose {#transpose}

::: function
void ca_mat_transpose(ca_mat_t B, const ca_mat_t A, ca_ctx_t ctx)

Sets `B` to the transpose of `A`. Aliasing is allowed for square
matrices. Wrapper of `gr_mat_transpose`.
:::

::: function
void ca_mat_conj(ca_mat_t B, const ca_mat_t A, ca_ctx_t ctx)
void ca_mat_conj_transpose(ca_mat_t mat1, const ca_mat_t mat2, ca_ctx_t ctx)

Sets the output to the entrywise complex conjugate of the input, or to
its conjugate (Hermitian) transpose. Entries are conjugated with
`ca_conj`.
:::

## Arithmetic {#arithmetic}

::: function
void ca_mat_neg(ca_mat_t dest, const ca_mat_t src, ca_ctx_t ctx)
void ca_mat_add(ca_mat_t res, const ca_mat_t mat1, const ca_mat_t mat2, ca_ctx_t ctx)
void ca_mat_sub(ca_mat_t res, const ca_mat_t mat1, const ca_mat_t mat2, ca_ctx_t ctx)

Sets the output to the negation of `src`, or to the sum or difference of
`mat1` and `mat2`.
:::

::: function
void ca_mat_add_ca(ca_mat_t y, const ca_mat_t a, const ca_t x, ca_ctx_t ctx)
void ca_mat_sub_ca(ca_mat_t y, const ca_mat_t a, const ca_t x, ca_ctx_t ctx)

Sets `y` to $a + x I$ or $a - x I$, that is, adds or subtracts the
scalar `x` on the main diagonal.
:::

::: function
void ca_mat_addmul_ca(ca_mat_t y, const ca_mat_t a, const ca_t x, ca_ctx_t ctx)
void ca_mat_submul_ca(ca_mat_t y, const ca_mat_t a, const ca_t x, ca_ctx_t ctx)

Sets `y` to $y + x a$ or $y - x a$.
:::

::: function
void ca_mat_mul_si(ca_mat_t B, const ca_mat_t A, slong c, ca_ctx_t ctx)
void ca_mat_mul_fmpz(ca_mat_t B, const ca_mat_t A, const fmpz_t c, ca_ctx_t ctx)
void ca_mat_mul_fmpq(ca_mat_t B, const ca_mat_t A, const fmpq_t c, ca_ctx_t ctx)
void ca_mat_mul_ca(ca_mat_t B, const ca_mat_t A, const ca_t c, ca_ctx_t ctx)
void ca_mat_div_si(ca_mat_t B, const ca_mat_t A, slong c, ca_ctx_t ctx)
void ca_mat_div_fmpz(ca_mat_t B, const ca_mat_t A, const fmpz_t c, ca_ctx_t ctx)
void ca_mat_div_fmpq(ca_mat_t B, const ca_mat_t A, const fmpq_t c, ca_ctx_t ctx)
void ca_mat_div_ca(ca_mat_t B, const ca_mat_t A, const ca_t c, ca_ctx_t ctx)

Sets `B` to `A` multiplied or divided by the scalar `c`, entry by entry.
Division by zero gives special values as entries.
:::

::: function
void ca_mat_mul_classical(ca_mat_t C, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)
void ca_mat_mul_same_nf(ca_mat_t C, const ca_mat_t A, const ca_mat_t B, ca_field_t K, ca_ctx_t ctx)
void ca_mat_mul(ca_mat_t C, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)

Sets `C` to the matrix product $AB$. The dimensions must be compatible;
otherwise an exception is raised.

::: details
The *classical* version computes each entry as a dot product (wrapper of
`gr_mat_mul_classical`). The *same_nf* version requires, without
checking, that all entries of $A$ and $B$ are rational numbers or
elements of the number field `K` = $\mathbb{Q}(a)$; it writes the
entries as polynomials in $a$ over common denominators, multiplies the
integer polynomial matrices with `fmpz_poly_mat_mul` and reduces modulo
the minimal polynomial of $a$. The default version converts matrices
with rational entries (and at least 3 rows of $B$) to [](fmpz_mat) or
[](fmpq_mat), uses the *same_nf* version when all entries lie in one
number field and the dimensions are at least about 3, and otherwise
the classical algorithm.
:::
:::

::: function
void ca_mat_sqr(ca_mat_t res, const ca_mat_t A, ca_ctx_t ctx)

Sets `res` to $A^2$.
:::

::: function
void ca_mat_pow_ui_binexp(ca_mat_t B, const ca_mat_t A, ulong exp, ca_ctx_t ctx)

Sets `B` to $A^{\mathit{exp}}$ for a square matrix `A` (not checked),
by binary exponentiation.
:::

## Polynomial evaluation and trace {#evaluation}

::: function
void _ca_mat_ca_poly_evaluate(ca_mat_t y, ca_srcptr poly, slong len, const ca_mat_t x, ca_ctx_t ctx)
void ca_mat_ca_poly_evaluate(ca_mat_t res, const ca_poly_t f, const ca_mat_t a, ca_ctx_t ctx)

Sets the output to $f(A)$ for the polynomial $f$ (given as an array of
`len` coefficients for the underscore version) and the square matrix
$A$, using the Paterson–Stockmeyer algorithm, which needs about
$2\sqrt{n}$ matrix multiplications for a polynomial of degree $n$.
Wrappers of `_gr_mat_gr_poly_evaluate`.
:::

::: function
void ca_mat_trace(ca_t trace, const ca_mat_t mat, ca_ctx_t ctx)

Sets `trace` to the sum of the diagonal entries of the square matrix
`mat`. Wrapper of `gr_mat_trace`.
:::

## Gaussian elimination and LU decomposition {#lu}

::: function
truth_t ca_mat_find_pivot(slong * pivot_row, ca_mat_t mat, slong start_row, slong end_row, slong column, ca_ctx_t ctx)

Looks for a nonzero entry in column `column` of `mat` with row index in
$[\mathit{start\_row}, \mathit{end\_row})$, which must be a nonempty
range. Returns `T_TRUE` and sets `pivot_row` if one is found, `T_FALSE`
if all these entries are proved zero, and `T_UNKNOWN` otherwise.

::: details
Among the entries that are not rational zeros, the one with the
simplest representation (smallest in the order of `ca_cmp_repr`) is
tested first, since it is most likely to be proved nonzero quickly; the
other entries are tested if this fails. The function modifies `mat`:
entries that are proved zero are replaced by an exact zero.
:::
:::

::: function
void _ca_mat_swap_rows(ca_mat_t mat, slong * perm, slong r, slong s)

Swaps rows `r` and `s` of `mat` (by swapping entries) and, if `perm` is
not `NULL`, the corresponding entries of the permutation array `perm`.
:::

::: function
int ca_mat_lu_classical(slong * rank, slong * P, ca_mat_t LU, const ca_mat_t A, int rank_check, ca_ctx_t ctx)
int ca_mat_lu_recursive(slong * rank, slong * P, ca_mat_t LU, const ca_mat_t A, int rank_check, ca_ctx_t ctx)
int ca_mat_lu(slong * rank, slong * P, ca_mat_t LU, const ca_mat_t A, int rank_check, ca_ctx_t ctx)

Computes a generalised LU decomposition $A = PLU$ of the $m \times n$
matrix `A`, writes the rank to `rank` and the row permutation to the
array `P` of length $m$, and returns 1. Returns 0 if a pivot can neither
be proved zero nor nonzero, in which case the outputs are meaningless.

::: details
If $A$ is square and nonsingular, `LU` holds a unit lower triangular
$L$ (its diagonal not stored) and an upper triangular $U$. If $A$ has
rank $r$, $U$ is in row echelon form with $r$ nonzero rows and $L$ is
lower triangular, truncated to $r$ columns with implicit ones on the
first $r$ diagonal entries. If `rank_check` is nonzero, the computation
stops as soon as $A$ is found to be rank-deficient, setting the rank to
0 and leaving `LU` undefined.

These are wrappers of `gr_mat_lu_classical` (Gaussian elimination),
`gr_mat_lu_recursive` (block recursive elimination, which benefits from
fast matrix multiplication) and `gr_mat_lu_generic` (an automatic
choice).
:::
:::

::: function
int ca_mat_fflu(slong * rank, slong * P, ca_mat_t LU, ca_t den, const ca_mat_t A, int rank_check, ca_ctx_t ctx)

Like `ca_mat_lu`, but computes a fraction-free LU decomposition with
Bareiss's algorithm and writes its denominator (the last pivot) to
`den`. Although no divisions by pivots that do not cancel exactly are
performed, fractions can appear in the entries since the symbolic
simplification of exact quotients is not complete.
:::

::: function
truth_t ca_mat_nonsingular_lu(slong * P, ca_mat_t LU, const ca_mat_t A, ca_ctx_t ctx)
truth_t ca_mat_nonsingular_fflu(slong * P, ca_mat_t LU, ca_t den, const ca_mat_t A, ca_ctx_t ctx)

For a square matrix `A`: returns `T_TRUE` and an LU (or fraction-free
LU) decomposition $A = PLU$ if $A$ is proved nonsingular, `T_FALSE` if
it is proved singular, and `T_UNKNOWN` otherwise. Unless the result is
`T_TRUE`, the outputs are meaningless.
:::

## Solving and inverse {#solving}

::: function
truth_t ca_mat_nonsingular_solve_adjugate(ca_mat_t X, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)
truth_t ca_mat_nonsingular_solve_fflu(ca_mat_t X, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)
truth_t ca_mat_nonsingular_solve_lu(ca_mat_t X, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)
truth_t ca_mat_nonsingular_solve(ca_mat_t X, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)

For a square matrix `A`: if $A$ is proved nonsingular, sets $X =
A^{-1} B$ and returns `T_TRUE`; returns `T_FALSE` if $A$ is proved
singular and `T_UNKNOWN` if this cannot be decided.

::: details
The *adjugate* version computes $X = \operatorname{adj}(A) B / \det(A)$
(see `ca_mat_adjugate`), which needs only one zero test, of the
determinant. The *fflu* and *lu* versions use fraction-free or ordinary
LU decomposition followed by triangular solving; the *lu* version is a
wrapper of `gr_mat_nonsingular_solve_lu`. The default version currently
uses the *lu* version.
:::
:::

::: function
truth_t ca_mat_inv(ca_mat_t X, const ca_mat_t A, ca_ctx_t ctx)

For a square matrix `A`: if $A$ is proved nonsingular, sets $X =
A^{-1}$ and returns `T_TRUE`; returns `T_FALSE` if $A$ is proved
singular and `T_UNKNOWN` otherwise. Matrices of size at most 4, and
matrices whose entries do not all lie in $\mathbb{Q}$ or one number
field, are inverted through the adjugate matrix; the others by LU
decomposition.
:::

::: function
void ca_mat_solve_tril_classical(ca_mat_t X, const ca_mat_t L, const ca_mat_t B, int unit, ca_ctx_t ctx)
void ca_mat_solve_tril_recursive(ca_mat_t X, const ca_mat_t L, const ca_mat_t B, int unit, ca_ctx_t ctx)
void ca_mat_solve_tril(ca_mat_t X, const ca_mat_t L, const ca_mat_t B, int unit, ca_ctx_t ctx)
void ca_mat_solve_triu_classical(ca_mat_t X, const ca_mat_t U, const ca_mat_t B, int unit, ca_ctx_t ctx)
void ca_mat_solve_triu_recursive(ca_mat_t X, const ca_mat_t U, const ca_mat_t B, int unit, ca_ctx_t ctx)
void ca_mat_solve_triu(ca_mat_t X, const ca_mat_t U, const ca_mat_t B, int unit, ca_ctx_t ctx)

Solves $LX = B$ for a lower triangular `L`, or $UX = B$ for an upper
triangular `U`, by substitution (*classical*), block recursively
(*recursive*), or with an automatic choice. The matrix must be
triangular with nonzero diagonal entries; this is not checked. If
`unit` is nonzero, the diagonal is taken to consist of ones and is not
read. Wrappers of the `gr_mat_nonsingular_solve_tril` and `_triu`
functions.
:::

::: function
void ca_mat_solve_lu_precomp(ca_mat_t X, const slong * perm, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)
void ca_mat_solve_fflu_precomp(ca_mat_t X, const slong * perm, const ca_mat_t A, const ca_t den, const ca_mat_t B, ca_ctx_t ctx)

Solves $AX = B$ given a nonsingular LU decomposition (permutation `perm`
and the matrix `A` holding $L$ and $U$), or a fraction-free one with
denominator `den`, as computed by `ca_mat_nonsingular_lu` or
`ca_mat_nonsingular_fflu`. The decomposition is not checked. `X` may be
the same object as `B` but not as `A`.
:::

## Rank, echelon form and nullspace {#rank}

::: function
int ca_mat_rank(slong * rank, const ca_mat_t A, ca_ctx_t ctx)

Computes the rank of `A` by LU decomposition and returns 1, or returns 0
if a pivot cannot be decided. Wrapper of `gr_mat_rank_lu`.
:::

::: function
int ca_mat_rref_fflu(slong * rank, ca_mat_t R, const ca_mat_t A, ca_ctx_t ctx)
int ca_mat_rref_lu(slong * rank, ca_mat_t R, const ca_mat_t A, ca_ctx_t ctx)
int ca_mat_rref(slong * rank, ca_mat_t R, const ca_mat_t A, ca_ctx_t ctx)

Sets `R` to the reduced row echelon form of `A`, writes the rank to
`rank` and returns 1; returns 0, with meaningless outputs, if the rank
cannot be determined. The *fflu* version starts from a fraction-free LU
decomposition, the *lu* version (a wrapper of `gr_mat_rref_lu`) from an
ordinary one, and the default version is a wrapper of `gr_mat_rref`.
:::

::: function
int ca_mat_right_kernel(ca_mat_t X, const ca_mat_t A, ca_ctx_t ctx)

Sets `X` to a matrix whose columns form a basis of the right kernel
(nullspace) of `A` and returns 1, or returns 0 on failure. `X` is
reinitialised to have as many columns as the nullity of `A`. Wrapper of
`gr_mat_nullspace`.
:::

## Determinant and adjugate {#det}

::: function
void ca_mat_det_berkowitz(ca_t det, const ca_mat_t A, ca_ctx_t ctx)
int ca_mat_det_lu(ca_t det, const ca_mat_t A, ca_ctx_t ctx)
int ca_mat_det_bareiss(ca_t det, const ca_mat_t A, ca_ctx_t ctx)
void ca_mat_det_cofactor(ca_t det, const ca_mat_t A, ca_ctx_t ctx)
void ca_mat_det(ca_t det, const ca_mat_t A, ca_ctx_t ctx)

Sets `det` to the determinant of the square matrix `A`.

::: details
- The *berkowitz* version uses the division-free Berkowitz algorithm
  with $O(n^4)$ operations. It needs no zero tests and always succeeds.
- The *cofactor* version uses cofactor expansion; it is meant for
  matrices of size up to 4.
- The *lu* and *bareiss* versions use LU decomposition and fraction-free
  LU decomposition with $O(n^3)$ operations. They need zero tests for
  pivoting and return 0 if one fails (1 on success).
- The default version computes the determinant of a rational matrix
  of size at least 3 over $\Z$ or $\mathbb{Q}$ ([](fmpz_mat),
  [](fmpq_mat)), uses cofactor expansion for $n \le 4$, LU decomposition
  if all entries lie in one number field, and the Berkowitz algorithm
  otherwise.

Different algorithms can give different symbolic representations of
the same determinant, and which one is fastest depends strongly on the
structure of the matrix. All but the Bareiss version are wrappers of
the corresponding `gr_mat_det` functions.
:::
:::

::: function
void ca_mat_adjugate_cofactor(ca_mat_t adj, ca_t det, const ca_mat_t A, ca_ctx_t ctx)
void ca_mat_adjugate_charpoly(ca_mat_t adj, ca_t det, const ca_mat_t A, ca_ctx_t ctx)
void ca_mat_adjugate(ca_mat_t adj, ca_t det, const ca_mat_t A, ca_ctx_t ctx)

Sets `adj` to the adjugate matrix of `A` and `det` to its determinant,
computed together, by cofactor expansion or from the characteristic
polynomial $p$ (using $\operatorname{adj}(A) = \pm q(A)$ where $p(x) =
x q(x) + p(0)$, a wrapper of `gr_mat_adjugate_charpoly`). Neither needs
zero tests. The default version uses cofactors for $n \le 5$.
:::

## Characteristic polynomial {#charpoly}

::: function
void _ca_mat_charpoly_berkowitz(ca_ptr cp, const ca_mat_t mat, ca_ctx_t ctx)
void ca_mat_charpoly_berkowitz(ca_poly_t cp, const ca_mat_t mat, ca_ctx_t ctx)
int _ca_mat_charpoly_danilevsky(ca_ptr cp, const ca_mat_t mat, ca_ctx_t ctx)
int ca_mat_charpoly_danilevsky(ca_poly_t cp, const ca_mat_t mat, ca_ctx_t ctx)
void _ca_mat_charpoly(ca_ptr cp, const ca_mat_t mat, ca_ctx_t ctx)
void ca_mat_charpoly(ca_poly_t cp, const ca_mat_t mat, ca_ctx_t ctx)

Sets `cp` to the characteristic polynomial $\det(xI - A)$ of the square
matrix `mat`. The underscore versions write $n + 1$ coefficients.

::: details
The *berkowitz* version is division-free, uses $O(n^4)$ operations and
always succeeds. The *danilevsky* version uses $O(n^3)$ operations with
divisions; it needs zero tests and returns 0 if one fails (1 on
success). Both are wrappers of the corresponding `gr_mat` functions.
The default version currently always uses the Berkowitz algorithm.
:::
:::

::: function
int ca_mat_companion(ca_mat_t A, const ca_poly_t poly, ca_ctx_t ctx)

Sets `A` to the companion matrix of `poly` and returns 1, or returns 0
if the leading coefficient of `poly` is not proved nonzero or `A` does
not have size $\deg(\mathit{poly})$. Wrapper of `gr_mat_companion`.
:::

## Eigenvalues and Jordan form {#eigenvalues}

::: function
int ca_mat_eigenvalues(ca_vec_t lambda, ulong * exp, const ca_mat_t mat, ca_ctx_t ctx)

Computes the distinct eigenvalues of the square matrix `mat`, in no
particular order, with their algebraic multiplicities in `exp` (room for
$n$ entries), and returns 1; returns 0 on failure, leaving the outputs
unspecified. The eigenvalues are the roots of the characteristic
polynomial computed by `ca_poly_roots`, so the limitations of that
function apply: in particular, a matrix with rational entries always
succeeds, with eigenvalues that are algebraic numbers.
:::

::: function
truth_t ca_mat_diagonalization(ca_mat_t D, ca_mat_t P, const ca_mat_t A, ca_ctx_t ctx)

Attempts to compute a diagonal matrix `D` and an invertible matrix `P`
with $A = P D P^{-1}$. Returns `T_TRUE` on success, `T_FALSE` if $A$ is
proved not to be diagonalisable (a geometric multiplicity is smaller
than the algebraic one), and `T_UNKNOWN` if the eigenvalues or the
eigenspaces cannot be computed. The columns of $P$ are bases of the
kernels of $A - \lambda I$.
:::

::: function
int ca_mat_jordan_blocks(ca_vec_t lambda, slong * num_blocks, slong * block_lambda, slong * block_size, const ca_mat_t A, ca_ctx_t ctx)

Computes the block structure of the Jordan canonical form of `A` and
returns 1, or returns 0 on failure. Sets `lambda` to the distinct
eigenvalues, `num_blocks` to the number of Jordan blocks, and for block
$i$, `block_lambda[i]` to the index in `lambda` of its eigenvalue and
`block_size[i]` to its size. The arrays must have room for $n$ entries.
The order of the blocks is unspecified. Wrapper of
`gr_mat_jordan_blocks`.
:::

::: function
void ca_mat_set_jordan_blocks(ca_mat_t mat, const ca_vec_t lambda, slong num_blocks, slong * block_lambda, slong * block_size, ca_ctx_t ctx)

Sets `mat` to the block diagonal matrix of the Jordan blocks described
by `lambda`, `num_blocks`, `block_lambda` and `block_size` (see
`ca_mat_jordan_blocks`).
:::

::: function
int ca_mat_jordan_transformation(ca_mat_t mat, const ca_vec_t lambda, slong num_blocks, slong * block_lambda, slong * block_size, const ca_mat_t A, ca_ctx_t ctx)

Given the Jordan block structure of `A`, computes a matrix $P$ with $A
= P J P^{-1}$, where $J$ is the Jordan form, writes it to `mat` and
returns 1; returns 0 on failure. Wrapper of
`gr_mat_jordan_transformation`.
:::

::: function
int ca_mat_jordan_form(ca_mat_t J, ca_mat_t P, const ca_mat_t A, ca_ctx_t ctx)

Computes the Jordan decomposition $A = P J P^{-1}$ and returns 1, or
returns 0 on failure. `P` may be `NULL` if only $J$ is wanted. This
combines `ca_mat_jordan_blocks`, `ca_mat_set_jordan_blocks` and
`ca_mat_jordan_transformation`, which give direct access to the
eigenvalues and the block structure. Wrapper of `gr_mat_jordan_form`.
:::

## Matrix functions {#functions}

::: function
int ca_mat_exp(ca_mat_t res, const ca_mat_t A, ca_ctx_t ctx)

Sets `res` to the matrix exponential $e^A$ and returns 1, or returns 0
if the Jordan decomposition of $A$, on which the computation is based,
cannot be computed. Wrapper of `gr_mat_exp`.
:::

::: function
truth_t ca_mat_log(ca_mat_t res, const ca_mat_t A, ca_ctx_t ctx)

Sets `res` to a matrix logarithm of `A` and returns `T_TRUE`; returns
`T_FALSE` if $A$ is proved to have no logarithm (it is singular), and
`T_UNKNOWN` if this cannot be decided or the Jordan decomposition fails.
The branch is defined by the principal logarithms of the eigenvalues.
Wrapper of `gr_mat_log`.
:::

## Internal functions {#internal}

::: function
ca_field_ptr _ca_mat_same_field(const ca_mat_t A, ca_ctx_t ctx)
ca_field_ptr _ca_mat_same_field2(const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)

If all entries of `A` (and `B`) are rational numbers or elements of one
common field $K$, returns $K$ (the field $\mathbb{Q}$ if all entries are
rational); returns `NULL` if the entries lie in different fields or one
is a special value.
:::
