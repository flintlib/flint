---
title: fmpz_lll.h
header: fmpz_lll.h
description: LLL reduction of integer lattices given by a basis or a Gram matrix, with floating-point Gram–Schmidt data (L² algorithm), certified reducedness tests and removals of long vectors.
category: Matrices
order: 2
---

This module reduces integer lattices with the LLL algorithm. A lattice is
given either by a basis, the rows of an `fmpz_mat_t`, or by its Gram
matrix, and is reduced in place; the row operations can also be applied
to a second matrix to obtain the unimodular transformation. The
implementation is the floating-point L² algorithm of Nguyen and Stehlé,
which computes the Gram–Schmidt data in `double` precision when possible
and in multiprecision floating point otherwise, checks the result with a
certified reducedness test, and for large entries first reduces truncated
versions of the basis. Variants with *removals* discard vectors whose
Gram–Schmidt norms exceed a bound, as needed by van Hoeij's polynomial
factoring algorithm.

Most users only need `fmpz_lll_context_init_default` (or
`fmpz_lll_context_init`) and `fmpz_lll`; see [](fmpz_lll/lll). The other
functions expose the individual strategies and their building blocks. The
exact (rational arithmetic) LLL algorithms and reducedness tests are in
[](fmpz_mat#lll).

::: summary algorithms/linear-algebra/lll :::

## Types and parameters {#types}

::: type
fmpz_lll_struct
fmpz_lll_t

An LLL context: the fields `delta` and `eta` (the parameters $\delta$ and
$\eta$, as `double`), `rt` (the representation of the input, a
`rep_type`) and `gt` (the kind of Gram matrix used during the reduction,
a `gram_type`). An `fmpz_lll_t` is an array of one `fmpz_lll_struct`. A
context allocates no memory and needs no clearing.
:::

::: type
rep_type

The representation of the lattice: `Z_BASIS` if the rows of the matrix
are a basis (or generating set) of the lattice, `GRAM` if the matrix is
the Gram matrix $BB^T$ of a basis $B$, which must then be symmetric
positive semidefinite.
:::

::: type
gram_type

For `Z_BASIS` input, how the Gram matrix entries $\langle b_i, b_j
\rangle$ are obtained during the reduction: `APPROX` computes them in
floating point from approximations of the basis vectors, `EXACT` keeps
them exactly as integers. Ignored for `GRAM` input.
:::

::: type
fmpz_gram_union
fmpz_gram_t

A union holding the Gram matrix used during a reduction: `appSP` (a
`d_mat_t` of approximations), `appSP2` (a `gr_mat_t` of multiprecision
approximations) or `exactSP` (an exact `fmpz_mat_t`). Used by the Babai
procedures below.
:::

The input must be consistent with `rt` and `gt`; this is not checked. The
parameters must satisfy $0.25 < \delta < 1$ and
$0.5 \le \eta < \sqrt{\delta}$ (for floating-point reduction $\eta > 0.5$
is needed in practice). The output is $(\delta, \eta)$-LLL-reduced in the
sense of [](algorithms/linear-algebra/lll#def-lll).

::: macro
CPU_SIZE_1
FMPZ_LLL_MAX_LONG

The number of mantissa bits of a `double` used by the reduction (53 on
64-bit machines, 31 on 32-bit machines), and the largest integer
$2^{\mathrm{CPU\_SIZE\_1}} - 1$ handled as an exact `double`.
:::

::: macro
SIZE_RED_FAILURE_THRESH

After 20 iterations of a size reduction, the number of bits (5) by which
the largest Gram–Schmidt coefficient must have decreased; otherwise the
size reduction is deemed to fail for lack of precision.
:::

## Context initialisation {#context}

::: function
void fmpz_lll_context_init_default(fmpz_lll_t fl)

Sets `fl` to the default parameters $\delta = 0.99$, $\eta = 0.51$, with
representation `Z_BASIS` and Gram type `APPROX`.
:::

::: function
void fmpz_lll_context_init(fmpz_lll_t fl, double delta, double eta, rep_type rt, gram_type gt)

Sets the fields of `fl` to the given values. Nothing is checked.
:::

::: function
void fmpz_lll_randtest(fmpz_lll_t fl, flint_rand_t state)

Sets `fl` to random valid parameters, $\delta \in (0.25, 1)$ and
$\eta \in (0.5, \sqrt{\delta})$, and a random representation and Gram
type. For testing.
:::

## Main functions {#main}

In all reduction functions, $B$ is reduced in place. $U$ may be `NULL`;
otherwise it must have as many rows as $B$, and every row operation
applied to $B$ is also applied to $U$, so that with $U$ initialised to
the identity matrix ($d \times d$, or $d \times k$ to track only some
columns) it ends as the transformation $U$ with
$B_{\text{new}} = U B_{\text{old}}$ (for `Z_BASIS` input). An exception is
raised if $U$ has the wrong number of rows.

::: function
void fmpz_lll(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl)

Reduces $B$ according to `fl`. For `Z_BASIS` input the rows of $B$ may be
linearly dependent; the zero vectors of the result then come first. See
[](fmpz_lll/lll).

::: details
Calls `fmpz_lll_with_removal_ulll` with truncation to 250 bits and no
removal bound: entries larger than that are first reduced in truncated
form, and the final reduction is done by `fmpz_lll_wrapper_with_removal_knapsack`,
which tries doubles, heuristic dot products and multiprecision floats in
turn until the result is certified reduced.
:::
:::

::: function
int fmpz_lll_with_removal(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl)

Reduces $B$ like `fmpz_lll` and returns the number `newd` of leading
vectors to keep: the trailing vectors whose squared Gram–Schmidt norm
exceeds the bound `gs_B` are not counted. If `gs_B` is `NULL`, no
removal is done and the function returns 0.

::: details
The first `newd` rows are $(\delta, \eta)$-reduced and the squared
Gram–Schmidt norms of the others exceed `gs_B`; the floating-point
variants compare an estimate of half the squared norm with the bound.
:::
:::

## Reduction strategies {#strategies}

These functions implement the floating-point L² algorithm with different
arithmetic. The `double` versions are heuristic: they may return $-1$ or a
basis that is not reduced, which is why the wrappers check the result.
For the functions without removal, the return value is 0 on success and
$-1$ on failure. For the functions with removal, the return value is the
number of vectors kept (or 0 if `gs_B` is `NULL`), or $-1$ on failure.

::: function
int fmpz_lll_d(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl)
int fmpz_lll_d_with_removal(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl)

L² with Gram–Schmidt data in `double` precision, the basis vectors being
approximated by `double` vectors scaled by powers of two. Uses the fast
Babai procedure (`fmpz_lll_check_babai`), switching to the heuristic one
for one iteration after 50 failures. Fails, returning $-1$, if the
precision is insufficient or an iteration limit is reached.
:::

::: function
int fmpz_lll_d_heuristic(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl)
int fmpz_lll_d_heuristic_with_removal(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl)

As `fmpz_lll_d`, but always computes the approximate dot products with
`fmpz_lll_heuristic_dot`, which recomputes a dot product exactly when it
detects cancellation.
:::

::: function
int fmpz_lll_d_with_removal_knapsack(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl)

As `fmpz_lll_d_with_removal`, but size-reduces each new vector in advance
when it is first reached, which is faster for knapsack-type lattices,
and tolerates 150 instead of 50 failures of the fast Babai procedure.
:::

::: function
int fmpz_lll_mpf2(fmpz_mat_t B, fmpz_mat_t U, flint_bitcnt_t prec, const fmpz_lll_t fl)
int fmpz_lll_mpf2_with_removal(fmpz_mat_t B, fmpz_mat_t U, flint_bitcnt_t prec, const fmpz_t gs_B, const fmpz_lll_t fl)

L² with Gram–Schmidt data in floating point of `prec` bits. Returns $-1$
if the precision turns out to be insufficient.

::: details
The arithmetic is that of `nfloat` if `nfloat_ctx_init` accepts the
precision, and `arf` (`gr_ctx_init_real_float_arf`) otherwise. When the
entries fit in a few limbs, the basis is packed (`fmpz_lll_packed_t`) and
the Gram matrix entries are computed exactly and then rounded.
:::
:::

::: function
int fmpz_lll_mpf(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl)
int fmpz_lll_mpf_with_removal(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl)

Calls `fmpz_lll_mpf2` (respectively the version with removal) with
precisions $64, 128, 192, \ldots$ bits for the first 20 attempts and then
doubling the precision, until the result is reduced according to
`fmpz_lll_is_reduced` (respectively `fmpz_lll_is_reduced_with_removal`).
Returns the value of the last call. Since a precision linear in the
dimension suffices [@NS2009], this terminates.
:::

::: function
int fmpz_lll_wrapper(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl)
int fmpz_lll_wrapper_with_removal(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl)
int fmpz_lll_wrapper_with_removal_knapsack(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl)

Reduces $B$ with the cheapest strategy that succeeds: `fmpz_lll_d` (or
`fmpz_lll_d_with_removal`, `fmpz_lll_d_with_removal_knapsack`); if it
fails or the result is not reduced, `fmpz_lll_d_heuristic` (or its
removal version), but only for `Z_BASIS` input with Gram type `APPROX`;
and finally `fmpz_lll_mpf` (or `fmpz_lll_mpf_with_removal`). Returns the
value of the last strategy run.

::: details
Reducedness is checked with `fmpz_lll_is_reduced` at 53 bits
(respectively `fmpz_lll_is_reduced_with_removal` at 120 bits), which is
always conclusive, so the result of the wrappers is reduced.
:::
:::

## Truncation {#ulll}

::: function
int fmpz_lll_with_removal_ulll(fmpz_mat_t FM, fmpz_mat_t UM, slong new_size, const fmpz_t gs_B, const fmpz_lll_t fl)

Reduces `FM` with the truncation strategy: while the entries have more
than `new_size` bits, the lattice $(I \mid \lfloor \mathit{FM}/2^e
\rfloor)$, with $e$ chosen to leave `new_size` bits, is reduced with
`fmpz_lll_wrapper_with_removal_knapsack` and its transformation applied to
`FM`; then `FM` itself is reduced. Returns the number of vectors kept, as
`fmpz_lll_with_removal`. `UM` captures the transformation as in the other
functions. For `GRAM` input no truncation is done.

::: details
The truncation is repeated while it is effective: until the
transformation found is the identity or the entries shrink by less than
`new_size`/4 bits. See [](algorithms/linear-algebra/lll#ulll).
:::
:::

::: function
void fmpz_lll_storjohann_ulll(fmpz_mat_t FM, slong new_size, const fmpz_lll_t fl)

The truncation strategy of `fmpz_lll_with_removal_ulll` with
`fmpz_mat_lll_storjohann` as the reduction, without removals and without
transformation matrix. For `GRAM` input it calls
`fmpz_lll_with_removal_ulll` without removal bound.

::: details
This function is not covered by FLINT's test suite.
:::
:::

## Reducedness tests {#is-reduced}

::: function
int fmpz_lll_is_reduced_d(const fmpz_mat_t B, const fmpz_lll_t fl)
int fmpz_lll_is_reduced_mpfr(const fmpz_mat_t B, const fmpz_lll_t fl, flint_bitcnt_t prec)
int fmpz_lll_is_reduced_d_with_removal(const fmpz_mat_t B, const fmpz_lll_t fl, const fmpz_t gs_B, int newd)
int fmpz_lll_is_reduced_mpfr_with_removal(const fmpz_mat_t B, const fmpz_lll_t fl, const fmpz_t gs_B, int newd, flint_bitcnt_t prec)

Returns 1 if $B$ is proven to be reduced (in the sense of
`fmpz_mat_is_reduced`, `fmpz_mat_is_reduced_gram` or their removal
versions, according to `fl->rt`), and 0 if the test is inconclusive. The
`_d` versions work in `double` arithmetic, the `_mpfr` versions with
`nfloat` numbers of `prec` bits; both use directed rounding to obtain
rigorous error bounds.

::: details
For a basis, these follow Villard's certification of the QR factor
[@Villard2007]: an approximate QR decomposition and an approximate inverse
of $R$, with rigorous bounds on their errors computed with upward and
downward rounding, give intervals for the Gram–Schmidt coefficients on
which the conditions are checked. For a Gram matrix, an approximate
Cholesky decomposition is checked in the same way. The `_d` versions
require `fenv.h` and return 0 without it; they also return 0 if an entry
exceeds the range of `double`. The removal versions take `gs_B` and
`newd` as in `fmpz_mat_is_reduced_with_removal`.
:::
:::

::: function
int fmpz_lll_is_reduced(const fmpz_mat_t B, const fmpz_lll_t fl, flint_bitcnt_t prec)
int fmpz_lll_is_reduced_with_removal(const fmpz_mat_t B, const fmpz_lll_t fl, const fmpz_t gs_B, int newd, flint_bitcnt_t prec)

Returns 1 if $B$ is reduced and 0 otherwise, conclusively. Leading zero
rows of a basis are ignored. The fast tests `fmpz_lll_is_reduced_d` and
`fmpz_lll_is_reduced_mpfr` (with precision `prec`) are tried first; if
both are inconclusive, the exact test `fmpz_mat_is_reduced` (or
`fmpz_mat_is_reduced_gram`, or the removal versions) decides. If `gs_B`
is `NULL`, `fmpz_lll_is_reduced_with_removal` is `fmpz_lll_is_reduced`.
:::

## Babai procedures {#babai}

These functions perform the size reduction of one basis vector, the inner
step of the L² algorithm, and maintain the Gram–Schmidt data. They are
building blocks of the reduction functions above.

In each of them, row `kappa` of `B` is size-reduced against the previous
rows; `U` (if not `NULL`) receives the same operations; `mu` and `r` hold
the Gram–Schmidt coefficients $\mu_{ij}$ and $r_{ij}$; on return `s[i]`
is the squared norm of the projection of row `kappa` orthogonally to
rows $0, \ldots, i - 1$; `appB` holds approximations of the rows scaled
by $2^{-\mathit{expo}[i]}$; `A` is the Gram matrix (approximate or exact
according to `fl`); `a` is the first index to reduce against; `zeros` is
the number of leading zero vectors; `kappamax` the largest index reached
so far; `n` the number of columns to consider. They return 0 on success
and $-1$ if the size reduction does not converge (insufficient
precision).

::: function
int fmpz_lll_check_babai(int kappa, fmpz_mat_t B, fmpz_mat_t U, d_mat_t mu, d_mat_t r, double * s, d_mat_t appB, int * expo, fmpz_gram_t A, int a, int zeros, int kappamax, int n, const fmpz_lll_t fl)
int fmpz_lll_check_babai_heuristic_d(int kappa, fmpz_mat_t B, fmpz_mat_t U, d_mat_t mu, d_mat_t r, double * s, d_mat_t appB, int * expo, fmpz_gram_t A, int a, int zeros, int kappamax, int n, const fmpz_lll_t fl)

Size reduction in `double` precision, with plain floating-point dot
products, respectively with `fmpz_lll_heuristic_dot`.
:::

::: function
int fmpz_lll_advance_check_babai(int cur_kappa, int kappa, fmpz_mat_t B, fmpz_mat_t U, d_mat_t mu, d_mat_t r, double * s, d_mat_t appB, int * expo, fmpz_gram_t A, int a, int zeros, int kappamax, int n, const fmpz_lll_t fl)
int fmpz_lll_advance_check_babai_heuristic_d(int cur_kappa, int kappa, fmpz_mat_t B, fmpz_mat_t U, d_mat_t mu, d_mat_t r, double * s, d_mat_t appB, int * expo, fmpz_gram_t A, int a, int zeros, int kappamax, int n, const fmpz_lll_t fl)

As `fmpz_lll_check_babai` and `fmpz_lll_check_babai_heuristic_d`, but
reduce row `kappa` only against the rows before `cur_kappa`, which are
already LLL-reduced. Used for the early size reductions of the knapsack
variant.
:::

::: function
int fmpz_lll_check_babai_heuristic(int kappa, fmpz_mat_t B, fmpz_mat_t U, gr_mat_t mu, gr_mat_t r, gr_ptr s, gr_mat_t appB, fmpz_gram_t A, int a, int zeros, int kappamax, int n, gr_ptr tmp, gr_ptr rtmp, gr_ctx_t ctx, const fmpz_lll_t fl)

Size reduction with Gram–Schmidt data in the floating-point type of the
context `ctx` (`nfloat` or `arf`), as used by `fmpz_lll_mpf2`; `s` is a
vector and `tmp`, `rtmp` are scratch variables, all of that type.
:::

::: function
double fmpz_lll_heuristic_dot(const double * vec1, const double * vec2, slong len2, const fmpz_mat_t B, slong k, slong j, slong exp_adj)

Returns the dot product of the first `len2` entries of the `double`
vectors `vec1` and `vec2`, which approximate rows $k$ and $j$ of $B$
scaled by powers of two. If cancellation is detected (the square of the
result is at most $2^{-70}$ times the product of the squared norms), the
dot product of the rows of $B$ is computed exactly instead and returned
scaled by $2^{-\mathit{exp\_adj}}$, where `exp_adj` is the sum of the
scaling exponents of the two rows.
:::

::: function
int fmpz_lll_shift(const fmpz_mat_t B)

Returns the largest distance, over the rows $i$ of $B$, between column
$i$ and the last nonzero entry of row $i$, computed incrementally. The
reduction functions only consider the first
$\min(\mathit{kappamax} + 1 + \mathit{shift}, n)$ columns when
computing dot products, which saves work for bases in triangular form.
:::

## Packed bases {#packed}

The floating-point reduction functions work on a packed copy of the
basis when its entries are small: every entry occupies the same number
$m$ of limbs in two's complement, and the entries of a row are
contiguous. Row operations are done limb-wise modulo $2^{m\beta}$, which
is exact because each row carries an upper bound on the bit length of its
entries, maintained by all operations; the number of limbs grows when a
bound would exceed $m\beta - 2$ bits and shrinks when all rows fit in
fewer limbs. This is much cheaper than operating on `fmpz` entries.

::: type
fmpz_lll_packed_struct
fmpz_lll_packed_t

A packed basis: `entries` ($d n m$ limbs), `rows` (pointers to the rows,
permuted by row moves), `bits` (a bound on the bit length of the entries
of each row), and the numbers `d` of rows, `n` of columns and `m` of
limbs per entry.
:::

::: macro
FMPZ_LLL_PACKED_MAX_LIMBS
FMPZ_LLL_PACKED_MAX_LIMBS_MPF
FMPZ_LLL_PACKED_MAX_SIZE

Packing is used if the entries need at most 32 limbs (in `fmpz_lll_d`
and its variants), respectively at most `FLINT_MPN_DOT_TAB_N` limbs (in
`fmpz_lll_mpf2`, where Gram entries are computed exactly from the packed
rows), and if the packed matrix has at most $2^{26}$ limbs.
:::

::: function
slong fmpz_lll_packed_limbs(const fmpz_mat_t B)

Returns the number of limbs per entry for packing $B$: enough for its
largest entry plus 18 bits of headroom.
:::

::: function
void fmpz_lll_packed_init(fmpz_lll_packed_t P, slong d, slong n, slong m)
void fmpz_lll_packed_clear(fmpz_lll_packed_t P)

Initialises $P$ as a packed $d \times n$ matrix with $m$ limbs per entry,
respectively frees it.
:::

::: function
void fmpz_lll_packed_set_fmpz_mat(fmpz_lll_packed_t P, const fmpz_mat_t B)
void fmpz_lll_packed_get_fmpz_mat(fmpz_mat_t B, const fmpz_lll_packed_t P)

Packs $B$, whose entries must fit in $m\beta - 2$ bits, into $P$,
respectively unpacks $P$ into $B$ (rows in their current order).
:::

::: function
void fmpz_lll_packed_tighten(fmpz_lll_packed_t P, slong i)
void fmpz_lll_packed_grow(fmpz_lll_packed_t P)
void fmpz_lll_packed_maybe_shrink(fmpz_lll_packed_t P)

Recomputes the exact bit bound of row $i$; increases the number of limbs
per entry by one; decreases it by one if all rows fit comfortably in one
limb less.
:::

::: function
void fmpz_lll_packed_move_row(fmpz_lll_packed_t P, slong i, slong j)

Moves row $i$ to position $j$, shifting the rows in between (by
permuting row pointers).
:::

::: function
void fmpz_lll_packed_row_sub(fmpz_lll_packed_t P, slong i, slong j)
void fmpz_lll_packed_row_add(fmpz_lll_packed_t P, slong i, slong j)
void fmpz_lll_packed_row_submul_si(fmpz_lll_packed_t P, slong i, slong j, slong x)
void fmpz_lll_packed_row_submul_fmpz(fmpz_lll_packed_t P, slong i, slong j, const fmpz_t x)
void fmpz_lll_packed_row_submul_si_2exp(fmpz_lll_packed_t P, slong i, slong j, slong x, ulong e)

Sets row $i$ to row $i$ minus row $j$, plus row $j$, minus $x$ times row
$j$, respectively minus $x 2^e$ times row $j$, exactly, growing the
representation if necessary.
:::

::: function
slong fmpz_lll_packed_get_d_vec_2exp(double * appv, fmpz_lll_packed_t P, slong i)

Sets `appv` to `double` approximations of the entries of row $i$ scaled
by $2^{-e}$ and returns $e$, as `_fmpz_vec_get_d_vec_2exp`; also sets the
bit bound of the row to its exact value.
:::

::: function
void fmpz_lll_packed_dot(fmpz_t res, const fmpz_lll_packed_t P, slong i, slong j, slong len)

Sets `res` to the exact dot product of the first `len` entries of rows
$i$ and $j$.
:::

::: function
double fmpz_lll_packed_heuristic_dot(const double * vec1, const double * vec2, slong len2, const fmpz_lll_packed_t P, slong k, slong j, slong exp_adj)

`fmpz_lll_heuristic_dot` for a packed basis.
:::

::: function
int fmpz_lll_check_babai_packed(int kappa, fmpz_lll_packed_t P, fmpz_mat_t U, d_mat_t mu, d_mat_t r, double * s, d_mat_t appB, int * expo, fmpz_gram_t A, int a, int zeros, int kappamax, int n, const fmpz_lll_t fl, int heuristic)
int fmpz_lll_advance_check_babai_packed(int cur_kappa, int kappa, fmpz_lll_packed_t P, fmpz_mat_t U, d_mat_t mu, d_mat_t r, double * s, d_mat_t appB, int * expo, fmpz_gram_t A, int a, int zeros, int kappamax, int n, const fmpz_lll_t fl, int heuristic)
int fmpz_lll_check_babai_heuristic_packed(int kappa, fmpz_lll_packed_t P, fmpz_mat_t U, gr_mat_t mu, gr_mat_t r, gr_ptr s, fmpz_gram_t A, int a, int zeros, int kappamax, int n, gr_ptr tmp, gr_ptr rtmp, gr_ctx_t ctx, const fmpz_lll_t fl)

The Babai procedures `fmpz_lll_check_babai` (or
`fmpz_lll_check_babai_heuristic_d` if `heuristic` is nonzero),
`fmpz_lll_advance_check_babai` and `fmpz_lll_check_babai_heuristic`
operating on a packed basis. In the multiprecision version the Gram
matrix entries are computed exactly.
:::
