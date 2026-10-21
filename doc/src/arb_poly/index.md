---
title: arb_poly.h
header: arb_poly.h
description: Dense univariate polynomials and truncated power series with real ball coefficients, with rigorous multiplication, division, composition, evaluation, elementary and special functions of power series, and root refinement.
category: Real and complex numbers
order: 6
---

This module implements dense univariate polynomials whose coefficients
are real balls (`arb_t`). An `arb_poly_t` with coefficient balls $[m_k
\pm r_k]$ represents the set of all real polynomials $\sum_k c_k x^k$
with $|c_k - m_k| \le r_k$, and every operation returns a ball
polynomial that contains the exact result for every choice of inputs in
the input balls. The same objects serve as truncated power series: the
functions with `_series` in their name compute a result modulo $x^n$ for
a given length $n$, which makes the module the basis of the Taylor
expansions (automatic differentiation) of special functions in FLINT.

Products are computed either by the classical algorithm with dot
products, or by splitting the polynomials into blocks of similarly
scaled coefficients and multiplying each pair of blocks exactly over
$\Z$, with the radii propagated by separate products of nonnegative
polynomials ([](algorithms/ball-arithmetic/polynomial-multiplication)).
Division, power series inversion, the exponential and the other
elementary functions of power series use Newton iteration and cost
$O(M(n))$ operations. Composition, reversion, multipoint evaluation,
interpolation and many transcendental functions are wrappers around
the generic implementations in [](gr_poly), called with the real ball
context. The complex analogue is [](acb_poly), real root isolation for
integer polynomials is in [](arb_fmpz_poly), and the guide
[](arb_poly/power-series) shows how the module is used for Taylor
expansions.

::: summary algorithms/ball-arithmetic/polynomial-multiplication :::

::: summary algorithms/ball-arithmetic/power-series :::

## Types and conventions {#conventions}

::: type
arb_poly_struct
arb_poly_t

A polynomial with real ball coefficients. The structure has the fields
`coeffs` (an array of `alloc` initialised `arb_struct`), `alloc` and
`length`; the coefficient of $x^i$ is `coeffs[i]`. An `arb_poly_t` is
an array of one `arb_poly_struct`, so it is passed by reference.
:::

- A polynomial is *normalised* if its length is zero or its leading
  coefficient is not the exact zero ball. A leading coefficient such as
  $[0 \pm 10^{-20}]$ is kept, so the length and the degree are only
  upper bounds for those of the exact polynomials represented. Functions
  without a leading underscore return normalised polynomials.
- Every function that rounds takes a precision `prec` in bits, used for
  the midpoint arithmetic as in [](arb). The result is always a valid
  enclosure; the precision only determines how tight it is.
- Functions without a leading underscore manage memory, accept any
  lengths (including zero) and allow every output to be the same object
  as any input, unless stated otherwise. Several outputs of one function
  (such as quotient and remainder) must be distinct objects.
- Functions with a leading underscore work on coefficient arrays `(poly,
  len)`, write a fixed number of output coefficients and leave memory
  management to the caller. Unless stated otherwise they require all
  lengths to be positive and do not allow the output to overlap an
  input.
- Functions described as wrappers of a `_gr_poly` function call it with
  the context `gr_ctx_init_real_arb(ctx, prec)` (see [](gr_domains));
  the generic code works directly on the `arb_struct` arrays. If it
  reports failure (for example because it would have to invert a ball
  that contains zero), the output is set to indeterminate balls
  $[\operatorname{NaN} \pm \infty]$. Conversely, polynomials over the
  generic real ball context use this module for multiplication,
  division, Taylor shifts and root finding.

## Memory management {#memory}

::: function
void arb_poly_init(arb_poly_t poly)
void arb_poly_init2(arb_poly_t poly, slong len)

Initialise `poly` as the zero polynomial. `arb_poly_init2` also
allocates space for `len` coefficients.
:::

::: function
void arb_poly_clear(arb_poly_t poly)

Frees all coefficients and the coefficient array of `poly`.
:::

::: function
void arb_poly_fit_length(arb_poly_t poly, slong len)

Ensures that `poly` has at least `len` allocated and initialised
coefficients. When it has to grow, the allocation is at least doubled.
The length is not changed.
:::

::: function
void _arb_poly_set_length(arb_poly_t poly, slong len)

Sets the length of `poly` to `len`, which must not exceed the
allocation. If the length decreases, the coefficients from `len` up to
the old length are set to zero, so that unused coefficients are always
zero.
:::

::: function
void _arb_poly_normalise(arb_poly_t poly)

Decreases the length of `poly` until the leading coefficient is not the
exact zero ball. A coefficient that merely contains zero is kept.
:::

::: function
void arb_poly_swap(arb_poly_t poly1, arb_poly_t poly2)

Swaps `poly1` and `poly2` efficiently.
:::

::: function
slong arb_poly_allocated_bytes(const arb_poly_t x)

Returns the number of bytes allocated on the heap by `x` (the
coefficient array and the limbs of the coefficients), excluding
`sizeof(arb_poly_struct)`.
:::

## Basic manipulation {#basic}

::: function
slong arb_poly_length(const arb_poly_t poly)
slong arb_poly_degree(const arb_poly_t poly)

Return the length of `poly` and its degree, the length minus one ($-1$
for the zero polynomial).

::: details
If leading coefficients are balls that contain zero, the exact
polynomials represented by `poly` may have smaller degree; the returned
degree is then only an upper bound.
:::
:::

::: function
int arb_poly_is_zero(const arb_poly_t poly)
int arb_poly_is_one(const arb_poly_t poly)
int arb_poly_is_x(const arb_poly_t poly)

Return whether `poly` is exactly the polynomial 0, 1 or $x$ (with exact
coefficients).
:::

::: function
void arb_poly_zero(arb_poly_t poly)
void arb_poly_one(arb_poly_t poly)

Set `poly` to the constant 0 or 1.
:::

::: function
void arb_poly_set(arb_poly_t dest, const arb_poly_t src)
void arb_poly_set_round(arb_poly_t dest, const arb_poly_t src, slong prec)

Set `dest` to a copy of `src`, either exactly or with the midpoints
rounded to `prec` bits (the rounding error is added to the radii).
:::

::: function
void arb_poly_set_trunc(arb_poly_t dest, const arb_poly_t src, slong n)
void arb_poly_set_trunc_round(arb_poly_t dest, const arb_poly_t src, slong n, slong prec)

Set `dest` to `src` truncated to length `n` (reduced modulo $x^n$),
exactly or rounded to `prec` bits.
:::

::: function
void arb_poly_set_coeff_si(arb_poly_t poly, slong n, slong c)
void arb_poly_set_coeff_arb(arb_poly_t poly, slong n, const arb_t c)

Set the coefficient of $x^n$ in `poly` to `c`, extending the
polynomial with zeros if necessary, and normalise. Require $n \ge 0$.
:::

::: function
void arb_poly_get_coeff_arb(arb_t v, const arb_poly_t poly, slong n)

Sets `v` to the coefficient of $x^n$ in `poly`, which is zero if $n$ is
at least the length. Requires $n \ge 0$.
:::

::: macro
arb_poly_get_coeff_ptr(poly, n)

Evaluates to a pointer to the coefficient of $x^n$ in `poly`, or to
`NULL` if $n$ is at least the length. Requires $n \ge 0$.
:::

::: function
void _arb_poly_reverse(arb_ptr res, arb_srcptr poly, slong len, slong n)

Sets `{res, n}` to the reversal of `{poly, len}` as a polynomial of
length $n \ge \mathit{len}$, that is, `res[i] = poly[n - 1 - i]` with
zeros where the index is out of range. Aliasing is allowed. A wrapper
of `_gr_poly_reverse`.
:::

::: function
void _arb_poly_shift_right(arb_ptr res, arb_srcptr poly, slong len, slong n)
void arb_poly_shift_right(arb_poly_t res, const arb_poly_t poly, slong n)

Set `res` to `poly` divided by $x^n$, discarding the coefficients of
$x^0, \ldots, x^{n-1}$. Require $n \ge 0$. The underscore version
writes $\mathit{len} - n$ coefficients, requires $n \le \mathit{len}$
and allows `res` to be `poly`. Wrappers of `_gr_poly_shift_right`.
:::

::: function
void _arb_poly_shift_left(arb_ptr res, arb_srcptr poly, slong len, slong n)
void arb_poly_shift_left(arb_poly_t res, const arb_poly_t poly, slong n)

Set `res` to `poly` multiplied by $x^n$. Require $n \ge 0$. The
underscore version writes $\mathit{len} + n$ coefficients and allows
`res` to be `poly`. Wrappers of `_gr_poly_shift_left`.
:::

::: function
void arb_poly_truncate(arb_poly_t poly, slong n)

Truncates `poly` to length at most `n` (reduces it modulo $x^n$) and
normalises it. Requires $n \ge 0$.
:::

::: function
slong arb_poly_valuation(const arb_poly_t poly)

Returns the index of the first coefficient of `poly` that is not the
exact zero ball, or $-1$ if `poly` is zero. A coefficient that merely
contains zero counts as nonzero, so the result is a lower bound for the
valuation of the exact polynomials represented.
:::

## Conversions {#conversions}

::: function
void arb_poly_set_fmpz_poly(arb_poly_t poly, const fmpz_poly_t src, slong prec)
void arb_poly_set_fmpq_poly(arb_poly_t poly, const fmpq_poly_t src, slong prec)

Set `poly` to the integer or rational polynomial `src`, rounding each
coefficient to `prec` bits. Integers of at most `prec` bits are
converted exactly.
:::

::: function
void arb_poly_set_si(arb_poly_t poly, slong src)
void arb_poly_set_arb(arb_poly_t poly, const arb_t c)

Set `poly` to the constant polynomial `src` or `c`, exactly. The result
has length zero if the constant is exactly zero.
:::

## Input and output {#io}

::: function
void arb_poly_printd(const arb_poly_t poly, slong digits)
void arb_poly_fprintd(FILE * file, const arb_poly_t poly, slong digits)

Print `poly` to standard output or to `file` as a list `[(c0), (c1),
...]` of its coefficients from the constant term upwards, each printed
with `arb_fprintd` to `digits` decimal digits.
:::

## Random generation {#random}

::: function
void arb_poly_randtest(arb_poly_t poly, flint_rand_t state, slong len, slong prec, slong mag_bits)

Sets `poly` to a random polynomial of length at most `len` whose
coefficients are generated by `arb_randtest` or `arb_randtest_precise`
with the parameters `prec` and `mag_bits`. For testing.
:::

## Comparisons {#comparisons}

::: function
int arb_poly_contains(const arb_poly_t poly1, const arb_poly_t poly2)
int arb_poly_contains_fmpz_poly(const arb_poly_t poly1, const fmpz_poly_t poly2)
int arb_poly_contains_fmpq_poly(const arb_poly_t poly1, const fmpq_poly_t poly2)

Return nonzero if every coefficient of `poly2` is contained in the
corresponding coefficient of `poly1`, where coefficients beyond the
length of `poly2` are zero; that is, if the set represented by `poly1`
contains `poly2`.
:::

::: function
int arb_poly_equal(const arb_poly_t A, const arb_poly_t B)

Returns nonzero if `A` and `B` have the same length and identical
coefficient balls (equal midpoints and equal radii). This compares
representations, not the polynomials represented.
:::

::: function
int _arb_poly_overlaps(arb_srcptr poly1, slong len1, arb_srcptr poly2, slong len2)
int arb_poly_overlaps(const arb_poly_t poly1, const arb_poly_t poly2)

Return nonzero if the coefficients of `poly1` and `poly2` overlap
pairwise (coefficients of the longer polynomial beyond the length of the
shorter one must contain zero), that is, if the two sets of polynomials
intersect. The underscore version requires $\mathit{len1} \ge
\mathit{len2}$.
:::

::: function
int arb_poly_get_unique_fmpz_poly(fmpz_poly_t z, const arb_poly_t x)

If every coefficient of `x` contains exactly one integer, sets `z` to
the integer polynomial with these coefficients and returns nonzero.
Otherwise returns zero, and `z` may have been modified.
:::

## Bounds {#bounds}

::: function
void _arb_poly_majorant(arb_ptr res, arb_srcptr poly, slong len, slong prec)
void arb_poly_majorant(arb_poly_t res, const arb_poly_t poly, slong prec)

Set `res` to the polynomial with exact nonnegative coefficients that
are upper bounds, rounded to `prec` bits, for the absolute values of the
coefficients of `poly`. For $|x| \le t$, its value at $t$ bounds the
absolute value at $x$ of every polynomial in `poly`.
:::

## Addition and scalar arithmetic {#arithmetic}

::: function
void _arb_poly_add(arb_ptr C, arb_srcptr A, slong lenA, arb_srcptr B, slong lenB, slong prec)
void _arb_poly_sub(arb_ptr C, arb_srcptr A, slong lenA, arb_srcptr B, slong lenB, slong prec)

Set `{C, max(lenA, lenB)}` to the sum or difference of `{A, lenA}` and
`{B, lenB}`. The output may be the same array as either input, and the
lengths may be zero.
:::

::: function
void arb_poly_add(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, slong prec)
void arb_poly_sub(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, slong prec)
void arb_poly_add_si(arb_poly_t C, const arb_poly_t A, slong B, slong prec)

Set `C` to $A + B$ or $A - B$.
:::

::: function
void arb_poly_add_series(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, slong len, slong prec)
void arb_poly_sub_series(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, slong len, slong prec)

Set `C` to $A + B$ or $A - B$ truncated to length `len`, without
computing the discarded coefficients.
:::

::: function
void arb_poly_neg(arb_poly_t C, const arb_poly_t A)

Sets `C` to $-A$, exactly.
:::

::: function
void arb_poly_scalar_mul_2exp_si(arb_poly_t C, const arb_poly_t A, slong c)

Sets `C` to $2^c A$, exactly.
:::

::: function
void arb_poly_scalar_mul(arb_poly_t C, const arb_poly_t A, const arb_t c, slong prec)
void arb_poly_scalar_mul_si(arb_poly_t C, const arb_poly_t A, slong c, slong prec)
void arb_poly_scalar_div(arb_poly_t C, const arb_poly_t A, const arb_t c, slong prec)

Set `C` to $c A$ or $A / c$. Division by a ball that contains zero gives
non-finite coefficients.
:::

## Multiplication {#multiplication}

::: function
void _arb_poly_mullow_classical(arb_ptr C, arb_srcptr A, slong lenA, arb_srcptr B, slong lenB, slong n, slong prec)
void _arb_poly_mullow_block(arb_ptr C, arb_srcptr A, slong lenA, arb_srcptr B, slong lenB, slong n, slong prec)
void _arb_poly_mullow(arb_ptr C, arb_srcptr A, slong lenA, arb_srcptr B, slong lenB, slong n, slong prec)

Set `{C, n}` to the product of `{A, lenA}` and `{B, lenB}` truncated to
length `n`. Require $\mathit{lenA}, \mathit{lenB} \ge 1$ and $1 \le n
\le \mathit{lenA} + \mathit{lenB} - 1$; the output must not overlap
either input. If `A` and `B` are the same pointer with the same length,
the square is computed.

::: details
The *classical* version computes each output coefficient as one
`arb_dot`. The *block* version scales $x \to 2^{-c} x$, splits each
operand into blocks of coefficients whose midpoints fit in a common
range of about $3\,\mathit{prec}$ bits, multiplies the blocks of
midpoints exactly with `_fmpz_poly_mulmid`, and bounds the propagated
error by products of nonnegative polynomials computed in double or
integer arithmetic. The default `_arb_poly_mullow` uses the classical
version if $n \le 7$ or either length is at most 7, or if $n$ or twice
the shorter length is at most a cutoff (110 for $\mathit{prec} \le 128$,
otherwise $10^4 / (\ln \mathit{prec})^3$ clamped to $[8, 60]$, raised
by 25% for squaring at 256 bits and again at 4096 bits), and the block
version otherwise. See [](arb_poly/mul).
:::
:::

::: function
void arb_poly_mullow_classical(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, slong n, slong prec)
void arb_poly_mullow_block(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, slong n, slong prec)
void arb_poly_mullow(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, slong n, slong prec)

Set `C` to $AB$ truncated to length `n`. If `A` and `B` are the same
object, the square is computed, which is faster and gives tighter
radii than a product of two different polynomials.
:::

::: function
void _arb_poly_mul(arb_ptr C, arb_srcptr A, slong lenA, arb_srcptr B, slong lenB, slong prec)
void arb_poly_mul(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, slong prec)

Set `C` to the full product $AB$. The underscore version writes
$\mathit{lenA} + \mathit{lenB} - 1$ coefficients, requires positive
lengths and no aliasing, and calls `_arb_poly_mullow`. Squaring is
detected as for `arb_poly_mullow`.
:::

::: function
void _arb_poly_mulmid_classical(arb_ptr z, arb_srcptr x, slong xlen, arb_srcptr y, slong ylen, slong nlo, slong nhi, slong prec)
void _arb_poly_mulmid_block(arb_ptr z, arb_srcptr x, slong xlen, arb_srcptr y, slong ylen, slong nlo, slong nhi, slong prec)
void _arb_poly_mulmid(arb_ptr z, arb_srcptr x, slong xlen, arb_srcptr y, slong ylen, slong nlo, slong nhi, slong prec)
void arb_poly_mulmid_classical(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, slong nlo, slong nhi, slong prec)
void arb_poly_mulmid_block(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, slong nlo, slong nhi, slong prec)
void arb_poly_mulmid(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, slong nlo, slong nhi, slong prec)

Set `z` (or `res`) to the coefficients of $x^{\mathit{nlo}}, \ldots,
x^{\mathit{nhi}-1}$ of the product, shifted down to start at $x^0$ (the
*middle product*). The underscore versions write $\mathit{nhi} -
\mathit{nlo}$ coefficients and require positive lengths, $0 \le
\mathit{nlo} < \mathit{nhi} \le \mathit{xlen} + \mathit{ylen} - 1$ and
no aliasing; the algorithms and the choice between them are those of
the `mullow` functions.

::: details
The block version first discards the coefficients of each operand that
can only contribute below $x^{\mathit{nlo}}$. Middle products are used
in the Newton iterations of `_arb_poly_inv_series`,
`_arb_poly_exp_series` and `_arb_poly_sqrt_series`, where the low part
of a product is known in advance
([](algorithms/polynomials/middle-product)).
:::
:::

::: function
void _arb_poly_mul_monic(arb_ptr res, arb_srcptr poly1, slong len1, arb_srcptr poly2, slong len2, slong prec)

Sets `{res, len1 + len2 - 1}` to the product of the monic polynomials
`{poly1, len1}` and `{poly2, len2}`: the leading coefficient is set to
exactly 1, and the others are computed with `_arb_poly_mullow`. No
aliasing.
:::

::: function
void _arb_poly_addmulmid_mag(mag_ptr z, mag_srcptr x, slong xlen, mag_srcptr y, slong ylen, slong nlo, slong nhi, const fmpz_t scale)

Adds to `{z, nhi - nlo}` upper bounds for the coefficients of
$x^{\mathit{nlo}}, \ldots, x^{\mathit{nhi}-1}$ of the product of the
polynomials `{x, xlen}` and `{y, ylen}` with nonnegative `mag_t`
coefficients. The operands are split into blocks after the substitution
$x \to 2^{-\mathit{scale}} x$, as the radii in `_arb_poly_mulmid_block`.
This bounds the radii in `_acb_poly_mulmid_fmpzi`.
:::

## Division {#division}

::: function
void _arb_poly_inv_series(arb_ptr Q, arb_srcptr A, slong Alen, slong len, slong prec)
void arb_poly_inv_series(arb_poly_t Q, const arb_poly_t A, slong n, slong prec)

Set `Q` to the power series inverse $1/A$ truncated to length `n`
(`len`). The result is finite only if the constant term of $A$ excludes
zero; it is indeterminate if $A$ is zero. The underscore version
requires positive lengths and no aliasing.

::: details
The first coefficients are computed by the recurrence $q_i = -q_0
\sum_{j=1}^{i} a_j q_{i-j}$ with `arb_dot` (all of them if $n \le 8$ or
$A$ has length 2, otherwise the first 4), and the rest by the Newton
iteration $Q \leftarrow Q - Q(AQ - 1)$, which doubles the length at each
step and computes only the unknown part of $AQ - 1$ by a middle product.
The recurrence is faster up to much larger lengths but gives worse
error bounds ([](algorithms/ball-arithmetic/power-series)).
:::
:::

::: function
void _arb_poly_div_series(arb_ptr Q, arb_srcptr A, slong Alen, arb_srcptr B, slong Blen, slong n, slong prec)
void arb_poly_div_series(arb_poly_t Q, const arb_poly_t A, const arb_poly_t B, slong n, slong prec)

Set `Q` to the power series quotient $A / B$ truncated to length `n`.
For $n \le 10$ or $B$ of length 2 the coefficients are computed by the
recurrence, otherwise as $A \cdot (1/B)$ with `_arb_poly_inv_series`.
The underscore version requires positive lengths and no aliasing.
:::

::: function
void _arb_poly_div(arb_ptr Q, arb_srcptr A, slong lenA, arb_srcptr B, slong lenB, slong prec)
void _arb_poly_rem(arb_ptr R, arb_srcptr A, slong lenA, arb_srcptr B, slong lenB, slong prec)
void _arb_poly_divrem(arb_ptr Q, arb_ptr R, arb_srcptr A, slong lenA, arb_srcptr B, slong lenB, slong prec)

Compute the quotient `{Q, lenA - lenB + 1}` and/or the remainder `{R,
lenB - 1}` of the division of `{A, lenA}` by `{B, lenB}`. Require
$\mathit{lenA} \ge \mathit{lenB} \ge 1$ and no aliasing. If the leading
coefficient of $B$ contains zero, the output is not finite.
:::

::: function
int arb_poly_divrem(arb_poly_t Q, arb_poly_t R, const arb_poly_t A, const arb_poly_t B, slong prec)

Computes $Q$ and $R$ with $A = BQ + R$ and $\deg R < \deg B$, and
returns nonzero. If $B$ is zero or its leading coefficient contains
zero, returns zero and leaves `Q` and `R` unchanged.

::: details
The quotient is computed by reversing the operands and dividing them as
power series ([](algorithms/polynomials/newton)), and the remainder as
$A - BQ$ from the low $\deg B$ coefficients of the product. The
enclosures contain the exact quotient and remainder for every pair of
polynomials in the input balls (all of which have the degrees of `A`
and `B`, since the leading coefficient of `B` excludes zero).
:::
:::

::: function
void _arb_poly_div_root(arb_ptr Q, arb_t R, arb_srcptr A, slong len, const arb_t c, slong prec)

Divides `{A, len}` by $x - c$ by Horner's scheme, writing the quotient
`{Q, len - 1}` and the remainder $R = A(c)$. Requires $\mathit{len} \ge
1$.
:::

## Composition {#composition}

::: function
void _arb_poly_taylor_shift(arb_ptr g, const arb_t c, slong n, slong prec)
void arb_poly_taylor_shift(arb_poly_t g, const arb_poly_t f, const arb_t c, slong prec)

Set `g` to the Taylor shift $f(x + c)$. The underscore version works in
place on `{g, n}`. Wrappers of `_gr_poly_taylor_shift`.

::: details
Over the real ball context the generic function dispatches to
`_gr_arb_poly_taylor_shift`, which uses Horner's scheme
(`_gr_poly_taylor_shift_horner`) if $n \le 30$, if $n \le 100$ and $c$
has fewer than $\mathit{prec}/100$ bits, or if $n \le 500$, $c$ has a
one-bit midpoint and $n < 30 + 3\sqrt{\mathit{prec}}$; otherwise the
convolution algorithm (`_gr_poly_taylor_shift_convolution`, one
polynomial product) if $\mathit{prec} > 2n$, and the divide-and-conquer
algorithm (`_gr_poly_taylor_shift_divconquer`) if not.
:::
:::

::: function
void _arb_poly_compose(arb_ptr res, arb_srcptr poly1, slong len1, arb_srcptr poly2, slong len2, slong prec)
void arb_poly_compose(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, slong prec)

Set `res` to the composition $f(g(x))$ of $f$ = `poly1` and $g$ =
`poly2`. The underscore version writes $(\mathit{len1} -
1)(\mathit{len2} - 1) + 1$ coefficients and does not allow aliasing.
Wrappers of `_gr_poly_compose` ([](algorithms/polynomials/composition)).
:::

::: function
void _arb_poly_compose_series(arb_ptr res, arb_srcptr poly1, slong len1, arb_srcptr poly2, slong len2, slong n, slong prec)
void arb_poly_compose_series(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, slong n, slong prec)

Set `res` to the power series composition $f(g(x)) \bmod x^n$. The
constant term of $g$ must be exactly zero; `arb_poly_compose_series`
throws an exception otherwise. The underscore version does not allow
aliasing. Wrappers of `_gr_poly_compose_series`.

::: details
The generic function uses Horner's scheme if $f$ or $n$ is short
(length at most 7), divide and conquer if $\mathit{len1}^2 < n$ or the
untruncated composition is shorter than $4n$, and otherwise the
Brent–Kung baby-step giant-step algorithm [@BrentKung1978] for $n$
below $\max(100, 15000/\sqrt{\max(64, \mathit{prec})})$ and the
algorithm of Kinoshita and Li [@KL2024] above
([](algorithms/polynomials/composition)). If some coefficients of the
inputs are not finite, only the coefficients that depend on finite
input are computed; the others are set to indeterminate balls.
:::
:::

::: function
void _arb_poly_revert_series(arb_ptr h, arb_srcptr f, slong flen, slong n, slong prec)
void arb_poly_revert_series(arb_poly_t h, const arb_poly_t f, slong n, slong prec)

Set `h` to the power series reversion of $f$, the expansion of the
compositional inverse $f^{-1}(x)$ modulo $x^n$. The constant term of $f$
must be exactly zero and its linear coefficient must exclude zero
(otherwise the result is indeterminate). The underscore version
requires $\mathit{flen} \ge 2$ and no aliasing. Wrappers of
`_gr_poly_revert_series`, which uses the fast Lagrange inversion of
[@Joh2015b] up to five times the Kinoshita–Li cutoff of the composition,
and Newton iteration with fast composition beyond.
:::

## Evaluation {#evaluation}

::: function
void _arb_poly_evaluate_horner(arb_t y, arb_srcptr f, slong len, const arb_t x, slong prec)
void arb_poly_evaluate_horner(arb_t y, const arb_poly_t f, const arb_t x, slong prec)
void _arb_poly_evaluate_rectangular(arb_t y, arb_srcptr f, slong len, const arb_t x, slong prec)
void arb_poly_evaluate_rectangular(arb_t y, const arb_poly_t f, const arb_t x, slong prec)
void _arb_poly_evaluate(arb_t y, arb_srcptr f, slong len, const arb_t x, slong prec)
void arb_poly_evaluate(arb_t y, const arb_poly_t f, const arb_t x, slong prec)

Set $y = f(x)$ using Horner's rule, rectangular splitting [@PS1973], or
an automatic choice. The underscore versions accept $\mathit{len} = 0$.

::: details
Rectangular splitting computes $x^2, \ldots, x^m$ for $m \approx
\sqrt{\mathit{len}}$ and evaluates the inner polynomials as dot
products, using about $2\sqrt{\mathit{len}}$ full multiplications and
$\mathit{len}$ multiplications by the (short) coefficients. The
automatic choice uses it when $\mathit{prec} \ge 1024$, $\mathit{len}
\ge 5 + 20000 / \mathit{prec}$ and the coefficients have at most
$\mathit{prec}/2$ bits, and Horner's rule otherwise. Horner's rule
usually gives the tighter enclosure when $x$ is a wide ball.
:::
:::

::: function
void _arb_poly_evaluate_acb_horner(acb_t y, arb_srcptr f, slong len, const acb_t x, slong prec)
void arb_poly_evaluate_acb_horner(acb_t y, const arb_poly_t f, const acb_t x, slong prec)
void _arb_poly_evaluate_acb_rectangular(acb_t y, arb_srcptr f, slong len, const acb_t x, slong prec)
void arb_poly_evaluate_acb_rectangular(acb_t y, const arb_poly_t f, const acb_t x, slong prec)
void _arb_poly_evaluate_acb(acb_t y, arb_srcptr f, slong len, const acb_t x, slong prec)
void arb_poly_evaluate_acb(acb_t y, const arb_poly_t f, const acb_t x, slong prec)

Set $y = f(x)$ for a complex ball $x$, using Horner's rule, rectangular
splitting, or the default, which is rectangular splitting.
:::

::: function
void _arb_poly_evaluate2_horner(arb_t y, arb_t z, arb_srcptr f, slong len, const arb_t x, slong prec)
void arb_poly_evaluate2_horner(arb_t y, arb_t z, const arb_poly_t f, const arb_t x, slong prec)
void _arb_poly_evaluate2_rectangular(arb_t y, arb_t z, arb_srcptr f, slong len, const arb_t x, slong prec)
void arb_poly_evaluate2_rectangular(arb_t y, arb_t z, const arb_poly_t f, const arb_t x, slong prec)
void _arb_poly_evaluate2(arb_t y, arb_t z, arb_srcptr f, slong len, const arb_t x, slong prec)
void arb_poly_evaluate2(arb_t y, arb_t z, const arb_poly_t f, const arb_t x, slong prec)

Set $y = f(x)$ and $z = f'(x)$ together, using Horner's rule,
rectangular splitting (which reuses the powers of $x$ for both values)
or the automatic choice of `_arb_poly_evaluate`.
:::

::: function
void _arb_poly_evaluate2_acb_horner(acb_t y, acb_t z, arb_srcptr f, slong len, const acb_t x, slong prec)
void arb_poly_evaluate2_acb_horner(acb_t y, acb_t z, const arb_poly_t f, const acb_t x, slong prec)
void _arb_poly_evaluate2_acb_rectangular(acb_t y, acb_t z, arb_srcptr f, slong len, const acb_t x, slong prec)
void arb_poly_evaluate2_acb_rectangular(acb_t y, acb_t z, const arb_poly_t f, const acb_t x, slong prec)
void _arb_poly_evaluate2_acb(acb_t y, acb_t z, arb_srcptr f, slong len, const acb_t x, slong prec)
void arb_poly_evaluate2_acb(acb_t y, acb_t z, const arb_poly_t f, const acb_t x, slong prec)

Set $y = f(x)$ and $z = f'(x)$ for a complex ball $x$.
:::

## Product trees {#trees}

::: function
void _arb_poly_product_roots(arb_ptr poly, arb_srcptr xs, slong n, slong prec)
void arb_poly_product_roots(arb_poly_t poly, arb_srcptr xs, slong n, slong prec)

Set `poly` to $(x - x_0)(x - x_1) \cdots (x - x_{n-1})$. The
underscore version writes $n + 1$ coefficients. Wrappers of
`_gr_poly_product_roots`, which multiplies the factors in a balanced
tree.
:::

::: function
void _arb_poly_product_roots_complex(arb_ptr poly, arb_srcptr r, slong rn, acb_srcptr c, slong cn, slong prec)
void arb_poly_product_roots_complex(arb_poly_t poly, arb_srcptr r, slong rn, acb_srcptr c, slong cn, slong prec)

Set `poly` to the real polynomial
$$\prod_{i=0}^{\mathit{rn}-1} (x - r_i) \prod_{i=0}^{\mathit{cn}-1} (x - c_i)(x - \bar c_i)$$
with $\mathit{rn}$ real roots and $\mathit{cn}$ pairs of complex
conjugate roots, of which one member of each pair is given. Either count
may be zero. Use `acb_poly_product_roots` for complex roots that are
not grouped in conjugate pairs.
:::

::: function
arb_ptr * _arb_poly_tree_alloc(slong len)
void _arb_poly_tree_free(arb_ptr * tree, slong len)
void _arb_poly_tree_build(arb_ptr * tree, arb_srcptr roots, slong len, slong prec)

Allocate, free and fill a subproduct tree for the factors $x - x_i$ of
`len` points, as used by fast multipoint evaluation and interpolation
([](algorithms/polynomials/subproduct-tree)). Wrappers of
`_gr_poly_tree_alloc`, `_gr_poly_tree_free` and `_gr_poly_tree_build`.
:::

## Multipoint evaluation and interpolation {#multipoint}

::: function
void _arb_poly_evaluate_vec_iter(arb_ptr ys, arb_srcptr poly, slong plen, arb_srcptr xs, slong n, slong prec)
void arb_poly_evaluate_vec_iter(arb_ptr ys, const arb_poly_t poly, arb_srcptr xs, slong n, slong prec)

Set `ys[i]` to the value of `poly` at `xs[i]` for $0 \le i < n$ by
calling `_arb_poly_evaluate` for each point.
:::

::: function
void _arb_poly_evaluate_vec_fast_precomp(arb_ptr vs, arb_srcptr poly, slong plen, arb_ptr * tree, slong len, slong prec)
void _arb_poly_evaluate_vec_fast(arb_ptr ys, arb_srcptr poly, slong plen, arb_srcptr xs, slong n, slong prec)
void arb_poly_evaluate_vec_fast(arb_ptr ys, const arb_poly_t poly, arb_srcptr xs, slong n, slong prec)

Evaluate `poly` at $n$ points by the remainder tree algorithm, given
the points or a precomputed subproduct tree. Wrappers of
`_gr_poly_evaluate_vec_fast_precomp`. The remainder tree is
numerically unstable: unless the precision grows with the number of
points, the radii become large.
:::

::: function
void _arb_poly_interpolate_newton(arb_ptr poly, arb_srcptr xs, arb_srcptr ys, slong n, slong prec)
void arb_poly_interpolate_newton(arb_poly_t poly, arb_srcptr xs, arb_srcptr ys, slong n, slong prec)

Set `poly` to the polynomial of length at most $n$ with
$\mathit{poly}(x_i) = y_i$, by computing divided differences (the
Newton form) and converting to the monomial basis. If two of the points
overlap as balls, the result is not finite.
:::

::: function
void _arb_poly_interpolate_barycentric(arb_ptr poly, arb_srcptr xs, arb_srcptr ys, slong n, slong prec)
void arb_poly_interpolate_barycentric(arb_poly_t poly, arb_srcptr xs, arb_srcptr ys, slong n, slong prec)

Interpolate as `arb_poly_interpolate_newton`, using the barycentric
form of the Lagrange interpolation formula.
:::

::: function
void _arb_poly_interpolation_weights(arb_ptr w, arb_ptr * tree, slong len, slong prec)
void _arb_poly_interpolate_fast_precomp(arb_ptr poly, arb_srcptr ys, arb_ptr * tree, arb_srcptr weights, slong len, slong prec)
void _arb_poly_interpolate_fast(arb_ptr poly, arb_srcptr xs, arb_srcptr ys, slong len, slong prec)
void arb_poly_interpolate_fast(arb_poly_t poly, arb_srcptr xs, arb_srcptr ys, slong n, slong prec)

Interpolate by the fast Lagrange formula on a subproduct tree.
`_arb_poly_interpolation_weights` computes the weights $w_i = 1 /
\prod_{j \ne i} (x_i - x_j)$ from a tree, and the *precomp* version
takes the tree and the weights as input. Wrappers of
`_gr_poly_interpolation_weights` and
`_gr_poly_interpolate_fast_precomp`.
:::

## Differentiation and integration {#calculus}

::: function
void _arb_poly_derivative(arb_ptr res, arb_srcptr poly, slong len, slong prec)
void arb_poly_derivative(arb_poly_t res, const arb_poly_t poly, slong prec)

Set `res` to the derivative of `poly`. The underscore version writes
$\mathit{len} - 1$ coefficients and allows aliasing. Wrappers of
`_gr_poly_derivative`.
:::

::: function
void _arb_poly_nth_derivative(arb_ptr res, arb_srcptr poly, ulong n, slong len, slong prec)
void arb_poly_nth_derivative(arb_poly_t res, const arb_poly_t poly, ulong n, slong prec)

Set `res` to the $n$-th derivative of `poly`. The underscore version
writes $\mathit{len} - n$ coefficients (nothing if $\mathit{len} \le
n$) and allows aliasing. Wrappers of `_gr_poly_nth_derivative`.
:::

::: function
void _arb_poly_integral(arb_ptr res, arb_srcptr poly, slong len, slong prec)
void arb_poly_integral(arb_poly_t res, const arb_poly_t poly, slong prec)

Set `res` to the integral of `poly` with zero constant term. The
underscore version reads `{poly, len - 1}`, writes `{res, len}` and
allows aliasing. Wrappers of `_gr_poly_integral`.
:::

## Transforms {#transforms}

::: function
void _arb_poly_borel_transform(arb_ptr res, arb_srcptr poly, slong len, slong prec)
void arb_poly_borel_transform(arb_poly_t res, const arb_poly_t poly, slong prec)
void _arb_poly_inv_borel_transform(arb_ptr res, arb_srcptr poly, slong len, slong prec)
void arb_poly_inv_borel_transform(arb_poly_t res, const arb_poly_t poly, slong prec)

Map $\sum_k a_k x^k$ to $\sum_k (a_k / k!) x^k$ (Borel transform) or to
$\sum_k a_k k!\, x^k$ (inverse Borel transform). The underscore
versions allow aliasing.
:::

::: function
void _arb_poly_binomial_transform_basecase(arb_ptr b, arb_srcptr a, slong alen, slong len, slong prec)
void arb_poly_binomial_transform_basecase(arb_poly_t b, const arb_poly_t a, slong len, slong prec)
void _arb_poly_binomial_transform_convolution(arb_ptr b, arb_srcptr a, slong alen, slong len, slong prec)
void arb_poly_binomial_transform_convolution(arb_poly_t b, const arb_poly_t a, slong len, slong prec)
void _arb_poly_binomial_transform(arb_ptr b, arb_srcptr a, slong alen, slong len, slong prec)
void arb_poly_binomial_transform(arb_poly_t b, const arb_poly_t a, slong len, slong prec)

Set `b` to the binomial transform of `a` truncated to length `len`,
$b_n = \sum_{k=0}^{n} (-1)^k \binom{n}{k} a_k$. The transform is the
power series composition $f(x) \mapsto (1-x)^{-1} f(x/(x-1))$ and is an
involution. The underscore versions require positive lengths and do not
allow aliasing.

::: details
The *basecase* version evaluates the sums directly, generating the
binomial coefficients by recurrence. The *convolution* version uses
$T(f(x)) = B^{-1}(e^x B(f(-x)))$, where $T$ is the binomial transform
and $B$ the Borel transform, at the cost of one polynomial
multiplication. The default uses the basecase version if `alen` or
`len` is below 10.
:::
:::

::: function
void _arb_poly_graeffe_transform(arb_ptr b, arb_srcptr a, slong len, slong prec)
void arb_poly_graeffe_transform(arb_poly_t b, const arb_poly_t a, slong prec)

Set `b` to the Graeffe transform $G$ of $P$ = `a`, defined by $G(x^2) =
(-1)^d P(x) P(-x)$ with $d = \deg P$, so that the roots of $G$ are the
squares of the roots of $P$. The underscore version reads `{a, len}`,
writes `len` coefficients and allows aliasing.
:::

## Powers and elementary functions {#elementary}

The power series functions in this section compute their result modulo
$x^n$ for a given length $n$ (or `len`). Unless stated otherwise, their
underscore versions require both lengths to be positive and allow the
input to be shorter or longer than the output. The output is finite only
where the function is analytic at the constant term of the input: for
example, `arb_poly_log_series` needs a constant term that excludes zero.
See [](algorithms/ball-arithmetic/power-series) for the algorithms.

::: function
void _arb_poly_pow_ui_trunc_binexp(arb_ptr res, arb_srcptr f, slong flen, ulong exp, slong len, slong prec)
void arb_poly_pow_ui_trunc_binexp(arb_poly_t res, const arb_poly_t poly, ulong exp, slong len, slong prec)

Set `res` to $f^{e}$ truncated to length `len`, by binary
exponentiation with truncated products. The underscore version requires
`len` to be at most the length of the untruncated power, positive
lengths and no aliasing. Wrappers of `_gr_poly_pow_series_ui_binexp`.
:::

::: function
void _arb_poly_pow_ui(arb_ptr res, arb_srcptr f, slong flen, ulong exp, slong prec)
void arb_poly_pow_ui(arb_poly_t res, const arb_poly_t poly, ulong exp, slong prec)

Set `res` to $f^{e}$. The underscore version writes $e(\mathit{flen} -
1) + 1$ coefficients, requires $\mathit{flen} \ge 1$ and no aliasing.
:::

::: function
void _arb_poly_pow_series(arb_ptr h, arb_srcptr f, slong flen, arb_srcptr g, slong glen, slong len, slong prec)
void arb_poly_pow_series(arb_poly_t h, const arb_poly_t f, const arb_poly_t g, slong len, slong prec)

Set `h` to the power series $f^g = \exp(g \log f)$ truncated to length
`len`. If $g$ is constant, `_arb_poly_pow_arb_series` is used. The
underscore version requires positive lengths with $\mathit{flen},
\mathit{glen} \le \mathit{len}$ and no aliasing.
:::

::: function
void _arb_poly_pow_arb_series(arb_ptr h, arb_srcptr f, slong flen, const arb_t g, slong len, slong prec)
void arb_poly_pow_arb_series(arb_poly_t h, const arb_poly_t f, const arb_t g, slong len, slong prec)

Set `h` to $f^g = \exp(g \log f)$ for a constant exponent $g$,
truncated to length `len`. The underscore version requires positive
lengths with $\mathit{flen} \le \mathit{len}$ and no aliasing.

::: details
Special cases: if $g$ is an exact integer with $|g| < 2^{63}$, binary
exponentiation is used (followed by an inversion if $g < 0$); if $f = a
+ bx^d$ is a binomial, the coefficients come from the binomial series
(`_arb_poly_binomial_pow_arb_series`); if $g = \pm 1/2$ exactly,
`_arb_poly_sqrt_series` or `_arb_poly_rsqrt_series` is used; otherwise
$\exp(g \log f)$.
:::
:::

::: function
void _arb_poly_sqrt_series(arb_ptr g, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_sqrt_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec)

Set `g` to the power series square root of $h$, with constant term
$\sqrt{h_0}$, truncated to length `n`. Computes the reciprocal square
root to half the length and finishes with one Karp–Markstein step
[@KM1997]. The underscore version does not allow aliasing.
:::

::: function
void _arb_poly_rsqrt_series(arb_ptr g, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_rsqrt_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec)

Set `g` to the power series $h^{-1/2}$ truncated to length `n`, by the
division-free Newton iteration $g \leftarrow g + g(1 - hg^2)/2$.
Wrappers of `_gr_poly_rsqrt_series_newton`. The underscore version does
not allow aliasing.
:::

::: function
void _arb_poly_log_series(arb_ptr res, arb_srcptr f, slong flen, slong n, slong prec)
void arb_poly_log_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec)

Set `res` to $\log f = \log f_0 + \int f'/f\, dx$ truncated to length
`n`. The underscore version allows aliasing. Wrappers of
`_gr_poly_log_series`.
:::

::: function
void _arb_poly_log1p_series(arb_ptr res, arb_srcptr f, slong flen, slong n, slong prec)
void arb_poly_log1p_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec)

Set `res` to $\log(1 + f)$, which is more accurate than
`arb_poly_log_series` when the constant term of $f$ is small. Wrappers
of `_gr_poly_log1p_series`.
:::

::: function
void _arb_poly_atan_series(arb_ptr res, arb_srcptr f, slong flen, slong n, slong prec)
void arb_poly_atan_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec)
void _arb_poly_asin_series(arb_ptr res, arb_srcptr f, slong flen, slong n, slong prec)
void arb_poly_asin_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec)
void _arb_poly_acos_series(arb_ptr res, arb_srcptr f, slong flen, slong n, slong prec)
void arb_poly_acos_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec)

Set `res` to $\operatorname{atan} f$, $\operatorname{asin} f$ or
$\operatorname{acos} f$ truncated to length `n`. The underscore versions
allow aliasing. Wrappers of `_gr_poly_atan_series`,
`_gr_poly_asin_series` and `_gr_poly_acos_series`.

::: details
The generic implementations integrate the derivative,
$$\operatorname{atan} f = \operatorname{atan} f_0 + \int \frac{f'}{1 + f^2}\, dx, \qquad
\operatorname{asin} f = \operatorname{asin} f_0 + \int \frac{f'}{\sqrt{1 - f^2}}\, dx,$$
and $\operatorname{acos} f = \operatorname{acos} f_0 - \int f' /
\sqrt{1 - f^2}\, dx$.
:::
:::

::: function
void _arb_poly_atanh_series(arb_ptr res, arb_srcptr f, slong flen, slong n, slong prec)
void arb_poly_atanh_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec)
void _arb_poly_asinh_series(arb_ptr res, arb_srcptr f, slong flen, slong n, slong prec)
void arb_poly_asinh_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec)
void _arb_poly_acosh_series(arb_ptr res, arb_srcptr f, slong flen, slong n, slong prec)
void arb_poly_acosh_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec)

Set `res` to the inverse hyperbolic tangent, sine or cosine of $f$
truncated to length `n`. The underscore versions allow aliasing.
Wrappers of the corresponding `_gr_poly` functions.
:::

::: function
void _arb_poly_exp_series_basecase(arb_ptr f, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_exp_series_basecase(arb_poly_t f, const arb_poly_t h, slong n, slong prec)
void _arb_poly_exp_series(arb_ptr f, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_exp_series(arb_poly_t f, const arb_poly_t h, slong n, slong prec)

Set `f` to $\exp h$ truncated to length `n`. The underscore versions
allow aliasing.

::: details
The *basecase* version uses the recurrence $k f_k = \sum_{j=1}^{k} j
h_j f_{k-j}$ (from $f' = h' f$) with $O(n \cdot \mathit{hlen})$
operations; for long dense input at high precision it first splits off
the high part of $h$ to save operations. The default handles $h = a +
bx^d$ directly, uses the basecase version if $\mathit{hlen}$ is at most
750 (for $\mathit{prec} \le 256$) or $10^5 / (\ln \mathit{prec})^3$, and
otherwise the Newton iteration of [@HZ2004] that updates $\exp(h)$ and
$\exp(-h)$ together, at a cost of $O(M(n))$.
:::
:::

::: function
void _arb_poly_sin_cos_series(arb_ptr s, arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_sin_cos_series(arb_poly_t s, arb_poly_t c, const arb_poly_t h, slong n, slong prec)

Set `s` and `c` to $\sin h$ and $\cos h$ truncated to length `n`,
computed together. The underscore version allows aliasing. For $n \ge 3$
and nonconstant $h$ this wraps `_gr_poly_sin_cos_series`.
:::

::: function
void _arb_poly_sin_series(arb_ptr s, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_sin_series(arb_poly_t s, const arb_poly_t h, slong n, slong prec)
void _arb_poly_cos_series(arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_cos_series(arb_poly_t c, const arb_poly_t h, slong n, slong prec)

Set `s` to $\sin h$ or `c` to $\cos h$ truncated to length `n`. The
underscore versions allow aliasing. Wrappers of `_gr_poly_sin_series`
and `_gr_poly_cos_series`.
:::

::: function
void _arb_poly_tan_series(arb_ptr g, arb_srcptr h, slong hlen, slong len, slong prec)
void arb_poly_tan_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec)

Set `g` to $\tan h$ truncated to length `n`. The underscore version
does not allow aliasing. Wrappers of `_gr_poly_tan_series`, which uses
Newton iteration on $\operatorname{atan}$ for long series.
:::

::: function
void _arb_poly_cot_series(arb_ptr g, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_cot_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec)

Set `g` to $\cot h$ truncated to length `n`. The underscore version
allows aliasing. Wrappers of `_gr_poly_cot_series`.
:::

::: function
void _arb_poly_sin_cos_pi_series(arb_ptr s, arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_sin_cos_pi_series(arb_poly_t s, arb_poly_t c, const arb_poly_t h, slong n, slong prec)
void _arb_poly_sin_pi_series(arb_ptr s, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_sin_pi_series(arb_poly_t s, const arb_poly_t h, slong n, slong prec)
void _arb_poly_cos_pi_series(arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_cos_pi_series(arb_poly_t c, const arb_poly_t h, slong n, slong prec)
void _arb_poly_tan_pi_series(arb_ptr g, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_tan_pi_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec)
void _arb_poly_cot_pi_series(arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_cot_pi_series(arb_poly_t c, const arb_poly_t h, slong n, slong prec)

Compute the trigonometric functions of $\pi h$ truncated to length `n`.
These are more accurate than multiplying by $\pi$ first when the
constant term of $h$ is close to an integer or a half-integer. Wrappers
of the corresponding `_gr_poly` functions.
:::

::: function
void _arb_poly_sinh_cosh_series_basecase(arb_ptr s, arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_sinh_cosh_series_basecase(arb_poly_t s, arb_poly_t c, const arb_poly_t h, slong n, slong prec)
void _arb_poly_sinh_cosh_series_exponential(arb_ptr s, arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_sinh_cosh_series_exponential(arb_poly_t s, arb_poly_t c, const arb_poly_t h, slong n, slong prec)
void _arb_poly_sinh_cosh_series(arb_ptr s, arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_sinh_cosh_series(arb_poly_t s, arb_poly_t c, const arb_poly_t h, slong n, slong prec)
void _arb_poly_sinh_series(arb_ptr s, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_sinh_series(arb_poly_t s, const arb_poly_t h, slong n, slong prec)
void _arb_poly_cosh_series(arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_cosh_series(arb_poly_t c, const arb_poly_t h, slong n, slong prec)

Set `s` and/or `c` to $\sinh h$ and $\cosh h$ truncated to length `n`.

::: details
The *basecase* version uses the coupled recurrences that follow from $s'
= h'c$ and $c' = h's$, with $O(n \cdot \mathit{hlen})$ operations. The
*exponential* version computes $e^h$ and $e^{-h}$ (one exponential and
one inversion) and takes half their sum and difference. The default
uses the basecase version for $\mathit{hlen}$ below 400 (for
$\mathit{prec} \le 128$) or $30000 / (\ln \mathit{prec})^3$.
:::
:::

::: function
void _arb_poly_tanh_series(arb_ptr g, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_tanh_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec)
void _arb_poly_coth_series(arb_ptr g, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_coth_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec)

Set `g` to $\tanh h$ or $\coth h$ truncated to length `n`. The
underscore versions allow aliasing. Wrappers of `_gr_poly_tanh_series`
and `_gr_poly_coth_series`.
:::

::: function
void _arb_poly_sinc_series(arb_ptr s, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_sinc_series(arb_poly_t s, const arb_poly_t h, slong n, slong prec)
void _arb_poly_sinc_pi_series(arb_ptr s, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_sinc_pi_series(arb_poly_t s, const arb_poly_t h, slong n, slong prec)

Set `s` to $\operatorname{sinc} h = \sin(h)/h$ or to
$\operatorname{sinc}(\pi h)$ truncated to length `n`. A constant term
that contains zero is handled through the Taylor series of
$\operatorname{sinc}$ on a ball around zero, so that the removable
singularity does not cause a division by zero.
:::

## Lambert W, gamma and zeta functions {#special}

::: function
void _arb_poly_lambertw_series(arb_ptr res, arb_srcptr z, slong zlen, int flags, slong len, slong prec)
void arb_poly_lambertw_series(arb_poly_t res, const arb_poly_t z, int flags, slong len, slong prec)

Set `res` to $W(z)$ truncated to length `len`: the principal branch
$W_0$ if `flags` is 0, the branch $W_{-1}$ if `flags` is 1. The
underscore version allows aliasing. The value at the constant term
comes from `arb_lambertw` [@Joh2017b], and the series from Newton
iteration on $w e^w = z$.
:::

::: function
void _arb_poly_gamma_series(arb_ptr res, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_gamma_series(arb_poly_t res, const arb_poly_t h, slong n, slong prec)
void _arb_poly_rgamma_series(arb_ptr res, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_rgamma_series(arb_poly_t res, const arb_poly_t h, slong n, slong prec)
void _arb_poly_lgamma_series(arb_ptr res, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_lgamma_series(arb_poly_t res, const arb_poly_t h, slong n, slong prec)
void _arb_poly_digamma_series(arb_ptr res, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_digamma_series(arb_poly_t res, const arb_poly_t h, slong n, slong prec)

Set `res` to $\Gamma(h)$, $1/\Gamma(h)$, $\log \Gamma(h)$ or $\psi(h)$
truncated to length `n`. The underscore versions allow aliasing.

::: details
These functions compute the Taylor series at the constant term $h_0$
and compose it with $h - h_0$ by `_arb_poly_compose_series`. The Taylor
coefficients come from the zeta values $\zeta(2), \zeta(3), \ldots$ when
$h_0$ is a small integer, and otherwise from Stirling's series after
shifting the argument by a rising factorial; the Stirling sum is
evaluated as a power series by binary splitting
([](algorithms/special-functions/gamma)).
:::
:::

::: function
void _arb_poly_rising_ui_series(arb_ptr res, arb_srcptr f, slong flen, ulong r, slong trunc, slong prec)
void arb_poly_rising_ui_series(arb_poly_t res, const arb_poly_t f, ulong r, slong trunc, slong prec)

Set `res` to the rising factorial $f (f+1) \cdots (f+r-1)$ truncated to
length `trunc`, by binary splitting. The underscore version requires
`flen`, $r$ and `trunc` to be positive and no aliasing.
:::

::: function
void _arb_poly_zeta_series(arb_ptr res, arb_srcptr h, slong hlen, const arb_t a, int deflate, slong len, slong prec)
void arb_poly_zeta_series(arb_poly_t res, const arb_poly_t s, const arb_t a, int deflate, slong n, slong prec)

Set `res` to the Hurwitz zeta function $\zeta(s, a)$ of the power
series $s$, for a constant $a > 0$, truncated to length `n`; $a = 1$
gives the Riemann zeta function. If $a$ contains a nonpositive number,
the result is indeterminate. If `deflate` is nonzero, $\zeta(s, a) + 1/(1
- s)$ is computed, which is finite at $s = 1$.

::: details
The series of $\zeta(s_0 + x, a)$ is computed by Euler–Maclaurin
summation (`_acb_poly_zeta_cpx_series`) and composed with $s - s_0$
([](algorithms/special-functions/zeta)). For $a = 1$ and a constant
term with negative midpoint the functional equation is used. With $s =
1 + x$ and `deflate` set, the result is $\sum_{k=0}^{n-1}
\frac{(-1)^k}{k!} \gamma_k(a) x^k$ in terms of the Stieltjes constants
$\gamma_k(a)$.
:::
:::

::: function
void _arb_poly_riemann_siegel_theta_series(arb_ptr res, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_riemann_siegel_theta_series(arb_poly_t res, const arb_poly_t h, slong n, slong prec)

Set `res` to the Riemann–Siegel theta function
$$\theta(h) = \arg \Gamma\!\left(\frac{2ih + 1}{4}\right) - \frac{\log \pi}{2} h$$
of the power series $h$, truncated to length `n`, where the argument of
the gamma function is the imaginary part of $\log \Gamma$ (continuous
in $h$). The underscore version does not allow aliasing.
:::

::: function
void _arb_poly_riemann_siegel_z_series(arb_ptr res, arb_srcptr h, slong hlen, slong n, slong prec)
void arb_poly_riemann_siegel_z_series(arb_poly_t res, const arb_poly_t h, slong n, slong prec)

Set `res` to the Riemann–Siegel function $Z(h) = e^{i\theta(h)}
\zeta(1/2 + ih)$ of the power series $h$, truncated to length `n`. The
real zeros of $Z$ are the ordinates of the zeros of $\zeta$ on the
critical line. The underscore version allows aliasing.
:::

## Root finding {#roots}

::: function
void _arb_poly_root_bound_fujiwara(mag_t bound, arb_srcptr poly, slong len)
void arb_poly_root_bound_fujiwara(mag_t bound, arb_poly_t poly)

Set `bound` to an upper bound for the absolute values of the complex
roots of every polynomial in `poly`, using Fujiwara's bound [@Fuj1916]
$$2 \max \left\{\left|\frac{a_{n-1}}{a_n}\right|,
\left|\frac{a_{n-2}}{a_n}\right|^{1/2}, \ldots,
\left|\frac{a_1}{a_n}\right|^{1/(n-1)},
\left|\frac{a_0}{2a_n}\right|^{1/n}\right\}.$$
The bound is infinite if the leading coefficient contains zero or the
length is at most 1.
:::

::: function
void _arb_poly_newton_convergence_factor(arf_t convergence_factor, arb_srcptr poly, slong len, const arb_t convergence_interval, slong prec)

Sets `convergence_factor` to an upper bound for $C = \sup_{t, u \in I}
\tfrac12 |f''(t)| / |f'(u)|$, where $I$ is `convergence_interval` and
$f$ is `{poly, len}`, obtained by evaluating $f'$ and $f''$ on $I$. The
bound is infinite unless $f'$ is bounded away from zero on $I$; for
polynomials with large coefficients, $I$ must be narrow.
:::

::: function
int _arb_poly_newton_step(arb_t xnew, arb_srcptr poly, slong len, const arb_t x, const arb_t convergence_interval, const arf_t convergence_factor, slong prec)

Performs one step of Newton's method for a ball $x = [m \pm r]$ that
contains exactly one root of $f$ = `{poly, len}`, given a ball $I
\supseteq x$ (`convergence_interval`) and the factor $C$ computed by
`_arb_poly_newton_convergence_factor`. Computes $m' = m - f(m)/f'(m)$
in ball arithmetic and adds $C r^2$ to its radius. If the new ball $x'$
lies in $I$ and has a smaller radius than $x$, sets `xnew` to $x'$ and
returns nonzero; otherwise sets `xnew` to $x$ and returns zero.

::: details
If $\xi \in x$ is the root, Taylor's theorem gives $m' - \xi =
\frac{f''(\eta)}{2 f'(m)} (m - \xi)^2$ for some $\eta \in I$, hence
$|m' - \xi| \le C r^2$
([](algorithms/ball-arithmetic/real-root-isolation)).
:::
:::

::: function
void _arb_poly_newton_refine_root(arb_t r, arb_srcptr poly, slong len, const arb_t start, const arb_t convergence_interval, const arf_t convergence_factor, slong eval_extra_prec, slong prec)

Refines the enclosure `start` of a simple root to about `prec` bits by
Newton steps at roughly doubling precisions, as in
`arb_calc_refine_root_newton`. `eval_extra_prec` is the number of guard
bits needed to evaluate the polynomial accurately near the root (about
the bit size of the coefficients if they are large and of alternating
signs).
:::

## Other special polynomials {#special-polynomials}

::: function
void _arb_poly_swinnerton_dyer_ui(arb_ptr poly, ulong n, slong trunc, slong prec)
void arb_poly_swinnerton_dyer_ui(arb_poly_t poly, ulong n, slong prec)
slong _arb_poly_swinnerton_dyer_ui_prec(ulong n)

Set `poly` to the Swinnerton-Dyer polynomial $S_n = \prod (x \pm
\sqrt{2} \pm \sqrt{3} \pm \cdots \pm \sqrt{p_n})$ of degree $2^n$, the
minimal polynomial over $\Q$ of the sum of the square roots of the
first $n$ primes. If `prec` is 0, the precision
`_arb_poly_swinnerton_dyer_ui_prec(n)` is used, which is chosen so that
`arb_poly_get_unique_fmpz_poly` should succeed.

::: details
The underscore version computes only the first `trunc` coefficients,
but `poly` must have room for $2^n + 1$ coefficients, which are used as
temporary space. The factors are multiplied in a balanced tree of
products of polynomials in $x^2$. The precision estimate is half of a
bound derived from the binomial theorem, plus 15 bits; it is a
heuristic. `fmpz_poly_swinnerton_dyer` uses these functions.
:::
:::

## Internal functions {#internal}

::: function
void _arb_poly_binomial_pow_arb_series(arb_ptr h, arb_srcptr f, slong flen, const arb_t g, slong len, slong prec)

Sets `{h, len}` to $(a + bx^d)^g$ truncated to length `len`, where `{f,
flen}` is $a + bx^d$ with $d = \mathit{flen} - 1$ and zero coefficients
in between (not checked), using the binomial series. No aliasing.
:::

::: function
int n_zerobits(ulong e)

Returns the number of zero bits below the most significant bit of $e$,
which together with the bit length gives the number of squarings and
multiplications in binary exponentiation.
:::

::: function
slong poly_pow_length(slong poly_len, ulong exp, slong trunc)

Returns $\min((\mathit{poly\_len} - 1) e + 1, \mathit{trunc})$, the
length of the $e$-th power of a polynomial of length $\mathit{poly\_len}
\ge 1$ truncated to length `trunc`, without overflow.
:::
