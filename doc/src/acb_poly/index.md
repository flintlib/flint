---
title: acb_poly.h
header: acb_poly.h
description: Dense univariate polynomials and truncated power series with complex ball coefficients, with rigorous arithmetic, evaluation, composition, elementary and special functions of power series, and certified complex root finding.
category: Real and complex numbers
order: 7
---

This module implements dense univariate polynomials whose coefficients
are complex balls (`acb_t`, rectangles $[a \pm r] + [b \pm s]i$). An
`acb_poly_t` represents the set of complex polynomials whose
coefficients lie in the coefficient balls, and every operation returns
an enclosure of the exact result for all inputs in that set. As in
[](arb_poly), whose interface this module follows function by function,
the same objects serve as truncated power series, and the `_series`
functions compute Taylor expansions of elementary and special functions
of a power series argument.

Products of well-scaled polynomials are computed with one product of
Gaussian integer polynomials (`fmpzi`), other products by four real
products of the real and imaginary parts with `_arb_poly_mullow`
([](algorithms/ball-arithmetic/polynomial-multiplication)). Division
and the elementary functions of power series use Newton iteration, and
composition, reversion, multipoint evaluation and many transcendental
functions are wrappers around the generic implementations in
[](gr_poly) over the complex ball context. The module also contains the
Euler–Maclaurin machinery for power series of the Hurwitz zeta
function and polylogarithms, and `acb_poly_find_roots`, which computes
certified enclosures of all complex roots of a polynomial
([](acb_poly/find_roots)). The guide [](acb_poly/taylor-series) shows how
to expand special functions with this module.

::: summary algorithms/ball-arithmetic/polynomial-roots :::

## Types and conventions {#conventions}

::: type
acb_poly_struct
acb_poly_t

A polynomial with complex ball coefficients. The structure has the
fields `coeffs` (an array of `alloc` initialised `acb_struct`), `alloc`
and `length`; the coefficient of $x^i$ is `coeffs[i]`. An `acb_poly_t`
is an array of one `acb_poly_struct`, so it is passed by reference.
:::

The conventions are those of [](arb_poly#conventions):

- A polynomial is *normalised* if its length is zero or its leading
  coefficient is not exactly zero; leading coefficients that merely
  contain zero are kept, so the length and degree are upper bounds.
- `prec` is the working precision in bits for the midpoint arithmetic;
  the output is always a valid enclosure.
- Functions without a leading underscore manage memory, accept any
  lengths and allow aliasing between outputs and inputs unless stated
  otherwise. Functions with a leading underscore work on coefficient
  arrays, require positive lengths and do not allow the output to
  overlap an input unless stated otherwise.
- Functions described as wrappers of a `_gr_poly` function call it with
  `gr_ctx_init_complex_acb(ctx, prec)` (see [](gr_domains)); if the
  generic code reports failure, the output is set to indeterminate
  balls.

## Memory management {#memory}

::: function
void acb_poly_init(acb_poly_t poly)
void acb_poly_init2(acb_poly_t poly, slong len)

Initialise `poly` as the zero polynomial. `acb_poly_init2` also
allocates space for `len` coefficients.
:::

::: function
void acb_poly_clear(acb_poly_t poly)

Frees all coefficients and the coefficient array of `poly`.
:::

::: function
void acb_poly_fit_length(acb_poly_t poly, slong len)

Ensures that `poly` has at least `len` allocated and initialised
coefficients. When it has to grow, the allocation is at least doubled.
The length is not changed.
:::

::: function
void _acb_poly_set_length(acb_poly_t poly, slong len)

Sets the length of `poly` to `len`, which must not exceed the
allocation. If the length decreases, the coefficients from `len` up to
the old length are set to zero, so that unused coefficients are always
zero.
:::

::: function
void _acb_poly_normalise(acb_poly_t poly)

Decreases the length of `poly` until the leading coefficient is not the
exact zero ball. A coefficient that merely contains zero is kept.
:::

::: function
void acb_poly_swap(acb_poly_t poly1, acb_poly_t poly2)

Swaps `poly1` and `poly2` efficiently.
:::

::: function
slong acb_poly_allocated_bytes(const acb_poly_t x)

Returns the number of bytes allocated on the heap by `x` (the
coefficient array and the limbs of the coefficients), excluding
`sizeof(acb_poly_struct)`.
:::

## Basic manipulation {#basic}

::: function
slong acb_poly_length(const acb_poly_t poly)
slong acb_poly_degree(const acb_poly_t poly)

Return the length of `poly` and its degree, the length minus one ($-1$
for the zero polynomial).

::: details
If leading coefficients are balls that contain zero, the exact
polynomials represented by `poly` may have smaller degree; the returned
degree is then only an upper bound.
:::
:::

::: function
int acb_poly_is_zero(const acb_poly_t poly)
int acb_poly_is_one(const acb_poly_t poly)
int acb_poly_is_x(const acb_poly_t poly)

Return whether `poly` is exactly the polynomial 0, 1 or $x$ (with exact
coefficients).
:::

::: function
void acb_poly_zero(acb_poly_t poly)
void acb_poly_one(acb_poly_t poly)

Set `poly` to the constant 0 or 1.
:::

::: function
void acb_poly_set(acb_poly_t dest, const acb_poly_t src)
void acb_poly_set_round(acb_poly_t dest, const acb_poly_t src, slong prec)

Set `dest` to a copy of `src`, either exactly or with the midpoints
rounded to `prec` bits (the rounding error is added to the radii).
:::

::: function
void acb_poly_set_trunc(acb_poly_t dest, const acb_poly_t src, slong n)
void acb_poly_set_trunc_round(acb_poly_t dest, const acb_poly_t src, slong n, slong prec)

Set `dest` to `src` truncated to length `n` (reduced modulo $x^n$),
exactly or rounded to `prec` bits.
:::

::: function
void acb_poly_set_coeff_si(acb_poly_t poly, slong n, slong c)
void acb_poly_set_coeff_acb(acb_poly_t poly, slong n, const acb_t c)

Set the coefficient of $x^n$ in `poly` to `c`, extending the
polynomial with zeros if necessary, and normalise. Require $n \ge 0$.
:::

::: function
void acb_poly_get_coeff_acb(acb_t v, const acb_poly_t poly, slong n)

Sets `v` to the coefficient of $x^n$ in `poly`, which is zero if $n$ is
at least the length. Requires $n \ge 0$.
:::

::: macro
acb_poly_get_coeff_ptr(poly, n)

Evaluates to a pointer to the coefficient of $x^n$ in `poly`, or to
`NULL` if $n$ is at least the length. Requires $n \ge 0$.
:::

::: function
void _acb_poly_reverse(acb_ptr res, acb_srcptr poly, slong len, slong n)

Sets `{res, n}` to the reversal of `{poly, len}` as a polynomial of
length $n \ge \mathit{len}$, that is, `res[i] = poly[n - 1 - i]` with
zeros where the index is out of range. Aliasing is allowed. A wrapper
of `_gr_poly_reverse`.
:::

::: function
void _acb_poly_shift_right(acb_ptr res, acb_srcptr poly, slong len, slong n)
void acb_poly_shift_right(acb_poly_t res, const acb_poly_t poly, slong n)

Set `res` to `poly` divided by $x^n$, discarding the coefficients of
$x^0, \ldots, x^{n-1}$. Require $n \ge 0$. The underscore version
writes $\mathit{len} - n$ coefficients, requires $n \le \mathit{len}$
and allows `res` to be `poly`. Wrappers of `_gr_poly_shift_right`.
:::

::: function
void _acb_poly_shift_left(acb_ptr res, acb_srcptr poly, slong len, slong n)
void acb_poly_shift_left(acb_poly_t res, const acb_poly_t poly, slong n)

Set `res` to `poly` multiplied by $x^n$. Require $n \ge 0$. The
underscore version writes $\mathit{len} + n$ coefficients and allows
`res` to be `poly`. Wrappers of `_gr_poly_shift_left`.
:::

::: function
void acb_poly_truncate(acb_poly_t poly, slong n)

Truncates `poly` to length at most `n` (reduces it modulo $x^n$) and
normalises it. Requires $n \ge 0$.
:::

::: function
slong acb_poly_valuation(const acb_poly_t poly)

Returns the index of the first coefficient of `poly` that is not the
exact zero ball, or $-1$ if `poly` is zero. A coefficient that merely
contains zero counts as nonzero, so the result is a lower bound for the
valuation of the exact polynomials represented.
:::

## Conversions {#conversions}

::: function
void acb_poly_set_fmpz_poly(acb_poly_t poly, const fmpz_poly_t re, slong prec)
void acb_poly_set_fmpq_poly(acb_poly_t poly, const fmpq_poly_t re, slong prec)
void acb_poly_set_arb_poly(acb_poly_t poly, const arb_poly_t re)

Set `poly` to the real polynomial `re`, rounding integer or rational
coefficients to `prec` bits; `acb_poly_set_arb_poly` copies the balls
exactly.
:::

::: function
void acb_poly_set2_fmpz_poly(acb_poly_t poly, const fmpz_poly_t re, const fmpz_poly_t im, slong prec)
void acb_poly_set2_fmpq_poly(acb_poly_t poly, const fmpq_poly_t re, const fmpq_poly_t im, slong prec)
void acb_poly_set2_arb_poly(acb_poly_t poly, const arb_poly_t re, const arb_poly_t im)

Set `poly` to $\mathit{re} + i \cdot \mathit{im}$, rounding integer or
rational coefficients to `prec` bits.
:::

::: function
void acb_poly_set_acb(acb_poly_t poly, const acb_t src)
void acb_poly_set_si(acb_poly_t poly, slong src)

Set `poly` to the constant polynomial `src`, exactly. The result has
length zero if `src` is exactly zero.
:::

## Input and output {#io}

::: function
void acb_poly_printd(const acb_poly_t poly, slong digits)
void acb_poly_fprintd(FILE * file, const acb_poly_t poly, slong digits)

Print `poly` to standard output or to `file` as a list `[(c0), (c1),
...]` of its coefficients from the constant term upwards, each printed
with `acb_fprintd` to `digits` decimal digits.
:::

## Random generation {#random}

::: function
void acb_poly_randtest(acb_poly_t poly, flint_rand_t state, slong len, slong prec, slong mag_bits)

Sets `poly` to a random polynomial of length at most `len` whose
coefficients are generated by `acb_randtest` or `acb_randtest_precise`
with the parameters `prec` and `mag_bits`. For testing.
:::

## Comparisons {#comparisons}

::: function
int acb_poly_equal(const acb_poly_t A, const acb_poly_t B)

Returns nonzero if `A` and `B` have the same length and identical
coefficient balls (equal midpoints and radii of the real and imaginary
parts).
:::

::: function
int acb_poly_contains(const acb_poly_t poly1, const acb_poly_t poly2)
int acb_poly_contains_fmpz_poly(const acb_poly_t poly1, const fmpz_poly_t poly2)
int acb_poly_contains_fmpq_poly(const acb_poly_t poly1, const fmpq_poly_t poly2)

Return nonzero if every coefficient of `poly2` is contained in the
corresponding coefficient of `poly1`, coefficients beyond the length of
`poly2` being zero.
:::

::: function
int _acb_poly_overlaps(acb_srcptr poly1, slong len1, acb_srcptr poly2, slong len2)
int acb_poly_overlaps(const acb_poly_t poly1, const acb_poly_t poly2)

Return nonzero if the coefficients of `poly1` and `poly2` overlap
pairwise (coefficients of the longer polynomial beyond the length of the
shorter one must contain zero). The underscore version requires
$\mathit{len1} \ge \mathit{len2}$.
:::

::: function
int acb_poly_get_unique_fmpz_poly(fmpz_poly_t z, const acb_poly_t x)

If, for every coefficient of `x`, the real part contains exactly one
integer and the imaginary part contains zero, sets `z` to the integer
polynomial with these coefficients and returns nonzero; `z` is then the
only integer polynomial contained in `x`. Otherwise returns zero, and
`z` may have been modified.
:::

::: function
int acb_poly_is_real(const acb_poly_t poly)

Returns nonzero if all coefficients of `poly` have imaginary part
exactly zero.
:::

## Bounds {#bounds}

::: function
void _acb_poly_majorant(arb_ptr res, acb_srcptr poly, slong len, slong prec)
void acb_poly_majorant(arb_poly_t res, const acb_poly_t poly, slong prec)

Set `res` to the real polynomial with exact nonnegative coefficients
that are upper bounds, rounded to `prec` bits, for the absolute values
of the coefficients of `poly`. For $|x| \le t$, its value at $t$ bounds
the absolute value at $x$ of every polynomial in `poly`.
:::

## Addition and scalar arithmetic {#arithmetic}

::: function
void _acb_poly_add(acb_ptr C, acb_srcptr A, slong lenA, acb_srcptr B, slong lenB, slong prec)
void _acb_poly_sub(acb_ptr C, acb_srcptr A, slong lenA, acb_srcptr B, slong lenB, slong prec)

Set `{C, max(lenA, lenB)}` to the sum or difference of `{A, lenA}` and
`{B, lenB}`. The output may be the same array as either input, and the
lengths may be zero.
:::

::: function
void acb_poly_add(acb_poly_t C, const acb_poly_t A, const acb_poly_t B, slong prec)
void acb_poly_sub(acb_poly_t C, const acb_poly_t A, const acb_poly_t B, slong prec)
void acb_poly_add_si(acb_poly_t C, const acb_poly_t A, slong B, slong prec)

Set `C` to $A + B$ or $A - B$.
:::

::: function
void acb_poly_add_series(acb_poly_t C, const acb_poly_t A, const acb_poly_t B, slong len, slong prec)
void acb_poly_sub_series(acb_poly_t C, const acb_poly_t A, const acb_poly_t B, slong len, slong prec)

Set `C` to $A + B$ or $A - B$ truncated to length `len`, without
computing the discarded coefficients.
:::

::: function
void acb_poly_neg(acb_poly_t C, const acb_poly_t A)

Sets `C` to $-A$, exactly.
:::

::: function
void acb_poly_scalar_mul_2exp_si(acb_poly_t C, const acb_poly_t A, slong c)

Sets `C` to $2^c A$, exactly.
:::

::: function
void acb_poly_scalar_mul(acb_poly_t C, const acb_poly_t A, const acb_t c, slong prec)
void acb_poly_scalar_mul_si(acb_poly_t C, const acb_poly_t A, slong c, slong prec)
void acb_poly_scalar_div(acb_poly_t C, const acb_poly_t A, const acb_t c, slong prec)

Set `C` to $c A$ or $A / c$. Division by a ball that contains zero gives
non-finite coefficients.
:::

## Multiplication {#multiplication}

::: function
void _acb_poly_mullow_classical(acb_ptr C, acb_srcptr A, slong lenA, acb_srcptr B, slong lenB, slong n, slong prec)
void _acb_poly_mullow_transpose(acb_ptr C, acb_srcptr A, slong lenA, acb_srcptr B, slong lenB, slong n, slong prec)
void _acb_poly_mullow_transpose_gauss(acb_ptr C, acb_srcptr A, slong lenA, acb_srcptr B, slong lenB, slong n, slong prec)
void _acb_poly_mullow(acb_ptr C, acb_srcptr A, slong lenA, acb_srcptr B, slong lenB, slong n, slong prec)

Set `{C, n}` to the product of `{A, lenA}` and `{B, lenB}` truncated to
length `n`. Require $\mathit{lenA}, \mathit{lenB} \ge 1$ and $1 \le n
\le \mathit{lenA} + \mathit{lenB} - 1$; the output must not overlap
either input. If `A` and `B` are the same pointer with the same length,
the square is computed.

::: details
The *classical* version computes each coefficient as one `acb_dot`.
The *transpose* version writes $A = a + bi$ and $B = c + di$ with real
ball polynomials and computes $(ac - bd) + (ad + bc)i$ with four calls
to `_arb_poly_mullow` (three for a square). The *transpose_gauss*
version uses Gauss's three-multiplication formula $ad + bc = (a+b)(c+d)
- ac - bd$; it is faster but gives much larger radii when the real and
imaginary parts differ greatly in magnitude, and it is not used by
default.

The default `_acb_poly_mullow` first asks `_acb_poly_mulmid_use_fmpzi`
whether a single Gaussian integer product is worthwhile and, if so,
tries `_acb_poly_mulmid_fmpzi`, which succeeds for well-scaled input.
Otherwise it uses the classical version for short operands and the
transposed version for long ones, with the cutoffs of
`_arb_poly_mullow`.
:::
:::

::: function
void acb_poly_mullow_classical(acb_poly_t C, const acb_poly_t A, const acb_poly_t B, slong n, slong prec)
void acb_poly_mullow_transpose(acb_poly_t C, const acb_poly_t A, const acb_poly_t B, slong n, slong prec)
void acb_poly_mullow_transpose_gauss(acb_poly_t C, const acb_poly_t A, const acb_poly_t B, slong n, slong prec)
void acb_poly_mullow(acb_poly_t C, const acb_poly_t A, const acb_poly_t B, slong n, slong prec)

Set `C` to $AB$ truncated to length `n`. If `A` and `B` are the same
object, the square is computed.
:::

::: function
void _acb_poly_mul(acb_ptr C, acb_srcptr A, slong lenA, acb_srcptr B, slong lenB, slong prec)
void acb_poly_mul(acb_poly_t C, const acb_poly_t A1, const acb_poly_t B2, slong prec)

Set `C` to the full product $AB$. The underscore version writes
$\mathit{lenA} + \mathit{lenB} - 1$ coefficients, requires positive
lengths and no aliasing, and calls `_acb_poly_mullow`.
:::

::: function
void _acb_poly_mulmid_classical(acb_ptr z, acb_srcptr x, slong xlen, acb_srcptr y, slong ylen, slong nlo, slong nhi, slong prec)
void _acb_poly_mulmid_transpose(acb_ptr z, acb_srcptr x, slong xlen, acb_srcptr y, slong ylen, slong nlo, slong nhi, slong prec)
void _acb_poly_mulmid(acb_ptr z, acb_srcptr x, slong xlen, acb_srcptr y, slong ylen, slong nlo, slong nhi, slong prec)
void acb_poly_mulmid_classical(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, slong nlo, slong nhi, slong prec)
void acb_poly_mulmid_transpose(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, slong nlo, slong nhi, slong prec)
void acb_poly_mulmid(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, slong nlo, slong nhi, slong prec)

Set `z` (or `res`) to the coefficients of $x^{\mathit{nlo}}, \ldots,
x^{\mathit{nhi}-1}$ of the product, shifted down to start at $x^0$. The
underscore versions write $\mathit{nhi} - \mathit{nlo}$ coefficients
and require positive lengths, $0 \le \mathit{nlo} < \mathit{nhi} \le
\mathit{xlen} + \mathit{ylen} - 1$ and no aliasing. The algorithms and
the choice between them are those of the `mullow` functions.
:::

::: function
int _acb_poly_mulmid_fmpzi(acb_ptr z, acb_srcptr x, slong xlen, acb_srcptr y, slong ylen, slong nlo, slong nhi, slong prec)

If the inputs are well-scaled, sets `{z, nhi - nlo}` to the
coefficients of $x^{\mathit{nlo}}, \ldots, x^{\mathit{nhi}-1}$ of the
product of `{x, xlen}` and `{y, ylen}` (the square if the pointers and
lengths are equal) and returns 1; otherwise returns 0 and leaves `z`
undefined. The window may extend beyond the product, and `z` may alias
an input.

::: details
After the substitution $x \to 2^{-c} x$, with $c$ chosen from the first
and last nonzero midpoints as in `_arb_poly_mulmid_block` (using the
larger of the real and imaginary parts), the real and imaginary
midpoints of each operand are written exactly as Gaussian integers
times one common power of two. The input is well-scaled if this takes at
most $3\,\mathit{prec} + 512$ bits for each operand, the height that
`_arb_poly_mulmid_block` allows for one block. The midpoints of the
result are then computed by one `_fmpzi_poly_mulmid` and rounded to
`prec` bits. The radii are bounded separately for the real and
imaginary parts: for $(a_r + a_i i)(b_r + b_i i)$ with radii $\rho_{a_r},
\ldots$, the real part gets
$|a_r| \rho_{b_r} + \rho_{a_r}(|b_r| + \rho_{b_r}) + |a_i| \rho_{b_i} + \rho_{a_i}(|b_i| + \rho_{b_i})$
and the imaginary part
$|a_r| \rho_{b_i} + \rho_{a_r}(|b_i| + \rho_{b_i}) + |a_i| \rho_{b_r} + \rho_{a_i}(|b_r| + \rho_{b_r})$,
summed over the convolution with up to eight products of nonnegative
polynomials (four for a square) by `_arb_poly_addmulmid_mag`.
:::
:::

::: function
int _acb_poly_mulmid_use_fmpzi(slong len1, slong len2, slong nlo, slong nhi, slong prec, int squaring)

Returns whether `_acb_poly_mulmid` and `_acb_poly_mullow` should try
`_acb_poly_mulmid_fmpzi` first: this is the case when the number $q$ of
coefficient products in the window per input coefficient satisfies $q
\ge 1$ and $q \cdot \mathit{prec}^{0.57} \ge 240$ (330 for squares).
Both lengths must be at least 2.
:::

::: function
void _acb_poly_mul_monic(acb_ptr res, acb_srcptr poly1, slong len1, acb_srcptr poly2, slong len2, slong prec)

Sets `{res, len1 + len2 - 1}` to the product of the monic polynomials
`{poly1, len1}` and `{poly2, len2}`, with leading coefficient exactly
1. No aliasing.
:::

## Division {#division}

::: function
void _acb_poly_inv_series(acb_ptr Q, acb_srcptr A, slong Alen, slong len, slong prec)
void acb_poly_inv_series(acb_poly_t Q, const acb_poly_t A, slong n, slong prec)

Set `Q` to the power series inverse $1/A$ truncated to length `n`
(`len`). The result is finite only if the constant term of $A$ excludes
zero; it is indeterminate if $A$ is zero. The underscore version
requires positive lengths and no aliasing.

::: details
The first coefficients are computed by the recurrence $q_i = -q_0
\sum_{j=1}^{i} a_j q_{i-j}$ with `acb_dot` (all of them if $n \le 8$ or
$A$ has length 2, otherwise the first 4), and the rest by the Newton
iteration $Q \leftarrow Q - Q(AQ - 1)$, which doubles the length at each
step and computes only the unknown part of $AQ - 1$ by a middle product.
The recurrence is faster up to much larger lengths but gives worse
error bounds ([](algorithms/ball-arithmetic/power-series)).
:::
:::

::: function
void _acb_poly_div_series(acb_ptr Q, acb_srcptr A, slong Alen, acb_srcptr B, slong Blen, slong n, slong prec)
void acb_poly_div_series(acb_poly_t Q, const acb_poly_t A, const acb_poly_t B, slong n, slong prec)

Set `Q` to the power series quotient $A / B$ truncated to length `n`.
For $n \le 10$ or $B$ of length 2 the coefficients are computed by the
recurrence, otherwise as $A \cdot (1/B)$ with `_acb_poly_inv_series`.
The underscore version requires positive lengths and no aliasing.
:::

::: function
void _acb_poly_div(acb_ptr Q, acb_srcptr A, slong lenA, acb_srcptr B, slong lenB, slong prec)
void _acb_poly_rem(acb_ptr R, acb_srcptr A, slong lenA, acb_srcptr B, slong lenB, slong prec)
void _acb_poly_divrem(acb_ptr Q, acb_ptr R, acb_srcptr A, slong lenA, acb_srcptr B, slong lenB, slong prec)

Compute the quotient `{Q, lenA - lenB + 1}` and/or the remainder `{R,
lenB - 1}` of the division of `{A, lenA}` by `{B, lenB}`. Require
$\mathit{lenA} \ge \mathit{lenB} \ge 1$ and no aliasing. If the leading
coefficient of $B$ contains zero, the output is not finite.
:::

::: function
int acb_poly_divrem(acb_poly_t Q, acb_poly_t R, const acb_poly_t A, const acb_poly_t B, slong prec)

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
void _acb_poly_div_root(acb_ptr Q, acb_t R, acb_srcptr A, slong len, const acb_t c, slong prec)

Divides `{A, len}` by $x - c$ by Horner's scheme, writing the quotient
`{Q, len - 1}` and the remainder $R = A(c)$. Requires $\mathit{len} \ge
1$.
:::

## Composition {#composition}

::: function
void _acb_poly_taylor_shift(acb_ptr g, const acb_t c, slong n, slong prec)
void acb_poly_taylor_shift(acb_poly_t g, const acb_poly_t f, const acb_t c, slong prec)

Set `g` to the Taylor shift $f(x + c)$. The underscore version works in
place on `{g, n}`. Wrappers of `_gr_poly_taylor_shift`.

::: details
Over the complex ball context the generic function dispatches to
`_gr_acb_poly_taylor_shift`, which uses Horner's scheme
(`_gr_poly_taylor_shift_horner`) if $n \le 30$, if $n \le 100$ and $c$
has fewer than $\mathit{prec}/100$ bits, or if $n \le 500$, $c$ has a
one-bit midpoint and $n < 30 + 3\sqrt{\mathit{prec}}$; otherwise the
convolution algorithm (`_gr_poly_taylor_shift_convolution`, one
polynomial product) if $\mathit{prec} > 2n$, and the divide-and-conquer
algorithm (`_gr_poly_taylor_shift_divconquer`) if not.
:::
:::

::: function
void _acb_poly_compose(acb_ptr res, acb_srcptr poly1, slong len1, acb_srcptr poly2, slong len2, slong prec)
void acb_poly_compose(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, slong prec)

Set `res` to the composition $f(g(x))$ of $f$ = `poly1` and $g$ =
`poly2`. The underscore version writes $(\mathit{len1} -
1)(\mathit{len2} - 1) + 1$ coefficients and does not allow aliasing.
Wrappers of `_gr_poly_compose_divconquer`
([](algorithms/polynomials/composition)).
:::

::: function
void _acb_poly_compose_series(acb_ptr res, acb_srcptr poly1, slong len1, acb_srcptr poly2, slong len2, slong n, slong prec)
void acb_poly_compose_series(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, slong n, slong prec)

Set `res` to the power series composition $f(g(x)) \bmod x^n$. The
constant term of $g$ must be exactly zero; `acb_poly_compose_series`
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
void _acb_poly_revert_series(acb_ptr h, acb_srcptr f, slong flen, slong n, slong prec)
void acb_poly_revert_series(acb_poly_t h, const acb_poly_t f, slong n, slong prec)

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
void _acb_poly_evaluate_horner(acb_t y, acb_srcptr f, slong len, const acb_t x, slong prec)
void acb_poly_evaluate_horner(acb_t y, const acb_poly_t f, const acb_t x, slong prec)
void _acb_poly_evaluate_rectangular(acb_t y, acb_srcptr f, slong len, const acb_t x, slong prec)
void acb_poly_evaluate_rectangular(acb_t y, const acb_poly_t f, const acb_t x, slong prec)
void _acb_poly_evaluate(acb_t y, acb_srcptr f, slong len, const acb_t x, slong prec)
void acb_poly_evaluate(acb_t y, const acb_poly_t f, const acb_t x, slong prec)

Set $y = f(x)$ using Horner's rule, rectangular splitting [@PS1973], or
an automatic choice: rectangular splitting when $\mathit{prec} \ge
1024$, $\mathit{len} \ge 5 + 20000/\mathit{prec}$ and the coefficients
have at most $\mathit{prec}/2$ bits, Horner's rule otherwise. The
underscore versions accept $\mathit{len} = 0$.
:::

::: function
void _acb_poly_evaluate2_horner(acb_t y, acb_t z, acb_srcptr f, slong len, const acb_t x, slong prec)
void acb_poly_evaluate2_horner(acb_t y, acb_t z, const acb_poly_t f, const acb_t x, slong prec)
void _acb_poly_evaluate2_rectangular(acb_t y, acb_t z, acb_srcptr f, slong len, const acb_t x, slong prec)
void acb_poly_evaluate2_rectangular(acb_t y, acb_t z, const acb_poly_t f, const acb_t x, slong prec)
void _acb_poly_evaluate2(acb_t y, acb_t z, acb_srcptr f, slong len, const acb_t x, slong prec)
void acb_poly_evaluate2(acb_t y, acb_t z, const acb_poly_t f, const acb_t x, slong prec)

Set $y = f(x)$ and $z = f'(x)$ together, using Horner's rule,
rectangular splitting (which reuses the powers of $x$ for both values)
or the automatic choice of `_acb_poly_evaluate`.
:::

## Product trees {#trees}

::: function
void _acb_poly_product_roots(acb_ptr poly, acb_srcptr xs, slong n, slong prec)
void acb_poly_product_roots(acb_poly_t poly, acb_srcptr xs, slong n, slong prec)

Set `poly` to $(x - x_0)(x - x_1) \cdots (x - x_{n-1})$. The underscore
version writes $n + 1$ coefficients. Wrappers of
`_gr_poly_product_roots`, which multiplies the factors in a balanced
tree. For a real polynomial with roots in conjugate pairs,
`arb_poly_product_roots_complex` avoids the spurious imaginary parts.
:::

::: function
acb_ptr * _acb_poly_tree_alloc(slong len)
void _acb_poly_tree_free(acb_ptr * tree, slong len)
void _acb_poly_tree_build(acb_ptr * tree, acb_srcptr roots, slong len, slong prec)

Allocate, free and fill a subproduct tree for the factors $x - x_i$ of
`len` points ([](algorithms/polynomials/subproduct-tree)). Wrappers of
`_gr_poly_tree_alloc`, `_gr_poly_tree_free` and `_gr_poly_tree_build`.
:::

## Multipoint evaluation and interpolation {#multipoint}

::: function
void _acb_poly_evaluate_vec_iter(acb_ptr ys, acb_srcptr poly, slong plen, acb_srcptr xs, slong n, slong prec)
void acb_poly_evaluate_vec_iter(acb_ptr ys, const acb_poly_t poly, acb_srcptr xs, slong n, slong prec)

Set `ys[i]` to the value of `poly` at `xs[i]` for $0 \le i < n$ by
calling `_acb_poly_evaluate` for each point.
:::

::: function
void _acb_poly_evaluate_vec_fast_precomp(acb_ptr vs, acb_srcptr poly, slong plen, acb_ptr * tree, slong len, slong prec)
void _acb_poly_evaluate_vec_fast(acb_ptr ys, acb_srcptr poly, slong plen, acb_srcptr xs, slong n, slong prec)
void acb_poly_evaluate_vec_fast(acb_ptr ys, const acb_poly_t poly, acb_srcptr xs, slong n, slong prec)

Evaluate `poly` at $n$ points by the remainder tree algorithm, given
the points or a precomputed subproduct tree. Wrappers of
`_gr_poly_evaluate_vec_fast_precomp`. The remainder tree is
numerically unstable: unless the precision grows with the number of
points, the radii become large.
:::

::: function
void _acb_poly_interpolate_newton(acb_ptr poly, acb_srcptr xs, acb_srcptr ys, slong n, slong prec)
void acb_poly_interpolate_newton(acb_poly_t poly, acb_srcptr xs, acb_srcptr ys, slong n, slong prec)

Set `poly` to the polynomial of length at most $n$ with
$\mathit{poly}(x_i) = y_i$, by computing divided differences (the
Newton form) and converting to the monomial basis. If two of the points
overlap as balls, the result is not finite.
:::

::: function
void _acb_poly_interpolate_barycentric(acb_ptr poly, acb_srcptr xs, acb_srcptr ys, slong n, slong prec)
void acb_poly_interpolate_barycentric(acb_poly_t poly, acb_srcptr xs, acb_srcptr ys, slong n, slong prec)

Interpolate as `acb_poly_interpolate_newton`, using the barycentric
form of the Lagrange interpolation formula.
:::

::: function
void _acb_poly_interpolation_weights(acb_ptr w, acb_ptr * tree, slong len, slong prec)
void _acb_poly_interpolate_fast_precomp(acb_ptr poly, acb_srcptr ys, acb_ptr * tree, acb_srcptr weights, slong len, slong prec)
void _acb_poly_interpolate_fast(acb_ptr poly, acb_srcptr xs, acb_srcptr ys, slong len, slong prec)
void acb_poly_interpolate_fast(acb_poly_t poly, acb_srcptr xs, acb_srcptr ys, slong n, slong prec)

Interpolate by the fast Lagrange formula on a subproduct tree.
`_acb_poly_interpolation_weights` computes the weights $w_i = 1 /
\prod_{j \ne i} (x_i - x_j)$ from a tree, and the *precomp* version
takes the tree and the weights as input. Wrappers of
`_gr_poly_interpolation_weights` and
`_gr_poly_interpolate_fast_precomp`.
:::

## Differentiation and integration {#calculus}

::: function
void _acb_poly_derivative(acb_ptr res, acb_srcptr poly, slong len, slong prec)
void acb_poly_derivative(acb_poly_t res, const acb_poly_t poly, slong prec)

Set `res` to the derivative of `poly`. The underscore version writes
$\mathit{len} - 1$ coefficients and allows aliasing. Wrappers of
`_gr_poly_derivative`.
:::

::: function
void _acb_poly_nth_derivative(acb_ptr res, acb_srcptr poly, ulong n, slong len, slong prec)
void acb_poly_nth_derivative(acb_poly_t res, const acb_poly_t poly, ulong n, slong prec)

Set `res` to the $n$-th derivative of `poly`. The underscore version
writes $\mathit{len} - n$ coefficients (nothing if $\mathit{len} \le
n$) and allows aliasing. Wrappers of `_gr_poly_nth_derivative`.
:::

::: function
void _acb_poly_integral(acb_ptr res, acb_srcptr poly, slong len, slong prec)
void acb_poly_integral(acb_poly_t res, const acb_poly_t poly, slong prec)

Set `res` to the integral of `poly` with zero constant term. The
underscore version reads `{poly, len - 1}`, writes `{res, len}` and
allows aliasing. Wrappers of `_gr_poly_integral`.
:::

## Transforms {#transforms}

::: function
void _acb_poly_borel_transform(acb_ptr res, acb_srcptr poly, slong len, slong prec)
void acb_poly_borel_transform(acb_poly_t res, const acb_poly_t poly, slong prec)
void _acb_poly_inv_borel_transform(acb_ptr res, acb_srcptr poly, slong len, slong prec)
void acb_poly_inv_borel_transform(acb_poly_t res, const acb_poly_t poly, slong prec)

Map $\sum_k a_k x^k$ to $\sum_k (a_k / k!) x^k$ (Borel transform) or to
$\sum_k a_k k!\, x^k$ (inverse Borel transform). The underscore
versions allow aliasing.
:::

::: function
void _acb_poly_binomial_transform_basecase(acb_ptr b, acb_srcptr a, slong alen, slong len, slong prec)
void acb_poly_binomial_transform_basecase(acb_poly_t b, const acb_poly_t a, slong len, slong prec)
void _acb_poly_binomial_transform_convolution(acb_ptr b, acb_srcptr a, slong alen, slong len, slong prec)
void acb_poly_binomial_transform_convolution(acb_poly_t b, const acb_poly_t a, slong len, slong prec)
void _acb_poly_binomial_transform(acb_ptr b, acb_srcptr a, slong alen, slong len, slong prec)
void acb_poly_binomial_transform(acb_poly_t b, const acb_poly_t a, slong len, slong prec)

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
void _acb_poly_graeffe_transform(acb_ptr b, acb_srcptr a, slong len, slong prec)
void acb_poly_graeffe_transform(acb_poly_t b, const acb_poly_t a, slong prec)

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
example, `acb_poly_log_series` needs a constant term that excludes zero.
See [](algorithms/ball-arithmetic/power-series) for the algorithms.

::: function
void _acb_poly_pow_ui_trunc_binexp(acb_ptr res, acb_srcptr f, slong flen, ulong exp, slong len, slong prec)
void acb_poly_pow_ui_trunc_binexp(acb_poly_t res, const acb_poly_t poly, ulong exp, slong len, slong prec)

Set `res` to $f^{e}$ truncated to length `len`, by binary
exponentiation with truncated products. The underscore version requires
`len` to be at most the length of the untruncated power, positive
lengths and no aliasing. Wrappers of `_gr_poly_pow_series_ui_binexp`.
:::

::: function
void _acb_poly_pow_ui(acb_ptr res, acb_srcptr f, slong flen, ulong exp, slong prec)
void acb_poly_pow_ui(acb_poly_t res, const acb_poly_t poly, ulong exp, slong prec)

Set `res` to $f^{e}$. The underscore version writes $e(\mathit{flen} -
1) + 1$ coefficients, requires $\mathit{flen} \ge 1$ and no aliasing.
:::

::: function
void _acb_poly_pow_series(acb_ptr h, acb_srcptr f, slong flen, acb_srcptr g, slong glen, slong len, slong prec)
void acb_poly_pow_series(acb_poly_t h, const acb_poly_t f, const acb_poly_t g, slong len, slong prec)

Set `h` to the power series $f^g = \exp(g \log f)$ truncated to length
`len`. If $g$ is constant, `_acb_poly_pow_acb_series` is used. The
underscore version requires positive lengths with $\mathit{flen},
\mathit{glen} \le \mathit{len}$ and no aliasing.
:::

::: function
void _acb_poly_pow_acb_series(acb_ptr h, acb_srcptr f, slong flen, const acb_t g, slong len, slong prec)
void acb_poly_pow_acb_series(acb_poly_t h, const acb_poly_t f, const acb_t g, slong len, slong prec)

Set `h` to $f^g = \exp(g \log f)$ for a constant exponent $g$,
truncated to length `len`. The underscore version requires positive
lengths with $\mathit{flen} \le \mathit{len}$ and no aliasing.

::: details
Special cases: if $g$ is an exact real integer with $|g| < 2^{63}$, binary
exponentiation is used (followed by an inversion if $g < 0$); if $f = a
+ bx^d$ is a binomial, the coefficients come from the binomial series
(`_acb_poly_binomial_pow_acb_series`); if $g = \pm 1/2$ exactly,
`_acb_poly_sqrt_series` or `_acb_poly_rsqrt_series` is used; otherwise
$\exp(g \log f)$.
:::
:::

::: function
void _acb_poly_sqrt_series(acb_ptr g, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_sqrt_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec)

Set `g` to the power series square root of $h$, with constant term
$\sqrt{h_0}$, truncated to length `n`. Computes the reciprocal square
root to half the length and finishes with one Karp–Markstein step
[@KM1997]. The underscore version does not allow aliasing.
:::

::: function
void _acb_poly_rsqrt_series(acb_ptr g, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_rsqrt_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec)

Set `g` to the power series $h^{-1/2}$ truncated to length `n`, by the
division-free Newton iteration $g \leftarrow g + g(1 - hg^2)/2$.
Wrappers of `_gr_poly_rsqrt_series_newton`. The underscore version does
not allow aliasing.
:::

::: function
void _acb_poly_log_series(acb_ptr res, acb_srcptr f, slong flen, slong n, slong prec)
void acb_poly_log_series(acb_poly_t res, const acb_poly_t f, slong n, slong prec)

Set `res` to $\log f = \log f_0 + \int f'/f\, dx$ truncated to length
`n`. The underscore version allows aliasing. Wrappers of
`_gr_poly_log_series`.
:::

::: function
void _acb_poly_log1p_series(acb_ptr res, acb_srcptr f, slong flen, slong n, slong prec)
void acb_poly_log1p_series(acb_poly_t res, const acb_poly_t f, slong n, slong prec)

Set `res` to $\log(1 + f)$, which is more accurate than
`acb_poly_log_series` when the constant term of $f$ is small. Wrappers
of `_gr_poly_log1p_series`.
:::

::: function
void _acb_poly_atan_series(acb_ptr res, acb_srcptr f, slong flen, slong n, slong prec)
void acb_poly_atan_series(acb_poly_t res, const acb_poly_t f, slong n, slong prec)
void _acb_poly_asin_series(acb_ptr res, acb_srcptr f, slong flen, slong n, slong prec)
void acb_poly_asin_series(acb_poly_t res, const acb_poly_t f, slong n, slong prec)
void _acb_poly_acos_series(acb_ptr res, acb_srcptr f, slong flen, slong n, slong prec)
void acb_poly_acos_series(acb_poly_t res, const acb_poly_t f, slong n, slong prec)

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
void _acb_poly_atanh_series(acb_ptr res, acb_srcptr f, slong flen, slong n, slong prec)
void acb_poly_atanh_series(acb_poly_t res, const acb_poly_t f, slong n, slong prec)
void _acb_poly_asinh_series(acb_ptr res, acb_srcptr f, slong flen, slong n, slong prec)
void acb_poly_asinh_series(acb_poly_t res, const acb_poly_t f, slong n, slong prec)
void _acb_poly_acosh_series(acb_ptr res, acb_srcptr f, slong flen, slong n, slong prec)
void acb_poly_acosh_series(acb_poly_t res, const acb_poly_t f, slong n, slong prec)

Set `res` to the inverse hyperbolic tangent, sine or cosine of $f$
truncated to length `n`. The underscore versions allow aliasing.
Wrappers of the corresponding `_gr_poly` functions.
:::

::: function
void _acb_poly_exp_series_basecase(acb_ptr f, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_exp_series_basecase(acb_poly_t f, const acb_poly_t h, slong n, slong prec)
void _acb_poly_exp_series(acb_ptr f, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_exp_series(acb_poly_t f, const acb_poly_t h, slong n, slong prec)

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
void _acb_poly_exp_pi_i_series(acb_ptr f, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_exp_pi_i_series(acb_poly_t f, const acb_poly_t h, slong n, slong prec)

Set `f` to $\exp(\pi i h)$ truncated to length `n`, using
`acb_exp_pi_i` for the constant term (which is accurate when $h_0$ is
close to a rational number with small denominator). The underscore
version allows aliasing.
:::

::: function
void _acb_poly_sin_cos_series(acb_ptr s, acb_ptr c, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_sin_cos_series(acb_poly_t s, acb_poly_t c, const acb_poly_t h, slong n, slong prec)

Set `s` and `c` to $\sin h$ and $\cos h$ truncated to length `n`,
computed together. The underscore version allows aliasing. For $n \ge 3$
and nonconstant $h$ this wraps `_gr_poly_sin_cos_series`.
:::

::: function
void _acb_poly_sin_series(acb_ptr s, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_sin_series(acb_poly_t s, const acb_poly_t h, slong n, slong prec)
void _acb_poly_cos_series(acb_ptr c, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_cos_series(acb_poly_t c, const acb_poly_t h, slong n, slong prec)

Set `s` to $\sin h$ or `c` to $\cos h$ truncated to length `n`. The
underscore versions allow aliasing. Wrappers of `_gr_poly_sin_series`
and `_gr_poly_cos_series`.
:::

::: function
void _acb_poly_tan_series(acb_ptr g, acb_srcptr h, slong hlen, slong len, slong prec)
void acb_poly_tan_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec)

Set `g` to $\tan h$ truncated to length `n`. The underscore version
does not allow aliasing. Wrappers of `_gr_poly_tan_series`, which uses
Newton iteration on $\operatorname{atan}$ for long series.
:::

::: function
void _acb_poly_cot_series(acb_ptr g, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_cot_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec)

Set `g` to $\cot h$ truncated to length `n`. The underscore version
allows aliasing. Wrappers of `_gr_poly_cot_series`.
:::

::: function
void _acb_poly_sin_cos_pi_series(acb_ptr s, acb_ptr c, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_sin_cos_pi_series(acb_poly_t s, acb_poly_t c, const acb_poly_t h, slong n, slong prec)
void _acb_poly_sin_pi_series(acb_ptr s, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_sin_pi_series(acb_poly_t s, const acb_poly_t h, slong n, slong prec)
void _acb_poly_cos_pi_series(acb_ptr c, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_cos_pi_series(acb_poly_t c, const acb_poly_t h, slong n, slong prec)
void _acb_poly_tan_pi_series(acb_ptr g, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_tan_pi_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec)
void _acb_poly_cot_pi_series(acb_ptr c, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_cot_pi_series(acb_poly_t c, const acb_poly_t h, slong n, slong prec)

Compute the trigonometric functions of $\pi h$ truncated to length `n`.
These are more accurate than multiplying by $\pi$ first when the
constant term of $h$ is close to an integer or a half-integer. Wrappers
of the corresponding `_gr_poly` functions.
:::

::: function
void _acb_poly_sinh_cosh_series_basecase(acb_ptr s, acb_ptr c, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_sinh_cosh_series_basecase(acb_poly_t s, acb_poly_t c, const acb_poly_t h, slong n, slong prec)
void _acb_poly_sinh_cosh_series_exponential(acb_ptr s, acb_ptr c, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_sinh_cosh_series_exponential(acb_poly_t s, acb_poly_t c, const acb_poly_t h, slong n, slong prec)
void _acb_poly_sinh_cosh_series(acb_ptr s, acb_ptr c, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_sinh_cosh_series(acb_poly_t s, acb_poly_t c, const acb_poly_t h, slong n, slong prec)
void _acb_poly_sinh_series(acb_ptr s, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_sinh_series(acb_poly_t s, const acb_poly_t h, slong n, slong prec)
void _acb_poly_cosh_series(acb_ptr c, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_cosh_series(acb_poly_t c, const acb_poly_t h, slong n, slong prec)

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
void _acb_poly_tanh_series(acb_ptr g, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_tanh_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec)
void _acb_poly_coth_series(acb_ptr g, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_coth_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec)

Set `g` to $\tanh h$ or $\coth h$ truncated to length `n`. The
underscore versions allow aliasing. Wrappers of `_gr_poly_tanh_series`
and `_gr_poly_coth_series`.
:::

::: function
void _acb_poly_sinc_series(acb_ptr s, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_sinc_series(acb_poly_t s, const acb_poly_t h, slong n, slong prec)
void _acb_poly_sinc_pi_series(acb_ptr s, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_sinc_pi_series(acb_poly_t s, const acb_poly_t h, slong n, slong prec)

Set `s` to $\operatorname{sinc} h = \sin(h)/h$ or to
$\operatorname{sinc}(\pi h)$ truncated to length `n`. A constant term
that contains zero is handled through the Taylor series of
$\operatorname{sinc}$ on a ball around zero, so that the removable
singularity does not cause a division by zero.
:::

## Lambert W and gamma functions {#special}

::: function
void _acb_poly_lambertw_series(acb_ptr res, acb_srcptr z, slong zlen, const fmpz_t k, int flags, slong len, slong prec)
void acb_poly_lambertw_series(acb_poly_t res, const acb_poly_t z, const fmpz_t k, int flags, slong len, slong prec)

Set `res` to branch $k$ of the Lambert W function of the power series
$z$, truncated to length `len`. The value at the constant term is
computed by `acb_lambertw` with the same $k$ and `flags` (which select
the branch cut conventions, see there) [@Joh2017b], and the remaining
coefficients by Newton iteration on $w e^w = z$. The underscore version
allows aliasing.
:::

::: function
void _acb_poly_gamma_series(acb_ptr res, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_gamma_series(acb_poly_t res, const acb_poly_t h, slong n, slong prec)
void _acb_poly_rgamma_series(acb_ptr res, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_rgamma_series(acb_poly_t res, const acb_poly_t h, slong n, slong prec)
void _acb_poly_lgamma_series(acb_ptr res, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_lgamma_series(acb_poly_t res, const acb_poly_t h, slong n, slong prec)
void _acb_poly_digamma_series(acb_ptr res, acb_srcptr h, slong hlen, slong n, slong prec)
void acb_poly_digamma_series(acb_poly_t res, const acb_poly_t h, slong n, slong prec)

Set `res` to $\Gamma(h)$, $1/\Gamma(h)$, $\log \Gamma(h)$ or $\psi(h)$
truncated to length `n`. The underscore versions allow aliasing.

::: details
If all coefficients of $h$ are real, the corresponding `arb_poly`
function is used. Otherwise the Taylor series at the constant term is
computed from Stirling's series after shifting the argument by a rising
factorial (and using the reflection formula in the left half-plane),
with the parameters of `acb_hypgeom_gamma_stirling_choose_param`, and
composed with $h - h_0$
([](algorithms/special-functions/gamma)). `acb_poly_lgamma_series`
gives the principal branch of the log-gamma function, continuous away
from the negative real axis.
:::
:::

::: function
void _acb_poly_rising_ui_series(acb_ptr res, acb_srcptr f, slong flen, ulong r, slong trunc, slong prec)
void acb_poly_rising_ui_series(acb_poly_t res, const acb_poly_t f, ulong r, slong trunc, slong prec)

Set `res` to the rising factorial $f (f+1) \cdots (f+r-1)$ truncated to
length `trunc`, by binary splitting. The underscore version requires
`flen`, $r$ and `trunc` to be positive and no aliasing.
:::

## Power sums {#power-sums}

The functions in this section compute
$$S(s, a, q, n) = \sum_{k=0}^{n-1} \frac{q^k}{(k+a)^{s+t}}$$
as a power series in $t$ truncated to length `len`, the main part of
the Euler–Maclaurin evaluation of the Hurwitz zeta function and of the
direct summation of polylogarithms.

::: function
void _acb_poly_powsum_series_naive(acb_ptr z, const acb_t s, const acb_t a, const acb_t q, slong n, slong len, slong prec)
void _acb_poly_powsum_series_naive_threaded(acb_ptr z, const acb_t s, const acb_t a, const acb_t q, slong n, slong len, slong prec)

Set `{z, len}` to $S(s, a, q, n)$ by summing the terms one by one; the
series of $(k+a)^{-s-t}$ is $(k+a)^{-s} \sum_j (-\log(k+a))^j t^j/j!$.
The *threaded* version splits the sum over `flint_get_num_threads()`
threads.
:::

::: function
void _acb_poly_powsum_series_tree(acb_ptr z, const acb_t s, const acb_t a, const acb_t q, slong n, slong len, slong prec)

Sets `{z, len}` to $S(s, a, q, n)$ using a transposed multipoint
evaluation, in $\tilde O(\max(n, \mathit{len}))$ arithmetic operations
instead of $O(n \cdot \mathit{len})$.

::: details
With $L_k = \log(a+k)$, $w_k = q^k (a+k)^{-s}$ and $b_k = c - L_k$ for
a suitable centre $c$, we have $\sum_k w_k e^{-L_k t} = e^{-ct} \sum_j
p_j t^j / j!$, where the weighted power sums $p_j = \sum_k w_k b_k^j$
are the coefficients of the rational function $\sum_k w_k / (1 - b_k t)
= P(t)/Q(t)$. The numerator and denominator are computed with a
subproduct tree whose nodes are truncated to length `len`, followed by
one power series division. When $n$, `len` and `prec` are all
$\Theta(m)$, this costs $\tilde O(m^2)$ bit operations instead of
$O(m^3)$.

The algorithm is numerically unstable, so the working precision is
increased by the number of guard bits returned by
`_acb_poly_powsum_series_tree_guard_bits`, which is chosen so that the
output radii are comparable to those of the naive algorithm. The
computation is done at the midpoints of $s$, $a$ and $q$, and the error
propagated from inexact input is bounded separately term by term. It
requires $\operatorname{Re}(a) > 0$ and $q$ bounded away from zero, and
otherwise falls back to the naive algorithm. The tree is evaluated in
parallel when `flint_get_num_threads()` is greater than one.
:::
:::

::: function
slong _acb_poly_powsum_series_tree_guard_bits(const acb_t s, const acb_t a, const acb_t q, const acb_t c, slong n, slong len)

Returns the number of guard bits used by `_acb_poly_powsum_series_tree`
for the exact parameters $s$, $a$, $q$ and the centre $c$: roughly a
constant times $\min(n, \mathit{len})$ plus a term for the variation of
the weights, estimated in double precision from a sample of the terms.
:::

::: function
int _acb_poly_powsum_series_tree_is_faster(const acb_t s, const acb_t a, const acb_t q, slong n, slong len, slong prec, int naive_threaded)

Returns nonzero if `_acb_poly_powsum_series_tree` is predicted to be
faster than `_acb_poly_powsum_series_naive` (or the threaded naive
version if `naive_threaded` is set) for these parameters, using
empirical cost models and the guard bit estimate. This selects the
algorithm in the Euler–Maclaurin summation of the Hurwitz zeta function
and in the direct summation of the polylogarithm for small $|z|$.
:::

::: function
void _acb_poly_powsum_one_series_sieved(acb_ptr z, const acb_t s, slong n, slong len, slong prec)

Sets `{z, len}` to $\sum_{k=1}^{n} k^{-(s+t)}$ as a power series in $t$.
Powers of composite $k = ij$ are computed as products $i^{-(s+t)}
j^{-(s+t)}$ of stored powers, and the even terms are grouped as a
polynomial in $2^{-(s+t)}$. This takes about $n / \log n$ powers, $n/2$
multiplications and storage for $n/6$ power series, and is faster than
the naive algorithm only when `len` is small.
:::

## Zeta function {#zeta}

The functions with `zeta_em` in their name implement Euler–Maclaurin
summation for the Hurwitz zeta function and its derivatives, with the
error bounds of [](algorithms/special-functions/zeta#euler-maclaurin).

::: function
void _acb_poly_zeta_em_choose_param(mag_t bound, ulong * N, ulong * M, const acb_t s, const acb_t a, slong d, slong target, slong prec)

Chooses the parameters $N$ (number of terms summed directly) and $M$
(number of Bernoulli terms) of the Euler–Maclaurin formula for the
length-$d$ series of $\zeta(s, a)$, aiming at an error of about
$2^{-\mathit{target}}$ times an estimate $|a^{-s}|$ of the magnitude
(absolute error $2^{-\mathit{target}}$ for $a = 1$), and sets `bound`
to the corresponding error bound.
:::

::: function
void _acb_poly_zeta_em_bound1(mag_t bound, const acb_t s, const acb_t a, slong N, slong M, slong d, slong wp)
void _acb_poly_zeta_em_bound(arb_ptr vec, const acb_t s, const acb_t a, ulong N, ulong M, slong d, slong wp)

Compute bounds for the error of the Euler–Maclaurin evaluation with
parameters $N$, $M$: `_acb_poly_zeta_em_bound` bounds each of the $d$
coefficients, `_acb_poly_zeta_em_bound1` gives a single bound for all
of them. The formulas are those of [@Joh2013].
:::

::: function
void _acb_poly_zeta_em_tail_naive(acb_ptr z, const acb_t s, const acb_t Na, acb_srcptr Nasx, slong M, slong len, slong prec)
void _acb_poly_zeta_em_tail_bsplit(acb_ptr z, const acb_t s, const acb_t Na, acb_srcptr Nasx, slong M, slong len, slong prec)

Evaluate the tail (the Bernoulli number terms) of the Euler–Maclaurin
formula for $\zeta(s + x, a)$ as a power series, given $N + a$ and the
series of $(N + a)^{-(s+x)}$, by the naive recurrence or by binary
splitting.
:::

::: function
void _acb_poly_zeta_em_sum(acb_ptr z, const acb_t s, const acb_t a, int deflate, ulong N, ulong M, slong d, slong prec)

Sets `{z, d}` to the truncated Euler–Maclaurin sum of order $N$, $M$
for the series of $\zeta(s + x, a)$ (without the error term). If
`deflate` is nonzero, $\zeta(s + x, a) - 1/(s + x - 1)$ is computed,
which permits expansion at $s = 1$.
:::

::: function
void _acb_poly_zeta_cpx_series(acb_ptr z, const acb_t s, const acb_t a, int deflate, slong d, slong prec)

Sets `{z, d}` to the series of $\zeta(s + x, a)$ (or the deflated
function if `deflate` is nonzero), choosing the Euler–Maclaurin
parameters with `_acb_poly_zeta_em_choose_param` for the target `prec`
and adding the error bound to each coefficient
([](algorithms/special-functions/zeta)).
:::

::: function
void _acb_poly_zeta_series(acb_ptr res, acb_srcptr h, slong hlen, const acb_t a, int deflate, slong len, slong prec)
void acb_poly_zeta_series(acb_poly_t res, const acb_poly_t f, const acb_t a, int deflate, slong n, slong prec)

Set `res` to the Hurwitz zeta function $\zeta(s, a)$ of the power
series $s$ = `f` for a constant $a$, truncated to length `n`; $a = 1$
gives the Riemann zeta function. If `deflate` is nonzero, $\zeta(s, a)
+ 1/(1 - s)$ is computed, which is finite at $s = 1$; with $s = 1 + x$
this gives the Stieltjes constants, $\sum_{k=0}^{n-1} \frac{(-1)^k}{k!}
\gamma_k(a) x^k$.

::: details
The series at the constant term is computed by `acb_dirichlet_zeta_jet`
if $a = 1$ and by `_acb_poly_zeta_cpx_series` otherwise, and then
composed with the nonconstant part of $s$. `acb_dirichlet_zeta_jet` uses
the functional equation if the real part of the midpoint is negative,
the Riemann–Siegel formula for at most two coefficients at very large
imaginary parts (above about $60\,\mathit{prec}^{3/2}$), and
Euler–Maclaurin summation otherwise.
:::
:::

## Other special functions {#other-special}

::: function
void _acb_poly_polylog_cpx_small(acb_ptr w, const acb_t s, const acb_t z, slong len, slong prec)
void _acb_poly_polylog_cpx_zeta(acb_ptr w, const acb_t s, const acb_t z, slong len, slong prec)
void _acb_poly_polylog_cpx(acb_ptr w, const acb_t s, const acb_t z, slong len, slong prec)

Set `{w, len}` to the Taylor series in $x$ of the polylogarithm
$\operatorname{Li}_{s+x}(z)$ for complex constants $s$ and $z$; `len`
must be positive and `w` must not alias `s` or `z`. The *small* version
sums the power series in $z$ (for $|z| < 1$), the *zeta* version writes
the polylogarithm as a combination of two Hurwitz zeta functions, and
the default uses the *small* version for $|z| < 1/2$ and the *zeta*
version otherwise (and `acb_hypgeom_dilog` for $s = 2$ and length 1)
([](algorithms/special-functions/zeta)).
:::

::: function
void _acb_poly_polylog_series(acb_ptr w, acb_srcptr s, slong slen, const acb_t z, slong len, slong prec)
void acb_poly_polylog_series(acb_poly_t w, const acb_poly_t s, const acb_t z, slong len, slong prec)

Set `w` to $\operatorname{Li}_s(z)$ for a power series $s$ and a
constant $z$, truncated to length `len`. The underscore version requires
positive lengths and allows aliasing.
:::

::: function
void _acb_poly_erf_series(acb_ptr res, acb_srcptr z, slong zlen, slong n, slong prec)
void acb_poly_erf_series(acb_poly_t res, const acb_poly_t z, slong n, slong prec)

Set `res` to the error function of the power series $z$, truncated to
length `n`. These wrappers of `acb_hypgeom_erf_series` are kept for
backwards compatibility; see also `acb_hypgeom_erfc_series` and
`acb_hypgeom_erfi_series`.
:::

::: function
void _acb_poly_agm1_series(acb_ptr res, acb_srcptr z, slong zlen, slong len, slong prec)
void acb_poly_agm1_series(acb_poly_t res, const acb_poly_t z, slong n, slong prec)

Set `res` to the arithmetic-geometric mean $M(1, z)$ of 1 and the power
series $z$, truncated to length `n`: the derivatives at the constant
term come from `acb_agm1_cpx`, and the series is composed with the
nonconstant part of $z$.
:::

::: function
void _acb_poly_elliptic_k_series(acb_ptr res, acb_srcptr z, slong zlen, slong len, slong prec)
void acb_poly_elliptic_k_series(acb_poly_t res, const acb_poly_t z, slong n, slong prec)
void _acb_poly_elliptic_p_series(acb_ptr res, acb_srcptr z, slong zlen, const acb_t tau, slong len, slong prec)
void acb_poly_elliptic_p_series(acb_poly_t res, const acb_poly_t z, const acb_t tau, slong n, slong prec)

Wrappers, kept for backwards compatibility, of `acb_elliptic_k_series`
(the complete elliptic integral $K(z)$) and `acb_elliptic_p_series`
(the Weierstrass function $\wp(z, \tau)$) in [](acb_elliptic), which
should be used instead.
:::

## Root finding {#roots}

::: function
void _acb_poly_root_bound_fujiwara(mag_t bound, acb_srcptr poly, slong len)
void acb_poly_root_bound_fujiwara(mag_t bound, acb_poly_t poly)

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
void _acb_poly_root_inclusion(acb_t r, const acb_t m, acb_srcptr poly, acb_srcptr polyder, slong len, slong prec)

Given a point $m$, a nonconstant polynomial $f$ = `{poly, len}` of
degree $n$ and its derivative `polyder`, sets `r` to a box centred at
$m$ that contains a disk of radius $n |f(m)/f'(m)|$, which contains at
least one root of every polynomial in the ball $f$.

::: details
If no root $\zeta_i$ of $f$ lies within distance $\rho = n|f(m)/f'(m)|$
of $m$, then
$$\left|\frac{f'(m)}{f(m)}\right| = \left|\sum_{i=1}^{n} \frac{1}{m - \zeta_i}\right|
\le \sum_{i=1}^{n} \frac{1}{|m - \zeta_i|} < \frac{n}{\rho} = \left|\frac{f'(m)}{f(m)}\right|,$$
a contradiction [@Kob2010]. The values $f(m)$ and $f'(m)$ are evaluated
by Horner's rule with an error bound for the disk (rather than the
rectangle) containing each partial result, which matters for high
degree, and upper bounds for $|f(m)|$ and $1/|f'(m)|$ over the
coefficient balls are used, so the disk is valid for every polynomial in
the ball. If $f(m)$ is exactly zero, the radius is zero.
:::
:::

::: function
slong _acb_poly_validate_roots(acb_ptr roots, acb_srcptr poly, slong len, slong prec)

Given approximations `{roots, len - 1}` of the roots of `{poly, len}`,
replaces each by an inclusion box from `_acb_poly_root_inclusion`
around its midpoint, moves the boxes that do not overlap any other box
to the front, and returns their number. Every box contains a root; if
the return value equals the degree, each box contains exactly one root
and all roots are isolated.
:::

::: function
slong _acb_poly_find_roots(acb_ptr roots, acb_srcptr poly, acb_srcptr initial, slong len, slong maxiter, slong prec)
slong acb_poly_find_roots(acb_ptr roots, const acb_poly_t poly, acb_srcptr initial, slong maxiter, slong prec)

Attempts to compute all $n$ roots of the nonzero polynomial `poly` of
degree $n$, writing $n$ enclosures to `roots` and returning the number of
roots that are isolated; these come first. If the return value equals
$n$, every enclosure contains exactly one root. `initial` gives starting
values (or `NULL`), and `maxiter` limits the number of iterations (0
selects a default). The polynomial should be squarefree.

::: details
The roots are approximated numerically by the Weierstrass–Durand–Kerner
iteration in floating-point arithmetic, starting from points chosen
with the Newton polygon, and the approximations are then validated by
`_acb_poly_validate_roots`. See [](acb_poly/find_roots).
:::
:::

::: function
double _acb_poly_find_roots_double(acb_ptr roots, acb_srcptr poly, acb_srcptr initial, slong len, slong maxiter, slong prec)

Sets `{roots, len - 1}` to approximations of the roots of `{poly, len}`
computed by `cd_poly_find_roots` in double precision after rescaling
$x$ and the coefficients so that the coefficients are representable as
doubles, with the radii set to zero. Returns the largest relative step
of the last iteration. `prec` (if below 53) only serves as an early
termination criterion.

::: details
The scaling factors are chosen as the line that best fits the points
$(i, \log_2 |a_i|)$ in the minimax sense, computed with the convex hulls
of these points. Barring overflow and up to rounding errors, the
Gershgorin circle theorem places the roots in the union of the disks
centred at the outputs $z_k$ with radius $(\mathit{len} - 2) |z_k|$
times the returned correction.
:::
:::

::: function
int _acb_poly_validate_real_roots(acb_srcptr roots, acb_srcptr poly, slong len, slong prec)
int acb_poly_validate_real_roots(acb_srcptr roots, const acb_poly_t poly, slong prec)

Given a real polynomial `poly` (stored as an `acb_poly_t` with zero
imaginary parts) and isolating enclosures of all its complex roots,
returns nonzero if it can prove that the enclosures whose imaginary part
contains zero contain real roots. Then every such root is real (its
imaginary part can be set to zero), and every other enclosure contains a
nonreal root. Returns zero if this cannot be decided at the accuracy of
the input and the precision `prec`.

::: details
Let $r$ be the number of enclosures that meet the real axis. The test
requires $r \equiv \deg f \pmod 2$ (nonreal roots come in conjugate
pairs) and checks, at the midpoint between consecutive candidates (in
the order of their midpoints), that $f$ has the signs that alternate
from the sign at $-\infty$ to the sign at $+\infty$. Then $f$ changes
sign within each candidate interval, which therefore contains a real
root, and since each enclosure contains exactly one root, that root is
real.
:::
:::

## Root finding in double precision {#cd-poly}

These functions implement the Weierstrass–Durand–Kerner iteration for
complex polynomials with `double` coefficients, used by
`_acb_poly_find_roots` at precisions up to 53 bits. A polynomial of
length $n$ is an array of $2n$ doubles, the real and imaginary parts of
the coefficient of $x^k$ being at indices $2k$ and $2k + 1$; vectors of
$n - 1$ roots are stored in the same way. The internal functions work
on separate arrays `z_r`, `z_i` of real and imaginary parts.

::: function
double cd_poly_find_roots(double * roots, const double * poly, const double * initial, slong len, slong maxiter, double reltol)

Sets `roots` (room for $2\,\mathit{len}$ doubles) to approximations of
the $\mathit{len} - 1$ roots of `poly` and returns the largest relative
step $|\Delta z_k| / |z_k|$ of the last iteration (at least
$2^{-50}$). `initial` is `NULL` or holds $\mathit{len} - 1$ starting
values (zero entries are replaced by values from the Newton polygon).

::: details
Zero leading coefficients give roots at infinity (set to $2^{1023}(1 +
i)$) and zero trailing coefficients roots at zero. At most `maxiter`
iterations are done; the iteration stops when $d - 1$ times the
largest relative step is below `reltol`, where $d$ is the number of
finite nonzero roots, or when the step is below $2^{-50}$. Each step is
limited to an absolute size between $2^{-10}$ and $2^{-3}$, adapted to
the progress (the callers scale the polynomial so that the roots have
modulus near 1), and roots of modulus greater than 1 are updated
through the reversed polynomial to avoid overflow. The output is
sorted.
:::
:::

::: function
double cd_poly_refine_roots(double * z, const double * poly, slong len, double stepsize_bound)
double cd_poly_refine_roots_with_pivot(double * z, const double * poly, slong len, double stepsize_bound)

Update the $\mathit{len} - 1$ approximate roots in `z` with one damped
Weierstrass–Durand–Kerner step $z_k \leftarrow z_k - \delta_k$, where
$\delta_k = f(z_k) / (a_n \prod_{j \ne k} (z_k - z_j))$ is scaled down so
that $|\delta_k| \le$ `stepsize_bound`, and return the largest relative
step. The version with pivot evaluates large roots through the reversed
polynomial and orders the factors of the products to reduce the risk of
overflow; it does not preserve the order of the roots.
:::

::: function
void cd_poly_roots_initial_values(double * z_r, double * z_i, const double * p_r, const double * p_i, slong n, const double * z0, slong d)

Sets `z_r`, `z_i` to $n - 1$ initial values for the roots of the
polynomial with coefficients `p_r`, `p_i` of length $n$: the nonzero
entries of `z0` (of length $d$, or `NULL`), completed by points on
circles given by the Newton polygon of $(k, -\log|p_k|)$ [@Bin1996].
:::

::: function
void cd_poly_horner(double * results_r, double * results_i, const double * values_r, const double * values_i, slong n_start, slong n_end, const double * coefficients_r, const double * coefficients_i, slong n)

Evaluates the polynomial of length $n$ with the given coefficients at
the points with indices $n_{\mathrm{start}} \le k < n_{\mathrm{end}}$,
by Horner's rule vectorised over the points.
:::

::: function
void cd_poly_weierstrass(double * results_r, double * results_i, double lc_r, double lc_i, const double * values_r, const double * values_i, slong n_start, slong n_end, slong d)
void cd_poly_weierstrass_distinct_orders(double * results_r, double * results_i, double lc_r, double lc_i, const double * col_values_r, const double * col_values_i, slong d, const double * row_values_r, const double * row_values_i, slong n_start, slong n_end)

Compute the Weierstrass denominators $w_k = c \prod_{j \ne k} (z_k -
z_j)$, with $c$ = `lc_r + i lc_i`, for $n_{\mathrm{start}} \le k <
n_{\mathrm{end}}$. `cd_poly_weierstrass` uses the $d$ points `values`
for both factors and exploits the symmetry of the differences;
`cd_poly_weierstrass_distinct_orders` takes the factors $z_j$ from
`col_values` (in an order chosen to avoid overflow) and the points
$z_k$ from `row_values`, skipping a factor that is exactly zero.
:::

::: function
double cd_poly_wdk_update(double * z_r, double * z_i, const double * vp_r, const double * vp_i, const double * wdk_r, const double * wdk_i, slong n_start, slong n_end, double stepsize_bound)

Applies the damped updates $z_k \leftarrow z_k - \delta_k$ with
$\delta_k = v_k / w_k$ (values `vp`, denominators `wdk`, using Smith's
division), scaled so that $|\delta_k| \le$ `stepsize_bound`, for
$n_{\mathrm{start}} \le k < n_{\mathrm{end}}$, and returns the largest
relative step $|\delta_k|/|z_k|$ before damping.
:::

## Internal functions {#internal}

::: function
void _acb_poly_binomial_pow_acb_series(acb_ptr h, acb_srcptr f, slong flen, const acb_t g, slong len, slong prec)

Sets `{h, len}` to $(a + bx^d)^g$ truncated to length `len`, where `{f,
flen}` is $a + bx^d$ with $d = \mathit{flen} - 1$ and zero coefficients
in between (not checked), using the binomial series. No aliasing.
:::

::: function
void _acb_poly_acb_pow_cpx(acb_ptr w, const acb_t a, const acb_t b, slong len, slong prec)

Sets `{w, len}` to the power series $a^{b + x} = a^b \sum_k (\log
a)^k x^k / k!$ truncated to length `len`. Also available under the
name `_acb_poly_pow_cpx`.
:::

::: function
void _acb_poly_acb_invpow_cpx(acb_ptr res, const acb_t N, const acb_t c, slong trunc, slong prec)

Sets `{res, trunc}` to the power series $N^{-(c + x)}$ truncated to
length `trunc`.
:::

::: function
void _acb_poly_mullow_cpx(acb_ptr res, acb_srcptr src, slong len, const acb_t c, slong trunc, slong prec)

Sets `res` to the product of `{src, len}` and $c + x$, truncated to
length `trunc` (which must be `len` or `len + 1`). Aliasing of `res` and
`src` is allowed.
:::
