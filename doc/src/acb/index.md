---
title: acb.h
header: acb.h
description: Complex numbers as rectangles of real balls, with rigorous arithmetic and elementary and special functions
category: Real and complex numbers
order: 4
---

An `acb_t` is a complex number represented by a pair of real balls
(`arb_t`), its real and imaginary parts. It stands for the rectangle

$$
[m_1 \pm r_1] + [m_2 \pm r_2] i
$$

in the complex plane, and every function returns a rectangle that
contains the exact result for all points of the input rectangles. A
rectangle rather than a disk is used because most operations split
naturally into operations on the real and imaginary parts, and because it
can represent numbers whose real or imaginary part is exactly known (for
example exactly real numbers) while the other part is not.

The module provides arithmetic, elementary functions with the standard
branch cuts, the gamma, zeta and polylogarithm functions, the
arithmetic-geometric mean, the Lambert W function, and analytic
extensions of piecewise real functions for use in integration. Many more
special functions are in [acb_hypgeom](acb_hypgeom),
[acb_elliptic](acb_elliptic), [acb_modular](acb_modular) and
[acb_dirichlet](acb_dirichlet), and polynomials, power series and
matrices in [acb_poly](acb_poly) and [acb_mat](acb_mat) (Taylor
expansions with power series are explained in
[](arb_poly/power-series)). The semantics
of precision, predicates and printing are those of [arb](arb), applied to
each part; [](acb/getting-started) introduces the module and
[](acb/branch-cuts) explains how the multivalued functions are defined.
The tutorials [](tutorials/numerics/integrals) and
[](tutorials/numerics/polynomial-roots) use complex balls in complete
programs.

::: summary algorithms/ball-arithmetic/overview :::

## Types and conventions {#types}

::: type
acb_struct
acb_t
acb_ptr
acb_srcptr

An `acb_struct` has two `arb_struct` fields, the real and the imaginary
part. `acb_t` is an array of one `acb_struct`, so that it is passed by
reference; `acb_ptr` and `acb_srcptr` point to arrays of complex numbers.
:::

::: macro
acb_realref(x)
acb_imagref(x)

Pointers to the real and the imaginary part of `x` as `arb_t`, for use
with the functions of `arb`.
:::

::: function
arb_ptr acb_real_ptr(acb_t z)
arb_ptr acb_imag_ptr(acb_t z)

Function versions of `acb_realref` and `acb_imagref`, for language
bindings that cannot use macros.
:::

The conventions of [arb](arb) apply: functions with a `prec` argument
round their results to `prec` bits and include the rounding error in the
radii; outputs may be aliased with inputs; identical input pointers allow
algebraic simplification (for example `acb_mul(z, x, x, prec)` uses the
squaring formula); and predicates return 1 only when the property
certainly holds. A function that is not defined, or not analytic, on part
of its input returns an *indeterminate* number, whose real and imaginary
parts are $[\operatorname{NaN} \pm \infty]$. The branch cuts of the
multivalued functions are those of [@NIST2012]; see [](acb/branch-cuts).

## Memory management {#memory}

::: function
void acb_init(acb_t x)

Initialises `x` to the exact value 0.
:::

::: function
void acb_clear(acb_t x)

Frees the memory used by `x`.
:::

::: function
acb_ptr _acb_vec_init(slong n)
void _acb_vec_clear(acb_ptr v, slong n)

Allocate a vector of `n` complex numbers initialised to zero, and free
such a vector.
:::

::: function
void acb_swap(acb_t z, acb_t x)

Swaps `z` and `x` without copying.
:::

::: function
slong acb_allocated_bytes(const acb_t x)
slong _acb_vec_allocated_bytes(acb_srcptr vec, slong len)
double _acb_vec_estimate_allocated_bytes(slong len, slong prec)

Return the heap memory used by `x` (excluding `sizeof(acb_struct)`), the
total memory of a vector, and an estimate of the memory needed for a
vector of `len` numbers at precision `prec`; see
`_arb_vec_estimate_allocated_bytes`.
:::

## Assignment and rounding {#assignment}

::: function
void acb_zero(acb_t z)
void acb_one(acb_t z)
void acb_onei(acb_t z)

Set `z` to 0, 1 and $i$ respectively.
:::

::: function
void acb_indeterminate(acb_t x)

Sets `x` to $[\operatorname{NaN} \pm \infty] + [\operatorname{NaN} \pm
\infty] i$, an indeterminate value.
:::

::: function
void acb_set(acb_t z, const acb_t x)
void acb_set_ui(acb_t z, ulong c)
void acb_set_si(acb_t z, slong c)
void acb_set_d(acb_t z, double c)
void acb_set_fmpz(acb_t z, const fmpz_t c)
void acb_set_arb(acb_t z, const arb_t c)

Set `z` to the exact value of the input, with zero imaginary part for a
real input.
:::

::: function
void acb_set_si_si(acb_t z, slong x, slong y)
void acb_set_d_d(acb_t z, double x, double y)
void acb_set_fmpz_fmpz(acb_t z, const fmpz_t x, const fmpz_t y)
void acb_set_arb_arb(acb_t z, const arb_t x, const arb_t y)

Set `z` to $x + yi$ exactly.
:::

::: function
void acb_set_fmpq(acb_t z, const fmpq_t c, slong prec)
void acb_set_round(acb_t z, const acb_t x, slong prec)
void acb_set_round_fmpz(acb_t z, const fmpz_t y, slong prec)
void acb_set_round_arb(acb_t z, const arb_t x, slong prec)

Set `z` to a rectangle containing the input, with midpoints rounded to
`prec` bits.
:::

::: function
void acb_get_mid(acb_t res, const acb_t x)

Sets `res` to the exact midpoint of `x` (radii zero).
:::

::: function
void acb_add_error_arf(acb_t x, const arf_t err)
void acb_add_error_mag(acb_t x, const mag_t err)
void acb_add_error_arb(acb_t x, const arb_t err)

Add an upper bound for $|\mathit{err}|$ to the radii of both the real and
the imaginary part of `x`, so that `x` then contains every point within
distance $|\mathit{err}|$ of the old `x` in each coordinate.
:::

::: function
void acb_trim(acb_t z, const acb_t x)

Sets `z` to `x` with both parts trimmed by `arb_trim`.
:::

## Input and output {#io}

::: function
void acb_printn(const acb_t x, slong digits, ulong flags)
void acb_fprintn(FILE * fp, const acb_t z, slong digits, ulong flags)

Print `x` as `a + b*I` with the real and imaginary parts in the decimal
ball format of `arb_get_str` (omitting a zero part), to standard output
or to `fp`. The printed rectangle contains `x`, and the flags of
`arb_get_str` apply to both parts. This is the recommended way to print a
complex number.
:::

::: function
void acb_printd(const acb_t z, slong digits)
void acb_fprintd(FILE * file, const acb_t z, slong digits)

Print approximate decimal values of the midpoints with `digits` digits
and of the radii with 3 digits, in the form `(a + bj)  +/-  (r, sj)`. The
printed values are not guaranteed to enclose `z`.
:::

::: function
void acb_print(const acb_t x)
void acb_fprint(FILE * file, const acb_t x)

Print the exact binary representation of both parts, as `(a, b)` in the
format of `arb_print`.
:::

## Real and imaginary parts, absolute value and argument {#parts}

::: function
void acb_get_real(arb_t re, const acb_t z)
void acb_get_imag(arb_t im, const acb_t z)

Set `re` to the real part, respectively `im` to the imaginary part of `z`.
:::

::: function
void acb_abs(arb_t u, const acb_t z, slong prec)

Sets `u` to $|z|$, computed with `arb_hypot`.
:::

::: function
void acb_arg(arb_t r, const acb_t z, slong prec)

Sets `r` to the argument (phase) of `z`, $\operatorname{atan2}(b, a)$
for $z = a + bi$, with values in $(-\pi, \pi]$ and the branch cut on
$(-\infty, 0]$: $\arg(0) = 0$ and $\arg(a) = \pi$ for $a < 0$. A rectangle
that crosses the negative real axis gives an enclosure of $[-\pi, \pi]$.
:::

::: function
void acb_sgn(acb_t res, const acb_t z, slong prec)

Sets `res` to the complex sign $z / |z|$, and to 0 if $z$ is exactly 0.
:::

::: function
void acb_csgn(arb_t res, const acb_t z)

Sets `res` to the sign of the real part, extended to the imaginary axis
by the sign of the imaginary part: 1 in the right half plane, $-1$ in
the left half plane, $\operatorname{sgn}(\operatorname{Im} z)$ on the
imaginary axis, and 0 at $z = 0$. Equivalently $\operatorname{csgn}(z) =
z/\sqrt{z^2}$ for $z \ne 0$.
:::

## Bounds and accuracy {#bounds}

::: function
void acb_get_mag(mag_t z, const acb_t x)
void acb_get_mag_lower(mag_t z, const acb_t x)

Set `z` to an upper, respectively lower, bound for $|t|$ over all $t$ in
the rectangle `x`.
:::

::: function
void acb_get_abs_ubound_arf(arf_t u, const acb_t z, slong prec)
void acb_get_abs_lbound_arf(arf_t u, const acb_t z, slong prec)

Set `u` to an upper, respectively lower, bound for $|t|$, $t \in z$,
computed at precision `prec`.
:::

::: function
void acb_get_rad_ubound_arf(arf_t u, const acb_t z, slong prec)

Sets `u` to an upper bound for the distance from the midpoint to any
point of the rectangle, computed (not tightly) from the two radii.
:::

::: function
slong acb_rel_error_bits(const acb_t x)
slong acb_rel_accuracy_bits(const acb_t x)

Return the relative error, respectively accuracy, in bits, as
`arb_rel_error_bits` would for a real ball whose midpoint is the larger
of the two midpoints and whose radius is the larger of the two radii. An
exact `x` has accuracy `ARF_PREC_EXACT`, and an `x` with zero midpoint or
a non-finite part $-$`ARF_PREC_EXACT`.
:::

::: function
slong acb_rel_one_accuracy_bits(const acb_t x)

Returns the relative accuracy in bits of the rectangle with midpoint
magnitude replaced by $\max(1, |m|)$, measuring accuracy relative to 1
near the origin.
:::

::: function
slong acb_bits(const acb_t x)

Returns the larger of `arb_bits` of the two parts, the precision needed to
represent the midpoints exactly.
:::

::: function
void acb_union(acb_t res, const acb_t x, const acb_t y, slong prec)

Sets `res` to a rectangle containing both `x` and `y`.
:::

## Comparison and predicates {#comparison}

::: function
int acb_is_zero(const acb_t z)
int acb_is_one(const acb_t z)
int acb_equal_si(const acb_t x, slong y)

Return whether `z` is exactly 0, exactly 1, respectively exactly the
integer `y`.
:::

::: function
int acb_is_exact(const acb_t z)
int acb_is_finite(const acb_t x)
int acb_is_real(const acb_t x)

Return whether both radii are zero, whether both parts are finite, and
whether the imaginary part is exactly zero (without testing the real
part).
:::

::: function
int acb_is_int(const acb_t z)
int acb_is_int_2exp_si(const acb_t z, slong e)

Return whether `z` is an exact integer, respectively an exact integer
multiple of $2^e$ (with zero imaginary part).
:::

::: function
int acb_eq(const acb_t x, const acb_t y)
int acb_ne(const acb_t x, const acb_t y)

`acb_eq` returns 1 if `x` and `y` are certainly equal (`arb_eq` holds for
both parts), and `acb_ne` returns 1 if they are certainly different
(`arb_ne` holds for at least one part).
:::

::: function
int acb_equal(const acb_t x, const acb_t y)

Returns whether `x` and `y` are identical as rectangles, both parts being
equal as balls (`arb_equal`). This is a test of the representation, not
a proof that two numbers are equal.
:::

::: function
int acb_overlaps(const acb_t x, const acb_t y)

Returns whether the rectangles `x` and `y` have a point in common.
:::

::: function
int acb_contains(const acb_t x, const acb_t y)
int acb_contains_fmpq(const acb_t x, const fmpq_t y)
int acb_contains_fmpz(const acb_t x, const fmpz_t y)

Return whether `y` is contained in the rectangle `x`.
:::

::: function
int acb_contains_zero(const acb_t x)
int acb_contains_int(const acb_t x)

Return whether `x` contains zero, respectively an integer.
:::

::: function
int acb_contains_interior(const acb_t x, const acb_t y)

Returns whether `y` is contained in the interior of `x`.

::: details
The interior is that of the rectangle in the plane, so the test is always
false if `x` has an exact real or imaginary part: an imaginary part 0 is
not in the interior of the point interval 0. In particular it is false
for two real numbers. Such cases must be handled separately where a
different interpretation is intended.
:::
:::

::: function
int acb_get_unique_fmpz(fmpz_t z, const acb_t x)

If `x` contains exactly one integer (so its imaginary part contains 0
and its real part a unique integer), sets `z` to it and returns nonzero;
otherwise returns zero.
:::

## Arithmetic {#arithmetic}

[](acb/arithmetic) describes the error bounds and special
cases of these functions.

::: function
void acb_neg(acb_t z, const acb_t x)
void acb_neg_round(acb_t z, const acb_t x, slong prec)
void acb_conj(acb_t z, const acb_t x)

Set `z` to $-x$ (exactly, respectively rounded to `prec` bits) or to the
complex conjugate $\bar x$.
:::

::: function
void acb_add(acb_t z, const acb_t x, const acb_t y, slong prec)
void acb_add_ui(acb_t z, const acb_t x, ulong c, slong prec)
void acb_add_si(acb_t z, const acb_t x, slong c, slong prec)
void acb_add_fmpz(acb_t z, const acb_t x, const fmpz_t y, slong prec)
void acb_add_arb(acb_t z, const acb_t x, const arb_t y, slong prec)
void acb_sub(acb_t z, const acb_t x, const acb_t y, slong prec)
void acb_sub_ui(acb_t z, const acb_t x, ulong c, slong prec)
void acb_sub_si(acb_t z, const acb_t x, slong c, slong prec)
void acb_sub_fmpz(acb_t z, const acb_t x, const fmpz_t y, slong prec)
void acb_sub_arb(acb_t z, const acb_t x, const arb_t y, slong prec)

Set `z` to $x + y$, respectively $x - y$, by real additions of the parts.
:::

::: function
void acb_mul_onei(acb_t z, const acb_t x)
void acb_div_onei(acb_t z, const acb_t x)
void acb_mul_i_pow_si(acb_t z, const acb_t x, slong k)

Set `z` to $ix$, $x/i = -ix$, respectively $i^k x$. Exact.
:::

::: function
void acb_mul_2exp_si(acb_t z, const acb_t x, slong e)
void acb_mul_2exp_fmpz(acb_t z, const acb_t x, const fmpz_t c)

Set `z` to $x \cdot 2^e$. Exact.
:::

::: function
void acb_mul(acb_t z, const acb_t x, const acb_t y, slong prec)

Sets `z` to $xy$. If a part of `x` or `y` is exactly zero, two real
multiplications are used; if `x` and `y` are the same object, the squaring
formula $(a + bi)^2 = (a^2 - b^2) + 2abi$; otherwise the midpoints are
multiplied with `arf_complex_mul` and the radii bounded componentwise.
:::

::: function
void acb_mul_naive(acb_t z, const acb_t x, const acb_t y, slong prec)

Sets `z` to $xy$ using four real multiplications with `arb_mul`, as
$(ac - bd) + (ad + bc) i$. The result is usually slightly wider than that
of `acb_mul`; this version is kept for testing.
:::

::: function
void acb_mul_ui(acb_t z, const acb_t x, ulong y, slong prec)
void acb_mul_si(acb_t z, const acb_t x, slong y, slong prec)
void acb_mul_fmpz(acb_t z, const acb_t x, const fmpz_t y, slong prec)
void acb_mul_arf(acb_t z, const acb_t x, const arf_t y, slong prec)
void acb_mul_arb(acb_t z, const acb_t x, const arb_t y, slong prec)

Set `z` to $xy$ for a real `y`, multiplying both parts.
:::

::: function
void acb_sqr(acb_t res, const acb_t val, slong prec)
void acb_cube(acb_t y, const acb_t x, slong prec)

Set the output to $x^2$, respectively $x^3$. The cube is computed with two
real squarings, two real multiplications and some scalar operations.
:::

::: function
void acb_addmul(acb_t z, const acb_t x, const acb_t y, slong prec)
void acb_addmul_ui(acb_t z, const acb_t x, ulong y, slong prec)
void acb_addmul_si(acb_t z, const acb_t x, slong y, slong prec)
void acb_addmul_fmpz(acb_t z, const acb_t x, const fmpz_t y, slong prec)
void acb_addmul_arb(acb_t z, const acb_t x, const arb_t y, slong prec)
void acb_submul(acb_t z, const acb_t x, const acb_t y, slong prec)
void acb_submul_ui(acb_t z, const acb_t x, ulong y, slong prec)
void acb_submul_si(acb_t z, const acb_t x, slong y, slong prec)
void acb_submul_fmpz(acb_t z, const acb_t x, const fmpz_t y, slong prec)
void acb_submul_arb(acb_t z, const acb_t x, const arb_t y, slong prec)

Set `z` to $z + xy$, respectively $z - xy$.
:::

::: function
void acb_inv(acb_t z, const acb_t x, slong prec)

Sets `z` to $1/x$. If `x` contains zero, the result is indeterminate.

::: details
For $x = a + bi$ with both parts nonzero, the midpoint is computed as
$(a - bi)/(a^2 + b^2)$ with $a^2 + b^2$ rounded once (`arf_sosq`), and the
propagated error of each part is bounded directly from the exact
expression of the perturbation, for example
$$
\frac{a + x'}{(a + x')^2 + (b + y')^2} - \frac{a}{a^2 + b^2} =
\frac{(b^2 - a^2) x' - a(x'^2 + y'^2 + 2 b y')}{(a^2 + b^2)((a + x')^2
+ (b + y')^2)}
$$
for the real part, with $|x'| \le r_1$, $|y'| \le r_2$. If a part of `x`
is exactly zero, a real inversion is used.
:::
:::

::: function
void acb_div(acb_t z, const acb_t x, const acb_t y, slong prec)
void acb_div_ui(acb_t z, const acb_t x, ulong c, slong prec)
void acb_div_si(acb_t z, const acb_t x, slong c, slong prec)
void acb_div_fmpz(acb_t z, const acb_t x, const fmpz_t c, slong prec)
void acb_div_arb(acb_t z, const acb_t x, const arb_t c, slong prec)

Set `z` to $x / y$. If the divisor contains zero, the result is
indeterminate. If the divisor is real or purely imaginary, the parts of
`x` are divided directly; otherwise $x \cdot (1/y)$ is computed with
`acb_inv`.
:::

## Dot products {#dot}

::: function
void acb_dot(acb_t res, const acb_t initial, int subtract, acb_srcptr x, slong xstep, acb_srcptr y, slong ystep, slong len, slong prec)
void acb_dot_simple(acb_t res, const acb_t initial, int subtract, acb_srcptr x, slong xstep, acb_srcptr y, slong ystep, slong len, slong prec)
void acb_dot_precise(acb_t res, const acb_t initial, int subtract, acb_srcptr x, slong xstep, acb_srcptr y, slong ystep, slong len, slong prec)

Set `res` to $s + (-1)^{\mathit{subtract}} \sum_{i=0}^{\mathit{len}-1}
x_i y_i$, where $s$ is `initial` (0 if `NULL`), reading the entries with
strides `xstep` and `ystep`.

::: details
The conventions are those of `arb_dot`: `subtract` is 0 or 1, a negative
`len` counts as 0, strides may be negative, and `res` may alias `initial`
but not the entries of `x` and `y`. `acb_dot` accumulates the real and
imaginary parts in two fixed-point sums with a single rounding each; it
is the function to use. `acb_dot_simple` is a loop of `acb_addmul`, used
as a fallback for huge exponents, and `acb_dot_precise` rounds the exact
result once; both are mainly for testing. See
[](algorithms/ball-arithmetic/dot-products).
:::
:::

::: function
void acb_approx_dot(acb_t res, const acb_t initial, int subtract, acb_srcptr x, slong xstep, acb_srcptr y, slong ystep, slong len, slong prec)

Computes the same dot product without error bounds: the radii of the
inputs are ignored and only the midpoint of `res` is written.
:::

::: function
void acb_dot_ui(acb_t res, const acb_t initial, int subtract, acb_srcptr x, slong xstep, const ulong * y, slong ystep, slong len, slong prec)
void acb_dot_si(acb_t res, const acb_t initial, int subtract, acb_srcptr x, slong xstep, const slong * y, slong ystep, slong len, slong prec)
void acb_dot_uiui(acb_t res, const acb_t initial, int subtract, acb_srcptr x, slong xstep, const ulong * y, slong ystep, slong len, slong prec)
void acb_dot_siui(acb_t res, const acb_t initial, int subtract, acb_srcptr x, slong xstep, const ulong * y, slong ystep, slong len, slong prec)
void acb_dot_fmpz(acb_t res, const acb_t initial, int subtract, acb_srcptr x, slong xstep, const fmpz * y, slong ystep, slong len, slong prec)

Like `acb_dot` with a real integer vector `y`; the `uiui` and `siui`
versions take two-limb integers (unsigned, respectively signed in two's
complement) stored as pairs of `ulong`, low limb first.
:::

## Powers and roots {#powers}

::: function
void acb_sqrt(acb_t y, const acb_t x, slong prec)
void acb_sqrt_analytic(acb_t res, const acb_t z, int analytic, slong prec)

Set the output to the principal square root, with the branch cut on
$(-\infty, 0)$ and $\sqrt{-1} = i$ (continuous from above). If `analytic`
is set, the result is indeterminate when the input touches the branch
cut, including an exact negative real number.

::: details
If one part of $z = a + bi$ is exactly zero, a single real square root
suffices. In general $\sqrt{a + bi} = u/2 + ib/u$ with $u = \sqrt{2(|z| +
a)}$, which needs two real square roots; for $a < 0$ the roles of the
parts are exchanged to avoid cancellation.
:::
:::

::: function
void acb_rsqrt(acb_t y, const acb_t x, slong prec)
void acb_rsqrt_analytic(acb_t res, const acb_t z, int analytic, slong prec)

Set the output to the principal reciprocal square root $1/\sqrt{z}$.
`analytic` is as for `acb_sqrt_analytic`.

::: details
If one part is zero, a single real reciprocal square root suffices. In
general $1/\sqrt{a + bi} = ((a + r) - bi)/v$ with $r = |a + bi|$ and $v =
\sqrt{r |a + bi + r|^2}$, which needs one real square root and one real
reciprocal square root.
:::
:::

::: function
void acb_sqrts(acb_t y1, acb_t y2, const acb_t x, slong prec)

Sets `y1` and `y2` to the two square roots of `x`, in an unspecified
order. Unlike `acb_sqrt`, no accuracy is lost when `x` is near the branch
cut.
:::

::: function
void acb_root_ui(acb_t y, const acb_t x, ulong k, slong prec)

Sets `y` to the principal `k`-th root $\exp(\log(x)/k)$, using square
roots for $k = 2$ and $k = 4$ and `arb_root_ui` for nonnegative real
input. For $k = 0$ the result is indeterminate.
:::

::: function
void acb_quadratic_roots_fmpz(acb_t r1, acb_t r2, const fmpz_t a, const fmpz_t b, const fmpz_t c, slong prec)

Sets `r1` and `r2` to the roots of $ax^2 + bx + c$, $a \ne 0$. Both roots
are computed accurately, also when the quadratic formula would suffer
from cancellation.
:::

::: function
void acb_pow_ui(acb_t y, const acb_t b, ulong e, slong prec)
void acb_pow_si(acb_t y, const acb_t b, slong e, slong prec)
void acb_pow_fmpz(acb_t y, const acb_t b, const fmpz_t e, slong prec)

Set `y` to $b^e$ by binary exponentiation (with an inversion if $e <
0$). For huge exponents `acb_pow` can be faster.
:::

::: function
void acb_pow(acb_t r, const acb_t x, const acb_t y, slong prec)
void acb_pow_arb(acb_t z, const acb_t x, const arb_t y, slong prec)
void acb_pow_analytic(acb_t res, const acb_t z, const acb_t w, int analytic, slong prec)

Set the output to the principal power $x^y$: by binary exponentiation if
$y$ is a small exact integer, as $(x^{1/2})^{2y}$ if $y$ is a small exact
half-integer, and as $\exp(y \log x)$ otherwise. With `analytic` set, the
result is indeterminate if $x$ touches the branch cut of the logarithm,
unless $y$ is an integer.
:::

::: function
void acb_unit_root(acb_t res, ulong order, slong prec)

Sets `res` to the root of unity $\exp(2\pi i / \mathit{order})$.
:::

## Exponentials and logarithms {#exp-log}

::: function
void acb_exp(acb_t r, const acb_t z, slong prec)

Sets `r` to $e^z = e^a (\cos b + i \sin b)$ for $z = a + bi$.
:::

::: function
void acb_exp_pi_i(acb_t r, const acb_t z, slong prec)

Sets `r` to $\exp(\pi i z)$, using `arb_sin_cos_pi` so that the result is
exact or accurate at rational points.
:::

::: function
void acb_exp_invexp(acb_t r, acb_t s, const acb_t z, slong prec)

Sets `r` to $e^z$ and `s` to $e^{-z}$.
:::

::: function
void acb_expm1(acb_t r, const acb_t z, slong prec)

Sets `r` to $e^z - 1$, accurately also for $z$ near 0.
:::

::: function
void acb_log(acb_t r, const acb_t z, slong prec)
void acb_log_analytic(acb_t res, const acb_t z, int analytic, slong prec)

Set the output to the principal logarithm $\log z = \log|z| + i \arg z$,
computed with `arb_log_hypot` and `arb_atan2`, with the branch cut on
$(-\infty, 0]$ and imaginary part in $(-\pi, \pi]$. A rectangle crossing
the branch cut gives an imaginary part enclosing $[-\pi, \pi]$; with
`analytic` set the result is then indeterminate instead.
:::

::: function
void acb_log1p(acb_t r, const acb_t z, slong prec)

Sets `r` to $\log(1 + z)$, accurately also for $z$ near 0.
:::

## Trigonometric functions {#trig}

::: function
void acb_sin(acb_t r, const acb_t z, slong prec)
void acb_cos(acb_t r, const acb_t z, slong prec)
void acb_sin_cos(acb_t s, acb_t c, const acb_t z, slong prec)

Set the output to $\sin z$ and $\cos z$, computed as $\sin(a + bi) = \sin
a \cosh b + i \cos a \sinh b$ and $\cos(a + bi) = \cos a \cosh b - i \sin a
\sinh b$.
:::

::: function
void acb_tan(acb_t r, const acb_t z, slong prec)
void acb_cot(acb_t r, const acb_t z, slong prec)

Set `r` to $\tan z$, respectively $\cot z$. For large $|\operatorname{Im}
z|$ they are evaluated stably as $\pm i$ plus an exponentially small
correction.
:::

::: function
void acb_sin_pi(acb_t r, const acb_t z, slong prec)
void acb_cos_pi(acb_t r, const acb_t z, slong prec)
void acb_sin_cos_pi(acb_t s, acb_t c, const acb_t z, slong prec)
void acb_tan_pi(acb_t r, const acb_t z, slong prec)
void acb_cot_pi(acb_t r, const acb_t z, slong prec)
void acb_csc_pi(acb_t y, const acb_t x, slong prec)

Set the output to $\sin(\pi z)$, $\cos(\pi z)$, both, $\tan(\pi z)$,
$\cot(\pi z)$, respectively $\csc(\pi z)$, using `arb_sin_cos_pi` for the
real part so that the results are accurate near the zeros.
:::

::: function
void acb_sec(acb_t y, const acb_t x, slong prec)
void acb_csc(acb_t y, const acb_t x, slong prec)

Set `y` to $\sec x = 1/\cos x$, respectively $\csc x = 1/\sin x$.
:::

::: function
void acb_sinc(acb_t res, const acb_t z, slong prec)
void acb_sinc_pi(acb_t res, const acb_t z, slong prec)

Set `res` to $\sin(z)/z$, respectively $\sin(\pi z)/(\pi z)$, with the
value 1 at $z = 0$.
:::

## Inverse trigonometric functions {#inverse-trig}

::: function
void acb_asin(acb_t r, const acb_t z, slong prec)
void acb_acos(acb_t r, const acb_t z, slong prec)
void acb_atan(acb_t r, const acb_t z, slong prec)

Set `r` to $\operatorname{asin} z = -i \log(iz + \sqrt{1 - z^2})$,
$\operatorname{acos} z = \pi/2 - \operatorname{asin} z$, respectively
$\operatorname{atan} z = \tfrac{i}{2} (\log(1 - iz) - \log(1 + iz))$. The
branch cuts are $(-\infty, -1) \cup (1, \infty)$ for $\operatorname{asin}$
and $\operatorname{acos}$, and $(-i\infty, -i) \cup (i, i\infty)$ for
$\operatorname{atan}$.
:::

## Hyperbolic functions {#hyperbolic}

::: function
void acb_sinh(acb_t y, const acb_t x, slong prec)
void acb_cosh(acb_t y, const acb_t x, slong prec)
void acb_sinh_cosh(acb_t y, acb_t z, const acb_t x, slong prec)
void acb_tanh(acb_t y, const acb_t x, slong prec)
void acb_coth(acb_t y, const acb_t x, slong prec)

Compute $\sinh z = -i \sin(iz)$, $\cosh z = \cos(iz)$, both, $\tanh z =
-i \tan(iz)$, respectively $\coth z = i \cot(iz)$.
:::

::: function
void acb_sech(acb_t r, const acb_t z, slong prec)
void acb_csch(acb_t r, const acb_t z, slong prec)

Set `r` to $\operatorname{sech} z = 1/\cosh z$, respectively
$\operatorname{csch} z = 1/\sinh z$.
:::

## Inverse hyperbolic functions {#inverse-hyperbolic}

::: function
void acb_asinh(acb_t r, const acb_t z, slong prec)
void acb_acosh(acb_t r, const acb_t z, slong prec)
void acb_atanh(acb_t r, const acb_t z, slong prec)

Set `r` to $\operatorname{asinh} z = -i \operatorname{asin}(iz)$,
$\operatorname{acosh} z = \log(z + \sqrt{z + 1} \sqrt{z - 1})$,
respectively $\operatorname{atanh} z = -i \operatorname{atan}(iz)$.
:::

## Constants {#constants}

::: function
void acb_const_pi(acb_t x, slong prec)

Sets `x` to $\pi$ (see `arb_const_pi`).
:::

## Lambert W function {#lambertw}

::: macro
ACB_LAMBERTW_LEFT
ACB_LAMBERTW_MIDDLE

Values of the `flags` argument of `acb_lambertw` that select the
nonstandard branches $W_{\mathrm{left}|k}$ and $W_{\mathrm{middle}}$.
:::

::: function
void acb_lambertw(acb_t res, const acb_t z, const fmpz_t k, int flags, slong prec)

Sets `res` to the branch $W_k(z)$ of the Lambert W function, the
solution of $w e^w = z$; $k = 0$ is the principal branch. The branch cuts
are those of [@CGHJK1996]. The algorithm is described in [@Joh2017b]
and on [](algorithms/special-functions/lambert-w).

::: details
With `flags` equal to `ACB_LAMBERTW_LEFT`, the function computes
$W_{\mathrm{left}|k}(z)$, which is $W_k(z)$ in the upper half plane and
$W_{k+1}(z)$ in the lower half plane, connected continuously to the left
of the branch points; the branch cut on $(-\infty, 0)$ is rotated
counterclockwise to $(0, +\infty)$. (For $k = -1$ and $k = 0$ there is
also a branch cut on $(-1/e, 0)$, continuous from below instead of from
above to keep counterclockwise continuity.)

With `ACB_LAMBERTW_MIDDLE`, it computes $W_{\mathrm{middle}}(z)$, which is
$W_{-1}(z)$ in the upper half plane and $W_1(z)$ in the lower half plane,
connected continuously through $(-1/e, 0)$, with branch cuts on
$(-\infty, -1/e)$ and $(0, +\infty)$. It extends the real analytic
function $W_{-1}(x)$ on $(-1/e, 0)$ to a complex analytic function,
whereas the standard $W_{-1}(z)$ has a branch cut on that segment.
:::
:::

::: function
void acb_lambertw_asymp(acb_t res, const acb_t z, const fmpz_t k, slong L, slong M, slong prec)

Sets `res` to $W_k(z)$ computed from $L$ and $M$ terms of the bivariate
asymptotic series at 0 or $\infty$. The formula is valid everywhere, but
the error bound is finite only when $|\log z|$ is large enough.
:::

::: function
int acb_lambertw_check_branch(const acb_t w, const fmpz_t k, slong prec)

Returns 1 if `w` certainly lies in the image of the branch $W_k$. This is
used to verify that a computed approximation is on the intended branch;
it is false for points on or overlapping the branch cuts.
:::

::: function
void acb_lambertw_bound_deriv(mag_t res, const acb_t z, const acb_t ez1, const fmpz_t k)

Sets `res` to an upper bound for $|W_k'(z)|$, given `ez1` $= ez + 1$.
On the real line the directional derivative is taken, so the caller must
handle the discontinuity at the branch cut separately.
:::

## Gamma function and rising factorials {#gamma}

The algorithms are described in [](algorithms/special-functions/gamma).

::: function
void acb_gamma(acb_t y, const acb_t x, slong prec)
void acb_rgamma(acb_t y, const acb_t x, slong prec)

Set `y` to $\Gamma(x)$, respectively $1/\Gamma(x)$ (which is 0 at the
poles of $\Gamma$). Aliases of `acb_hypgeom_gamma` and
`acb_hypgeom_rgamma`.
:::

::: function
void acb_lgamma(acb_t y, const acb_t x, slong prec)

Sets `y` to $\log \Gamma(x)$, an alias of `acb_hypgeom_lgamma`.

::: details
The branch cut is on the negative real axis, so that $\log \Gamma(z) +
\log z = \log \Gamma(z + 1)$ holds for all $z$, whereas $\log \Gamma(z)
\ne \log(\Gamma(z))$ in general. In the left half plane, the reflection
formula is evaluated with the correct branch structure via
`acb_log_sin_pi`.
:::
:::

::: function
void acb_log_sin_pi(acb_t res, const acb_t z, slong prec)

Sets `res` to the logarithmic sine function $S(z) = \log \pi - \log
\Gamma(z) + \log \Gamma(1 - z)$, the branch of $\log \sin(\pi z)$ used by
the reflection formula of $\log \Gamma$.

::: details
Equivalently, $S(z) = \int_{1/2}^z \pi \cot(\pi t)\, dt$ with the path
through the upper half plane if $0 < \arg z \le \pi$ and through the lower
half plane otherwise, or

$$
S(z) = \log(\sin(\pi(z - n))) \mp n \pi i, \quad n = \lfloor
\operatorname{Re} z \rfloor,
$$

with the negative sign if $0 < \arg z \le \pi$ (the union of both cases
is computed when $\arg z$ is not determined). After subtracting $n$,
$0 \le \operatorname{Re} z < 1$, and $S(z) = \log \sin(\pi z)$ is used if
the imaginary part is small; otherwise $S(z) = i \pi (z - 1/2) + \log((1 +
e^{-2 \pi i z})/2)$ in the lower half plane and the conjugate expression
in the upper half plane, which avoid overflow. The function is evaluated
at the midpoint and the propagated error is bounded with $S'(z)$, so that
the result changes continuously when $z$ is not real and $n$ is not
determined.
:::
:::

::: function
void acb_digamma(acb_t y, const acb_t x, slong prec)

Sets `y` to the digamma function $\psi(x) = \Gamma'(x)/\Gamma(x)$.
:::

::: function
void acb_polygamma(acb_t res, const acb_t s, const acb_t z, slong prec)

Sets `res` to the generalised polygamma function $\psi(s, z)$.

::: details
For a nonnegative integer $s$ this is the $s$-th derivative of the
digamma function: $s = 0$ calls `acb_digamma` and integers $s \ge 1$ the
Hurwitz zeta function (for small $s$, reading the coefficients of
`acb_poly_digamma_series` can be faster). Other $s$ use the
generalisation of Espinosa and Moll [@EM2004],
$$
\psi(s, z) = \frac{\zeta'(s+1, z) + (\gamma + \psi(-s)) \zeta(s+1,
z)}{\Gamma(-s)}.
$$
:::
:::

::: function
void acb_barnes_g(acb_t res, const acb_t z, slong prec)
void acb_log_barnes_g(acb_t res, const acb_t z, slong prec)

Set `res` to the Barnes G-function, respectively its logarithm, which has
its branch cut on the negative real axis and satisfies $\log G(z + 1) =
\log \Gamma(z) + \log G(z)$ for all $z$.

::: details
Small integers use the recurrence $G(z + 1) = \Gamma(z) G(z)$ with $G(1)
= 1$. Otherwise $\log G(z) = (z - 1) \log \Gamma(z) - \zeta'(-1, z) +
\zeta'(-1)$.
:::
:::

::: function
void acb_rising_ui(acb_t z, const acb_t x, ulong n, slong prec)
void acb_rising(acb_t z, const acb_t x, const acb_t n, slong prec)

Set `z` to the rising factorial $x (x + 1) \cdots (x + n - 1)$. Aliases of
`acb_hypgeom_rising_ui` and `acb_hypgeom_rising`.
:::

::: function
void acb_rising2_ui(acb_t u, acb_t v, const acb_t x, ulong n, slong prec)

Sets `u` to the rising factorial $u(x)$ and `v` to $u'(x)$. Wrapper of
`acb_hypgeom_rising_ui_jet`.
:::

::: function
void acb_rising_ui_get_mag(mag_t bound, const acb_t s, ulong n)

Sets `bound` to an upper bound for $|s (s + 1) \cdots (s + n - 1)|$. Not
optimised for large $n$.
:::

## Zeta function, Bernoulli polynomials and polylogarithms {#zeta}

::: function
void acb_zeta(acb_t z, const acb_t s, slong prec)
void acb_hurwitz_zeta(acb_t z, const acb_t s, const acb_t a, slong prec)

Set `z` to $\zeta(s)$, respectively $\zeta(s, a)$. Wrappers of
`acb_dirichlet_zeta` and `acb_dirichlet_hurwitz`; for derivatives with
respect to $s$, use `acb_poly_zeta_series`.
:::

::: function
void acb_bernoulli_poly_ui(acb_t res, ulong n, const acb_t x, slong prec)

Sets `res` to the Bernoulli polynomial $B_n(x)$. This is fast only if $n$
or $x$ is a small integer; cached Bernoulli numbers are used but the cache
is not extended.
:::

::: function
void acb_polylog(acb_t w, const acb_t s, const acb_t z, slong prec)
void acb_polylog_si(acb_t w, slong s, const acb_t z, slong prec)

Set `w` to the polylogarithm $\operatorname{Li}_s(z)$, with the branch cut
on $(1, +\infty)$.
:::

## Arithmetic-geometric mean {#agm}

The algorithms are described in
[](algorithms/special-functions/elliptic-modular).

::: function
void acb_agm1(acb_t m, const acb_t z, slong prec)

Sets `m` to $M(z) = \operatorname{agm}(1, z)$, defined by always choosing
the optimal branch of the square root in the AGM iteration. $M$ is
continuous except for a branch cut on the negative real axis, where it is
continuous from above.
:::

::: function
void acb_agm1_cpx(acb_ptr m, const acb_t z, slong len, slong prec)

Sets `m` to the first `len` coefficients of the power series $M(z + x)
\in \mathbb{C}[[x]]$.
:::

::: function
void acb_agm(acb_t res, const acb_t a, const acb_t b, slong prec)

Sets `res` to the arithmetic-geometric mean of `a` and `b`, with the
square roots chosen to give the optimal AGM sequence. This is well
defined except when $a/b$ is a negative real number, where there are two
optimal sequences; then an arbitrary but consistent choice is made, or the
union of both when the choice cannot be decided.
:::

## Chebyshev polynomials {#chebyshev}

::: function
void acb_chebyshev_t_ui(acb_t a, ulong n, const acb_t x, slong prec)
void acb_chebyshev_u_ui(acb_t a, ulong n, const acb_t x, slong prec)

Set `a` to $T_n(x)$, respectively $U_n(x)$.
:::

::: function
void acb_chebyshev_t2_ui(acb_t a, acb_t b, ulong n, const acb_t x, slong prec)
void acb_chebyshev_u2_ui(acb_t a, acb_t b, ulong n, const acb_t x, slong prec)

Set $a = T_n(x)$, $b = T_{n-1}(x)$, respectively $a = U_n(x)$, $b =
U_{n-1}(x)$. The outputs and the input must be distinct objects.
:::

## Piecewise real functions {#piecewise}

These functions extend common piecewise real functions to piecewise
analytic functions of a complex variable, for use with rigorous
integration ([acb_calc](acb_calc)), which requires to know where the
integrand is analytic. If `analytic` is set, evaluation at a
discontinuity or a non-analytic point gives an indeterminate result.

::: function
void acb_real_abs(acb_t res, const acb_t z, int analytic, slong prec)

The absolute value, extended to $z$ in the right half plane and $-z$ in
the left half plane, with a discontinuity on $\operatorname{Re} z = 0$.
:::

::: function
void acb_real_sgn(acb_t res, const acb_t z, int analytic, slong prec)

The sign function, extended to 1 in the right half plane and $-1$ in the
left half plane, with a discontinuity on $\operatorname{Re} z = 0$.
Without `analytic`, it is the same as `acb_csgn`.
:::

::: function
void acb_real_heaviside(acb_t res, const acb_t z, int analytic, slong prec)

The Heaviside step function, extended to 1 in the right half plane and 0
in the left half plane, with a discontinuity on $\operatorname{Re} z =
0$.
:::

::: function
void acb_real_floor(acb_t res, const acb_t z, int analytic, slong prec)
void acb_real_ceil(acb_t res, const acb_t z, int analytic, slong prec)

The floor and ceiling functions, extended to the constants $n$,
respectively $n + 1$, on the strips $n < \operatorname{Re} z < n + 1$,
with discontinuities on the lines $\operatorname{Re} z = n$.
:::

::: function
void acb_real_max(acb_t res, const acb_t x, const acb_t y, int analytic, slong prec)
void acb_real_min(acb_t res, const acb_t x, const acb_t y, int analytic, slong prec)

$\max(x, y)$ extended to $x$ where $\operatorname{Re} x \ge
\operatorname{Re} y$ and $y$ elsewhere, and $\min(x, y)$ extended to $x$
where $\operatorname{Re} x \le \operatorname{Re} y$ and $y$ elsewhere,
with discontinuities where $\operatorname{Re} x = \operatorname{Re} y$.
:::

::: function
void acb_real_sqrtpos(acb_t res, const acb_t z, int analytic, slong prec)

The square root on $[0, +\infty)$, extended to the principal complex
square root. For a real `z` without `analytic`, only the nonnegative part
is used, as by `arb_sqrtpos`. This is the function to integrate
$\sqrt{f(x)}$ when $f(x) \ge 0$ is known: unlike `acb_sqrt_analytic`, it
creates no spurious imaginary parts $[\pm \varepsilon] i$ when an
enclosure of $f(x)$ straddles 0.
:::

## Random generation {#random}

::: function
void acb_randtest(acb_t z, flint_rand_t state, slong prec, slong mag_bits)
void acb_randtest_special(acb_t z, flint_rand_t state, slong prec, slong mag_bits)
void acb_randtest_precise(acb_t z, flint_rand_t state, slong prec, slong mag_bits)
void acb_randtest_param(acb_t z, flint_rand_t state, slong prec, slong mag_bits)

Set `z` to a random complex number for testing, with real and imaginary
parts generated like those of `arb_randtest`: in general, with
occasional NaN and infinities, with precise parts, and with a high
probability of integers and half-integers (useful for parameters of
special functions).
:::

::: function
void acb_urandom(acb_t z, flint_rand_t state, slong prec)

Sets `z` to a random number with precise parts, uniformly distributed in
the unit disk.
:::

## Vector functions {#vectors}

::: function
acb_ptr _acb_vec_entry_ptr(acb_ptr vec, slong i)

Returns a pointer to entry `i` of `vec`.
:::

::: function
void _acb_vec_zero(acb_ptr A, slong n)
void _acb_vec_indeterminate(acb_ptr vec, slong len)

Set all entries to zero, respectively to indeterminate values.
:::

::: function
int _acb_vec_is_zero(acb_srcptr vec, slong len)
int _acb_vec_is_real(acb_srcptr v, slong len)
int _acb_vec_is_finite(acb_srcptr vec, slong len)

Return whether all entries are exactly zero, have exactly zero imaginary
part, respectively are finite.
:::

::: function
int _acb_vec_equal(acb_srcptr vec1, acb_srcptr vec2, slong len)
int _acb_vec_overlaps(acb_srcptr vec1, acb_srcptr vec2, slong len)
int _acb_vec_contains(acb_srcptr vec1, acb_srcptr vec2, slong len)

Apply `acb_equal`, `acb_overlaps`, respectively `acb_contains` entrywise
and return whether it holds for all entries.
:::

::: function
void _acb_vec_set(acb_ptr res, acb_srcptr vec, slong len)
void _acb_vec_set_round(acb_ptr res, acb_srcptr vec, slong len, slong prec)
void _acb_vec_swap(acb_ptr res, acb_ptr vec, slong len)
void _acb_vec_neg(acb_ptr res, acb_srcptr vec, slong len)
void _acb_vec_trim(acb_ptr res, acb_srcptr vec, slong len)

Copy, copy with rounding, swap, negate, and apply `acb_trim` entrywise.
:::

::: function
void _acb_vec_get_real(arb_ptr re, acb_srcptr vec, slong len)
void _acb_vec_get_imag(arb_ptr im, acb_srcptr vec, slong len)
void _acb_vec_set_real_imag(acb_ptr vec, arb_srcptr re, arb_srcptr im, slong len)

Extract the real or imaginary parts of a vector, and build a vector from
its real and imaginary parts.
:::

::: function
void _acb_vec_add(acb_ptr res, acb_srcptr vec1, acb_srcptr vec2, slong len, slong prec)
void _acb_vec_sub(acb_ptr res, acb_srcptr vec1, acb_srcptr vec2, slong len, slong prec)
void _acb_vec_sqr(acb_ptr res, acb_srcptr vec, slong len, slong prec)

Entrywise sum, difference and square.
:::

::: function
void _acb_vec_scalar_mul(acb_ptr res, acb_srcptr vec, slong len, const acb_t c, slong prec)
void _acb_vec_scalar_mul_arb(acb_ptr res, acb_srcptr vec, slong len, const arb_t c, slong prec)
void _acb_vec_scalar_mul_arf(acb_ptr res, acb_srcptr vec, slong len, const arf_t c, slong prec)
void _acb_vec_scalar_mul_ui(acb_ptr res, acb_srcptr vec, slong len, ulong c, slong prec)
void _acb_vec_scalar_mul_si(acb_ptr res, acb_srcptr vec, slong len, slong c, slong prec)
void _acb_vec_scalar_mul_fmpz(acb_ptr res, acb_srcptr vec, slong len, const fmpz_t c, slong prec)
void _acb_vec_scalar_mul_2exp_si(acb_ptr res, acb_srcptr vec, slong len, slong c)
void _acb_vec_scalar_mul_onei(acb_ptr res, acb_srcptr vec, slong len)
void _acb_vec_scalar_div(acb_ptr res, acb_srcptr vec, slong len, const acb_t c, slong prec)
void _acb_vec_scalar_div_arb(acb_ptr res, acb_srcptr vec, slong len, const arb_t c, slong prec)
void _acb_vec_scalar_div_ui(acb_ptr res, acb_srcptr vec, slong len, ulong c, slong prec)
void _acb_vec_scalar_div_fmpz(acb_ptr res, acb_srcptr vec, slong len, const fmpz_t c, slong prec)
void _acb_vec_scalar_addmul(acb_ptr res, acb_srcptr vec, slong len, const acb_t c, slong prec)
void _acb_vec_scalar_submul(acb_ptr res, acb_srcptr vec, slong len, const acb_t c, slong prec)

Multiply each entry of `vec` by $c$ (by $2^c$ for `_2exp_si`, by $i$ for
`_onei`) or divide it by $c$, writing to `res`; the `addmul` and `submul`
versions add or subtract $c \cdot \mathit{vec}_i$ to $\mathit{res}_i$.
:::

::: function
slong _acb_vec_bits(acb_srcptr vec, slong len)

Returns the maximum of `acb_bits` over the entries.
:::

::: function
void _acb_vec_set_powers(acb_ptr xs, const acb_t x, slong len, slong prec)

Sets `xs` to $1, x, x^2, \ldots, x^{\mathit{len}-1}$.
:::

::: function
void _acb_vec_unit_roots(acb_ptr z, slong order, slong len, slong prec)

Sets `z` to $1, \omega, \omega^2, \ldots, \omega^{\mathit{len}-1}$ with
$\omega = \exp(2\pi i/\mathit{order})$; `order` may be negative. To avoid
error growth, the powers are not computed by repeated multiplication.
:::

::: function
void _acb_vec_add_error_arf_vec(acb_ptr res, arf_srcptr err, slong len)
void _acb_vec_add_error_mag_vec(acb_ptr res, mag_srcptr err, slong len)

Add $|\mathit{err}_i|$ to both radii of $\mathit{res}_i$.
:::

::: function
int _acb_vec_get_unique_fmpz_vec(fmpz * res, acb_srcptr vec, slong len)

Applies `acb_get_unique_fmpz` entrywise and returns nonzero if every
entry contains a unique integer.
:::

::: function
void _acb_vec_sort_pretty(acb_ptr vec, slong len)

Sorts the vector by real and imaginary parts, to make printed output
easier to read. This is not a rigorous ordering.
:::

::: function
void _acb_vec_printd(acb_srcptr vec, slong len, slong ndigits)
void _acb_vec_printn(acb_srcptr vec, slong len, slong ndigits, ulong flags)

Print the entries with `acb_printd`, respectively `acb_printn`.
:::
