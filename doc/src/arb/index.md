---
title: arb.h
header: arb.h
description: Real numbers as midpoint-radius balls with rigorous error propagation
category: Real and complex numbers
order: 3
---

An `arb_t` is a *ball* $[m \pm r] = [m - r, m + r]$ with an `arf_t`
midpoint $m$ and a `mag_t` radius $r \ge 0$. It represents a real number
known only up to the error $r$. Every function in this module satisfies
the *inclusion principle*: if each input ball contains some real number,
the output ball contains the exact result of the operation applied to
those numbers. Rounding errors and propagated input errors are added to
the radius automatically, so a computation that only calls `arb`
functions produces a rigorous enclosure without any error analysis by the
user. The output ball is usually not the smallest possible, but for the
basic operations it is close to optimal.

The module provides arithmetic, comparisons, conversions, elementary
functions, mathematical constants and a number of special functions of
real arguments (gamma, zeta, Bernoulli numbers, polylogarithms, Lambert
W, partitions). Complex balls are in `acb`, polynomials and power series
in `arb_poly` (see [](arb_poly/power-series) for Taylor expansions of
functions), matrices in `arb_mat`, and most other special functions
in `arb_hypgeom`. New users should start with [](arb/getting-started)
and the user guide [](guide/types/ball-arithmetic). The pages
[](arb/arithmetic), [](arb/predicates), [](arb/printing) and
[](arb/elementary) explain precision, comparisons, decimal conversion and
the elementary functions in detail. The tutorials
[](tutorials/numerics/constants), [](tutorials/numerics/series) and
[](tutorials/numerics/adaptive-precision) build complete programs with
real balls.

::: summary algorithms/ball-arithmetic/overview :::

## Types and conventions {#types}

::: type
arb_struct
arb_t
arb_ptr
arb_srcptr

An `arb_struct` has two fields, an `arf_struct` midpoint and a
`mag_struct` radius. `arb_t` is an array of one `arb_struct`, so that it is
passed by reference; `arb_ptr` and `arb_srcptr` point to arrays (vectors)
of balls.
:::

::: macro
arb_midref(x)
arb_radref(x)

Pointers to the midpoint (an `arf_t`) and the radius (a `mag_t`) of `x`,
for use with the functions of `arf` and `mag`.
:::

::: function
arf_ptr arb_mid_ptr(arb_t z)
mag_ptr arb_rad_ptr(arb_t z)

Function versions of `arb_midref` and `arb_radref`, for language bindings
that cannot use macros.
:::

::: macro
ARB_RND

The rounding mode used for midpoints, `ARF_RND_DOWN`. Rounding towards
zero is the cheapest mode, and the rounding error is added to the radius
anyway.
:::

*Precision.* A function with a `prec` argument rounds the midpoint of its
result to at most `prec` bits and adds the rounding error to the radius.
For arithmetic operations this is exactly the precision of the
corresponding `arf` operation, and an exact result (for example the sum of
two small integers) is returned exactly. For other functions `prec` is the
working precision: the implementation aims at a relative error of a few
units in the last place for well-conditioned input and may use more
precision internally, but the radius of the result is only guaranteed to
be a bound, not to be small. If the input is exact or computed with
increasing accuracy, the output converges as `prec` grows. [](arb/arithmetic) and the user guide [](guide/types/ball-arithmetic#precision)
discuss how to choose the precision.

*Special values.* The following balls with infinite or NaN parts may be
input or output of any function:

- $[+\infty \pm c]$ and $[-\infty \pm c]$ with $c$ finite represent the
  points $\pm\infty$ (they could be normalised to radius 0, which is not
  done automatically);
- $[c \pm \infty]$ for any $c$ that is not NaN represents the whole
  extended real line $[-\infty, +\infty]$; there is no way to represent a
  half-line such as $[0, +\infty]$;
- a ball with NaN midpoint is *indeterminate*: it may stand for any
  extended real number, or for a function evaluated outside its domain
  (for example at a point where the result would be complex).

*Aliasing.* Every function allows the output to be the same object as
any input of the same type, unless stated otherwise. Identical input
pointers allow, but do not require, algebraic simplification: the
implementation may treat `x` and `x` as the same number rather than two
independent numbers in the same ball. See
[](guide/types/ball-arithmetic#conventions) for the consequences.

## Memory management {#memory}

::: function
void arb_init(arb_t x)

Initialises `x` to the exact value 0.
:::

::: function
void arb_clear(arb_t x)

Frees the memory used by `x`.
:::

::: function
arb_ptr _arb_vec_init(slong n)
void _arb_vec_clear(arb_ptr v, slong n)

Allocate a vector of `n` balls initialised to zero, and free such a
vector.
:::

::: function
void arb_swap(arb_t x, arb_t y)

Swaps `x` and `y` without copying.
:::

::: function
slong arb_allocated_bytes(const arb_t x)

Returns the number of bytes allocated on the heap by `x`, excluding
`sizeof(arb_struct)`.
:::

::: function
slong _arb_vec_allocated_bytes(arb_srcptr vec, slong len)

Returns the number of bytes used by the vector `vec` of length `len`: the
array itself plus the heap memory of its entries.
:::

::: function
double _arb_vec_estimate_allocated_bytes(slong len, slong prec)

Estimates the number of bytes needed for a vector of `len` balls with
`prec`-bit midpoints, including the limb data. The return type is
`double` to avoid overflow.

::: details
This is only an estimate. Zeros and small integers use no heap memory
whatever the precision, the estimate assumes that exponents are small,
and the memory allocator and the operating system add their own overhead.
:::
:::

## Special values {#special}

::: function
void arb_zero(arb_t x)
void arb_one(arb_t f)

Set `x` to the exact value 0, respectively 1.
:::

::: function
void arb_pos_inf(arb_t x)
void arb_neg_inf(arb_t x)

Set `x` to $+\infty$, respectively $-\infty$, with radius 0.
:::

::: function
void arb_zero_pm_inf(arb_t x)

Sets `x` to $[0 \pm \infty]$, the whole extended real line.
:::

::: function
void arb_indeterminate(arb_t x)

Sets `x` to $[\operatorname{NaN} \pm \infty]$, an indeterminate value.
:::

::: function
void arb_zero_pm_one(arb_t x)
void arb_unit_interval(arb_t x)

Set `x` to $[0 \pm 1] = [-1, 1]$, respectively $[1/2 \pm 1/2] = [0, 1]$.
:::

## Assignment and rounding {#assignment}

::: function
void arb_set(arb_t y, const arb_t x)
void arb_set_arf(arb_t x, const arf_t y)
void arb_set_si(arb_t y, slong x)
void arb_set_ui(arb_t y, ulong x)
void arb_set_fmpz(arb_t y, const fmpz_t x)
void arb_set_d(arb_t y, double x)

Set the output to the exact value of the input, without rounding.

::: details
`arb_set_d(x, 1.1)` gives the exact ball $[1.100000000000000088817841970012523\ldots
\pm 0]$, the `double` nearest to 1.1, which does *not* contain $11/10$. To
enclose a decimal number, use `arb_set_str` or divide integers, as
explained in [](arb/printing).
:::
:::

::: function
void arb_set_fmpz_2exp(arb_t x, const fmpz_t y, const fmpz_t exp)

Sets `x` to $y \cdot 2^{\mathit{exp}}$ exactly.
:::

::: function
void arb_set_round(arb_t z, const arb_t x, slong prec)

Sets `z` to `x` with the midpoint rounded to `prec` bits, adding the
rounding error to the radius.
:::

::: function
void arb_set_round_fmpz(arb_t y, const fmpz_t x, slong prec)
void arb_set_round_fmpz_2exp(arb_t y, const fmpz_t x, const fmpz_t exp, slong prec)

Set `y` to a ball containing $x$, respectively $x \cdot 2^{\mathit{exp}}$,
with the midpoint rounded to `prec` bits.
:::

::: function
void arb_set_fmpq(arb_t y, const fmpq_t x, slong prec)

Sets `y` to a ball containing the rational number `x`, with the midpoint
rounded to `prec` bits.
:::

::: function
void arb_trim(arb_t y, const arb_t x)

Sets `y` to a ball containing `x` whose midpoint has only as many bits as
the relative accuracy of `x` justifies, plus 16 guard bits. This saves
memory and time when an inaccurate ball has a long midpoint. A ball with
negative accuracy becomes $[0 \pm (|m| + r)]$.
:::

## Decimal conversion and printing {#io}

::: function
int arb_set_str(arb_t res, const char * inp, slong prec)

Sets `res` to a ball containing the number described by the decimal
string `inp`, rounded to `prec` bits, and returns 0 on success. Accepted
forms are a decimal literal (`"25"`, `"0.001"`, `"7e+141"`,
`"-31.4159e-1"`), a ball `"[3.25 +/- 0.0001]"` (brackets optional) and
`"[+/- 10]"`, as well as `"inf"`, `"-inf"` and `"nan"`.

::: details
The conversion is rigorous: the error of the binary approximation of the
decimal midpoint and radius is added to the radius, so the result contains
the exact decimal value. A NaN midpoint gives an indeterminate ball. On a
parse error the function returns nonzero and sets `res` to an
indeterminate ball. The output of `arb_get_str` with default flags can be
read back with this function.
:::
:::

::: macro
ARB_STR_MORE
ARB_STR_NO_RADIUS
ARB_STR_CONDENSE

Flags for `arb_get_str` and `arb_printn`: print more digits than are
known to be correct, omit the radius, and (as a multiple
$m \cdot$`ARB_STR_CONDENSE`) abbreviate long runs of digits.
:::

::: function
char * arb_get_str(const arb_t x, slong n, ulong flags)

Returns a decimal string for `x`, to be freed with `flint_free`, with at
most `n` significant digits of the midpoint. By default it has the form
`[m +/- r]` with an $m$ that is correct to within one unit in its last
digit and an $r$ such that the decimal ball contains `x`; the string can
be read back with `arb_set_str`.

::: details
With `ARB_STR_MORE`, more (possibly incorrect) digits may be printed.
With `ARB_STR_NO_RADIUS` the radius is omitted; unless `ARB_STR_MORE`
is also set, the midpoint is still rounded so that it is correct to one
unit in the last place, and if no digit is significant the result is
shown as `0e+n`, meaning a number between $-10^n$ and $10^n$. Adding $m$
times `ARB_STR_CONDENSE` to `flags` condenses runs of more than $3m$
digits, printing the first and last $m$ digits and the number of omitted
digits in brackets.

The conversion is done exactly (through the `decimal` module): the
midpoint is correctly rounded to the printed number of digits (to
nearest, with ties away from zero), the number of digits is the largest
for which the printed radius is at most one unit in the last place, and
the printed radius is the true radius plus the exact rounding error,
rounded up to three significant digits. In particular, an exact midpoint
that is exactly representable with `n` digits is printed without a
radius (`arb_get_str` of the exact integer $10^{100}$ with $n = 31$
gives `1.000000000000000000000000000000e+100`), and a tiny radius is never
inflated ($[1 \pm 2^{-1000}]$ prints as
`[1.00000000000000000000000000000 +/- 9.34e-302]` with $n = 30$).
Arguments with huge exponents are handled by scaling with a power of ten
computed in ball arithmetic, which cannot decide whether the value is an
exact decimal number; an exact midpoint with a huge exponent is therefore
tested for exact representability directly.
:::
:::

::: function
void arb_printn(const arb_t x, slong digits, ulong flags)
void arb_fprintn(FILE * file, const arb_t x, slong digits, ulong flags)

Print the string of `arb_get_str` to standard output or to `file`.
This is the recommended way to print a ball.
:::

::: function
void arb_printd(const arb_t x, slong digits)
void arb_fprintd(FILE * file, const arb_t x, slong digits)

Print the midpoint and the radius as decimal floating-point numbers with
`digits` and 5 significant digits, in the form `m +/- r`. The decimal
values are approximations: the printed interval need not contain `x`.
:::

::: function
void arb_print(const arb_t x)
void arb_fprint(FILE * file, const arb_t x)

Print the exact binary representation of `x` as `(m * 2^e) +/- (r * 2^f)`.
:::

::: function
char * arb_dump_str(const arb_t x)
int arb_load_str(arb_t res, const char * data)

`arb_dump_str` returns a string, to be freed with `flint_free`, holding
the exact value of `x` as four hexadecimal integers (midpoint mantissa
and exponent, radius mantissa and exponent, with special encodings for
zero, infinities and NaN) separated by spaces. `arb_load_str` reads it
back exactly and returns nonzero on a format error.
:::

::: function
int arb_dump_file(FILE* stream, const arb_t x)
int arb_load_file(arb_t res, FILE *stream)

Write `x` to `stream` in the format of `arb_dump_str`, and read it back.
Both return nonzero on failure.

::: details
The data is delimited by whitespace or end of file, so values written in
sequence must be separated by whitespace:

```c
fp = fopen("data.txt", "w");
for (i = 0; i < n; i++)
{
    arb_dump_file(fp, vec + i);
    fprintf(fp, "\n");
}
fclose(fp);

fp = fopen("data.txt", "r");
for (i = 0; i < n; i++)
    arb_load_file(vec + i, fp);
fclose(fp);
```
:::
:::

## Midpoint, radius and bounds {#radius}

::: function
void arb_get_mid_arb(arb_t z, const arb_t x)
void arb_get_rad_arb(arb_t z, const arb_t x)

Set `z` to the midpoint, respectively the radius, of `x` as an exact ball.
:::

::: function
void arb_add_error_arf(arb_t x, const arf_t err)
void arb_add_error_mag(arb_t x, const mag_t err)
void arb_add_error(arb_t x, const arb_t error)

Add an upper bound for the absolute value of the error term to the radius
of `x`. This is how a truncation error, computed separately, is attached
to an approximation.
:::

::: function
void arb_add_error_2exp_si(arb_t x, slong err)
void arb_add_error_2exp_fmpz(arb_t x, const fmpz_t err)

Add $2^{\mathit{err}}$ to the radius of `x`.
:::

::: function
void arb_get_mag(mag_t z, const arb_t x)

Sets `z` to an upper bound for $|t|$ for all $t$ in `x`. A NaN midpoint
gives $+\infty$.
:::

::: function
void arb_get_mag_lower(mag_t z, const arb_t x)

Sets `z` to a lower bound for $|t|$ for all $t$ in `x`, which is 0 if `x`
contains zero or NaN.
:::

::: function
void arb_get_mag_lower_nonnegative(mag_t z, const arb_t x)

Sets `z` to a lower bound for all $t$ in `x`, or to 0 if `x` contains
nonpositive numbers or NaN.
:::

::: function
void arb_get_abs_ubound_arf(arf_t u, const arb_t x, slong prec)
void arb_get_abs_lbound_arf(arf_t u, const arb_t x, slong prec)

Set `u` to an upper, respectively lower, bound for $|t|$, $t \in x$,
rounded outward to `prec` bits. NaN gives NaN.
:::

::: function
void arb_get_ubound_arf(arf_t u, const arb_t x, slong prec)
void arb_get_lbound_arf(arf_t u, const arb_t x, slong prec)

Set `u` to the right endpoint $m + r$ rounded up, respectively the left
endpoint $m - r$ rounded down, to `prec` bits. NaN gives NaN.
:::

::: function
void arb_get_interval_arf(arf_t a, arf_t b, const arb_t x, slong prec)
void arb_get_interval_mpfr(mpfr_t a, mpfr_t b, const arb_t x)

Set $[a, b]$ to an interval containing `x`, with the endpoints rounded
outward to `prec` bits, respectively to the precision of the MPFR
variables.
:::

::: function
void arb_set_interval_arf(arb_t x, const arf_t a, const arf_t b, slong prec)
void arb_set_interval_mpfr(arb_t x, const mpfr_t a, const mpfr_t b, slong prec)
void arb_set_interval_mag(arb_t res, const mag_t a, const mag_t b, slong prec)

Set `x` to a ball containing the interval $[a, b]$, with the midpoint
rounded to `prec` bits. Require $a \le b$.
:::

::: function
void arb_set_interval_neg_pos_mag(arb_t res, const mag_t a, const mag_t b, slong prec)

Sets `res` to a ball containing the interval $[-a, b]$.
:::

::: function
void arb_union(arb_t z, const arb_t x, const arb_t y, slong prec)

Sets `z` to a ball containing both `x` and `y`. A NaN input gives an
indeterminate ball and an infinite radius gives $[0 \pm \infty]$.
:::

::: function
int arb_intersection(arb_t z, const arb_t x, const arb_t y, slong prec)

If `x` and `y` overlap (`arb_overlaps`), sets `z` to a ball containing
their intersection and returns nonzero; otherwise returns zero and leaves
`z` undefined. A NaN input gives an indeterminate `z` and a nonzero
return value.
:::

::: function
void arb_nonnegative_part(arb_t res, const arb_t x)

Sets `res` to a ball containing the intersection of `x` with $[0,
+\infty]$. A nonnegative `x` is copied; a finite `x` containing negative
numbers gives a ball $[r/2 \pm r/2]$ with no negative points; a negative
`x` gives 0.
:::

::: function
slong arb_rel_error_bits(const arb_t x)
slong arb_rel_accuracy_bits(const arb_t x)

`arb_rel_error_bits` returns $e_r - e_m + 1$, where $e_m$ and $e_r$ are
the exponents of the midpoint and radius ($2^{e-1} \le |\cdot| < 2^e$):
an estimate of $\log_2(r/|m|)$ rounded up. `arb_rel_accuracy_bits`
returns its negation, the number of correct bits. The value is clamped to
$\pm$`ARF_PREC_EXACT`: an exact ball has accuracy `ARF_PREC_EXACT`, and a
ball with a special midpoint or an infinite radius has accuracy
$-$`ARF_PREC_EXACT`.
:::

::: function
slong arb_rel_one_accuracy_bits(const arb_t x)

Returns the relative accuracy in bits of $[\max(1, |m|) \pm r]$: for balls
near zero, this measures the accuracy relative to 1 instead.
:::

::: function
slong arb_bits(const arb_t x)

Returns the number of bits of the odd mantissa of the midpoint
(`arf_bits`), the precision needed to represent it exactly.
:::

## Integers and rationals {#integers}

::: function
int arb_get_unique_fmpz(fmpz_t z, const arb_t x)

If `x` contains exactly one integer, sets `z` to it and returns nonzero;
otherwise returns zero.

::: details
This is the standard way to recover an integer that is known to be the
exact value of an expression once the ball is narrow enough. A unique
integer of huge size is computed exactly, which can exhaust memory; check
the magnitude of the midpoint first when it is not known to be
reasonable.
:::
:::

::: function
int arb_get_simplest_fmpq(fmpq_t res, const arb_t x)

Sets `res` to the rational number with smallest positive denominator in
the closed interval represented by `x`, and among these the one with
smallest absolute numerator (0 if `x` contains 0), and returns 1. Returns
0 and leaves `res` unchanged if `x` is not finite.

::: details
The rational is found from the continued fraction expansions of the
endpoints. As with `arb_get_unique_fmpz`, very large or very small
exponents of the midpoint or the endpoints can need much memory.
:::
:::

::: function
void arb_get_interval_fmpz_2exp(fmpz_t a, fmpz_t b, fmpz_t exp, const arb_t x)

Writes the exact interval of `x` as $[a, b] \cdot 2^{\mathit{exp}}$ with
integers $a, b$ without common trailing zero bits.

::: details
Raises an exception if `x` is infinite or NaN. The integers can have as
many bits as the difference of the exponents of the midpoint and the
radius, so this can exhaust memory for unusual balls.
:::
:::

::: function
void arb_floor(arb_t z, const arb_t x, slong prec)
void arb_ceil(arb_t z, const arb_t x, slong prec)
void arb_trunc(arb_t res, const arb_t x, slong prec)
void arb_nint(arb_t res, const arb_t x, slong prec)

Set the output to a ball containing $\lfloor t \rfloor$, $\lceil t
\rceil$, $\operatorname{trunc}(t)$, respectively the nearest integer to
$t$ (ties to even), for all $t$ in `x`, with the midpoint rounded to
`prec` bits.
:::

::: function
void arb_get_fmpz_mid_rad_10exp(fmpz_t mid, fmpz_t rad, fmpz_t exp, const arb_t x, slong n)

For finite nonzero `x`, computes integers with $x \subseteq
[\mathit{mid} - \mathit{rad}, \mathit{mid} + \mathit{rad}] \cdot
10^{\mathit{exp}}$ such that the larger of $|\mathit{mid}|$ and
$\mathit{rad}$ has at least `n` digits plus a few guard digits. All three
are set to 0 if `x` is infinite or exactly zero.
:::

::: function
int arb_can_round_arf(const arb_t x, slong prec, arf_rnd_t rnd)
int arb_can_round_mpfr(const arb_t x, slong prec, mpfr_rnd_t rnd)

Return nonzero if rounding the midpoint of `x` to `prec` bits in the
direction `rnd` is guaranteed to give the correctly rounded value of
every real number in `x`.

::: details
If the function returns nonzero, `arf_set_round`, `arf_get_mpfr` (with an
output variable of precision `prec`) or `arf_get_d` (with `prec` = 53)
applied to the midpoint give a correctly rounded result, and
`arf_get_mpfr` also returns the correct MPFR ternary value. This is how
a correctly rounded function can be built on ball arithmetic: evaluate
with a few guard bits, and repeat at higher precision while this test
fails. The test is fast and conservative: it may return 0 when rounding
is in fact possible. It returns 0 for non-finite `x`, even an exact
infinity. The two versions differ only in the type of the rounding
mode.
:::
:::

## Comparison and predicates {#comparison}

All predicates return 1 when the property certainly holds for every
number in the ball(s) and 0 otherwise, which includes the case where it
cannot be decided. [](arb/predicates) explains how to test the
negation of a property correctly.

::: function
int arb_is_zero(const arb_t x)
int arb_is_one(const arb_t f)
int arb_equal_si(const arb_t x, slong y)

Return whether `x` is exactly 0, exactly 1, respectively exactly the
integer `y` (midpoint equal and radius zero).
:::

::: function
int arb_is_exact(const arb_t x)
int arb_is_finite(const arb_t x)

Return whether the radius of `x` is zero, respectively whether both the
midpoint and the radius are finite (not infinite and not NaN).
:::

::: function
int arb_is_int(const arb_t x)
int arb_is_int_2exp_si(const arb_t x, slong e)

Return whether `x` is an exact integer, respectively an exact integer
multiple of $2^e$.
:::

::: function
int arb_is_nonzero(const arb_t x)
int arb_is_positive(const arb_t x)
int arb_is_nonnegative(const arb_t x)
int arb_is_negative(const arb_t x)
int arb_is_nonpositive(const arb_t x)

Return whether all points $t$ of `x` satisfy $t \ne 0$, $t > 0$, $t \ge
0$, $t < 0$, respectively $t \le 0$. They return 0 if `x` contains NaN.
:::

::: function
int arb_contains_zero(const arb_t x)
int arb_contains_negative(const arb_t x)
int arb_contains_nonpositive(const arb_t x)
int arb_contains_positive(const arb_t x)
int arb_contains_nonnegative(const arb_t x)

Return whether some point $t$ of `x` satisfies $t = 0$, $t < 0$, $t \le
0$, $t > 0$, respectively $t \ge 0$. They return 1 if `x` contains NaN.
:::

::: function
int arb_contains_int(const arb_t x)

Returns whether `x` contains an integer.
:::

::: function
int arb_eq(const arb_t x, const arb_t y)
int arb_ne(const arb_t x, const arb_t y)
int arb_lt(const arb_t x, const arb_t y)
int arb_le(const arb_t x, const arb_t y)
int arb_gt(const arb_t x, const arb_t y)
int arb_ge(const arb_t x, const arb_t y)

Return 1 if $t \;\mathrm{op}\; u$ holds for all $t \in x$ and $u \in y$,
where op is $=$, $\ne$, $<$, $\le$, $>$, respectively $\ge$, and 0
otherwise. They return 0 if either midpoint is NaN.

::: details
The balls are viewed as subsets of the extended real line, so formally
different balls can be equal: $[-\infty \pm 3] = [-\infty \pm 0]$, and
$[-\infty] \le [\infty \pm \infty]$. `arb_eq` is true only when both
balls are the same single point.
:::
:::

::: function
int arb_equal(const arb_t x, const arb_t y)

Returns whether `x` and `y` are identical as balls: the same midpoint
and the same radius.

::: details
This is a test of the representation, not of the numbers: two equal
inexact balls need not represent the same number, and two different balls
can. Use `arb_eq` to prove that two numbers are equal, and `arb_overlaps`
or `arb_contains` to test whether they might be.
:::
:::

::: function
int arb_overlaps(const arb_t x, const arb_t y)

Returns whether `x` and `y` have a point in common. Returns 1 if either
contains NaN.
:::

::: function
int arb_contains(const arb_t x, const arb_t y)
int arb_contains_arf(const arb_t x, const arf_t y)
int arb_contains_fmpq(const arb_t x, const fmpq_t y)
int arb_contains_fmpz(const arb_t x, const fmpz_t y)
int arb_contains_si(const arb_t x, slong y)
int arb_contains_mpfr(const arb_t x, const mpfr_t y)

Return whether the ball or number `y` is contained in `x`. If `x`
contains NaN the result is 1; otherwise, if `y` is NaN the result is 0.
:::

::: function
int arb_contains_interior(const arb_t x, const arb_t y)

Returns whether `y` is contained in the interior of `x`, that is,
contained in `x` without touching either endpoint.
:::

## Arithmetic {#arithmetic}

The basic operations compute the midpoint with one `arf` operation and
bound the propagated error with a few `mag` operations, so they are
nearly as fast as floating-point arithmetic of the same precision. The
error bounds are proved in [](algorithms/ball-arithmetic/error-propagation).

::: function
void arb_neg(arb_t y, const arb_t x)
void arb_neg_round(arb_t y, const arb_t x, slong prec)

Set `y` to $-x$, exactly, respectively with the midpoint rounded to
`prec` bits.
:::

::: function
void arb_abs(arb_t y, const arb_t x)

Sets `y` to $[|m| \pm r]$. This contains $|t|$ for all $t \in x$; if `x`
contains zero the result also contains negative numbers.
:::

::: function
void arb_nonnegative_abs(arb_t y, const arb_t x)

Sets `y` to a ball containing $|t|$ for all $t \in x$ that contains no
negative numbers: $[|m| \pm r]$ if $|m| \ge r$ and a ball $[s/2 \pm
s/2]$ with $s \ge |m| + r$ otherwise.
:::

::: function
void arb_sgn(arb_t res, const arb_t x)

Sets `res` to a ball containing the sign of every point of `x`: $\pm 1$
or 0 when the sign is determined, and $[0 \pm 1]$ otherwise.
:::

::: function
int arb_sgn_nonzero(const arb_t x)

Returns 1 if `x` is certainly positive, $-1$ if it is certainly negative,
and 0 otherwise.
:::

::: function
void arb_min(arb_t z, const arb_t x, const arb_t y, slong prec)
void arb_max(arb_t z, const arb_t x, const arb_t y, slong prec)
void arb_minmax(arb_t z1, arb_t z2, const arb_t x, const arb_t y, slong prec)

Set the output to a ball containing $\min(t, u)$, respectively $\max(t,
u)$, for all $t \in x$, $u \in y$; `arb_minmax` computes both.
:::

::: function
void arb_add(arb_t z, const arb_t x, const arb_t y, slong prec)
void arb_add_arf(arb_t z, const arb_t x, const arf_t y, slong prec)
void arb_add_ui(arb_t z, const arb_t x, ulong y, slong prec)
void arb_add_si(arb_t z, const arb_t x, slong y, slong prec)
void arb_add_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec)

Set `z` to $x + y$. The midpoint is $m_x + m_y$ rounded to `prec` bits
and the radius is $r_x + r_y$ plus the rounding error. `prec` may be
`ARF_PREC_EXACT` if the exact sum fits in memory.
:::

::: function
void arb_add_fmpz_2exp(arb_t z, const arb_t x, const fmpz_t man, const fmpz_t exp, slong prec)

Sets `z` to $x + \mathit{man} \cdot 2^{\mathit{exp}}$.
:::

::: function
void arb_sub(arb_t z, const arb_t x, const arb_t y, slong prec)
void arb_sub_arf(arb_t z, const arb_t x, const arf_t y, slong prec)
void arb_sub_ui(arb_t z, const arb_t x, ulong y, slong prec)
void arb_sub_si(arb_t z, const arb_t x, slong y, slong prec)
void arb_sub_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec)

Set `z` to $x - y$. The radius is $r_x + r_y$ plus the rounding error,
also when `x` and `y` are the same object.
:::

::: function
void arb_mul(arb_t z, const arb_t x, const arb_t y, slong prec)
void arb_mul_arf(arb_t z, const arb_t x, const arf_t y, slong prec)
void arb_mul_si(arb_t z, const arb_t x, slong y, slong prec)
void arb_mul_ui(arb_t z, const arb_t x, ulong y, slong prec)
void arb_mul_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec)

Set `z` to $xy$. The midpoint is $m_x m_y$ rounded to `prec` bits and the
radius is $|m_x| r_y + |m_y| r_x + r_x r_y$ plus the rounding error.
:::

::: function
void arb_sqr(arb_t res, const arb_t val, slong prec)

Sets `res` to $x^2$. Currently this is `arb_mul(res, val, val, prec)`,
which uses the same bound as a product of two independent balls.
:::

::: function
void arb_mul_2exp_si(arb_t y, const arb_t x, slong e)
void arb_mul_2exp_fmpz(arb_t y, const arb_t x, const fmpz_t e)

Set `y` to $x \cdot 2^e$. Exact.
:::

::: function
void arb_addmul(arb_t z, const arb_t x, const arb_t y, slong prec)
void arb_addmul_arf(arb_t z, const arb_t x, const arf_t y, slong prec)
void arb_addmul_si(arb_t z, const arb_t x, slong y, slong prec)
void arb_addmul_ui(arb_t z, const arb_t x, ulong y, slong prec)
void arb_addmul_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec)
void arb_submul(arb_t z, const arb_t x, const arb_t y, slong prec)
void arb_submul_arf(arb_t z, const arb_t x, const arf_t y, slong prec)
void arb_submul_si(arb_t z, const arb_t x, slong y, slong prec)
void arb_submul_ui(arb_t z, const arb_t x, ulong y, slong prec)
void arb_submul_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec)

Set `z` to $z + xy$, respectively $z - xy$, with a single rounding of
the midpoint.
:::

::: function
void arb_fma(arb_t res, const arb_t x, const arb_t y, const arb_t z, slong prec)
void arb_fma_arf(arb_t res, const arb_t x, const arf_t y, const arb_t z, slong prec)
void arb_fma_si(arb_t res, const arb_t x, slong y, const arb_t z, slong prec)
void arb_fma_ui(arb_t res, const arb_t x, ulong y, const arb_t z, slong prec)
void arb_fma_fmpz(arb_t res, const arb_t x, const fmpz_t y, const arb_t z, slong prec)

Set `res` to $xy + z$ with a single rounding of the midpoint; unlike
`arb_addmul`, the output may be a different object from `z`.
:::

::: function
void arb_inv(arb_t y, const arb_t x, slong prec)

Sets `y` to $1/x$. If `x` contains zero, the result is indeterminate.
:::

::: function
void arb_div(arb_t z, const arb_t x, const arb_t y, slong prec)
void arb_div_arf(arb_t z, const arb_t x, const arf_t y, slong prec)
void arb_div_si(arb_t z, const arb_t x, slong y, slong prec)
void arb_div_ui(arb_t z, const arb_t x, ulong y, slong prec)
void arb_div_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec)
void arb_fmpz_div_fmpz(arb_t z, const fmpz_t x, const fmpz_t y, slong prec)
void arb_ui_div(arb_t z, ulong x, const arb_t y, slong prec)

Set `z` to $x / y$. If `y` contains zero, the result is indeterminate
($[\operatorname{NaN} \pm \infty]$), even when $x = 0$.

::: details
The midpoint is computed with `arf_div` in the relaxed mode
`ARF_RND_FAST`, which allows an error of one unit in the last place and
uses Newton division at high precision. For balls $[m_x \pm r_x]$ and
$[m_y \pm r_y]$ with $r_y < |m_y|$, the propagated error is bounded by

$$
\frac{|m_x| r_y + |m_y| r_x}{|m_y| (|m_y| - r_y)},
$$

which follows from the triangle inequality in the numerator and the
reverse triangle inequality in the denominator of
$x/y - (x + a)/(y + b)$ with $|a| \le r_x$, $|b| \le r_y$. For wide balls
the endpoints are used instead to get a tighter enclosure.

The inexact flag of the relaxed division is reliable also on the Newton
path (from 12000 bits, or 5000 bits for $1/y$, when the divisor has
more than `prec`/2 bits), so one ulp is added to the radius whenever
the midpoint differs from the exact quotient, and an exact quotient of
exact input gets radius 0.
:::
:::

::: function
void arb_div_2expm1_ui(arb_t z, const arb_t x, ulong n, slong prec)

Sets `z` to $x / (2^n - 1)$, using the expansion $2^{-n} + 2^{-2n} +
\cdots$ for large $n$.
:::

## Dot products {#dot}

::: function
void arb_dot(arb_t res, const arb_t initial, int subtract, arb_srcptr x, slong xstep, arb_srcptr y, slong ystep, slong len, slong prec)
void arb_dot_simple(arb_t res, const arb_t initial, int subtract, arb_srcptr x, slong xstep, arb_srcptr y, slong ystep, slong len, slong prec)
void arb_dot_precise(arb_t res, const arb_t initial, int subtract, arb_srcptr x, slong xstep, arb_srcptr y, slong ystep, slong len, slong prec)

Set `res` to $s + (-1)^{\mathit{subtract}} \sum_{i=0}^{\mathit{len}-1}
x_i y_i$, where $s$ is `initial` (or 0 if `initial` is `NULL`) and the
entries are read with strides `xstep` and `ystep`.

::: details
`subtract` must be 0 or 1. A negative `len` counts as 0, and the strides
may be negative. `res` may alias `initial` but not an entry of `x` or `y`.

`arb_dot` is the function to use. It accumulates all products in one
fixed-point sum with a single final rounding and computes all radius
contributions with `mag` arithmetic in one pass, so it is both faster and
more accurate than a loop of `arb_addmul`: the internal sum carries
$4 + \log_2 \mathit{len}$ guard bits, so the midpoint error is close to
that of a single rounding of the exact dot product, and the midpoint is
not computed to more bits than the radii justify.
`arb_dot_simple` is that loop, kept for testing and used as a fallback
for exponents outside the range of the fast code. `arb_dot_precise`
computes the exact dot product of the midpoints before rounding, which
can be very slow; it is meant for testing. The algorithm is described in
[](algorithms/ball-arithmetic/dot-products).
:::
:::

::: function
void arb_approx_dot(arb_t res, const arb_t initial, int subtract, arb_srcptr x, slong xstep, arb_srcptr y, slong ystep, slong len, slong prec)

Computes the same dot product without error bounds: only the midpoints
of the inputs are read and only the midpoint of `res` is written, using
`arf_approx_dot`. Useful for approximate (floating-point) algorithms such
as iterative refinement.
:::

::: function
void arb_dot_ui(arb_t res, const arb_t initial, int subtract, arb_srcptr x, slong xstep, const ulong * y, slong ystep, slong len, slong prec)
void arb_dot_si(arb_t res, const arb_t initial, int subtract, arb_srcptr x, slong xstep, const slong * y, slong ystep, slong len, slong prec)
void arb_dot_uiui(arb_t res, const arb_t initial, int subtract, arb_srcptr x, slong xstep, const ulong * y, slong ystep, slong len, slong prec)
void arb_dot_siui(arb_t res, const arb_t initial, int subtract, arb_srcptr x, slong xstep, const ulong * y, slong ystep, slong len, slong prec)
void arb_dot_fmpz(arb_t res, const arb_t initial, int subtract, arb_srcptr x, slong xstep, const fmpz * y, slong ystep, slong len, slong prec)

Like `arb_dot` with an integer vector `y`. In the `uiui` and `siui`
versions each $y_i$ is a two-limb integer stored as two consecutive
`ulong` (low limb first), unsigned, respectively signed in two's
complement; `ystep` counts integers, not limbs.
:::

## Powers and roots {#powers}

::: function
void arb_sqrt(arb_t z, const arb_t x, slong prec)
void arb_sqrt_arf(arb_t z, const arf_t x, slong prec)
void arb_sqrt_fmpz(arb_t z, const fmpz_t x, slong prec)
void arb_sqrt_ui(arb_t z, ulong x, slong prec)

Set `z` to $\sqrt{x}$. If `x` contains negative numbers the result is
indeterminate; use `arb_sqrtpos` when the negative part is known to be
spurious.

::: details
For $x = [m \pm r]$ with $0 \le r \le m$ the propagated error is bounded
by $\sqrt{m} - \sqrt{m - r} \le \frac{\sqrt{m}}{2}(r/m + (r/m)^2)$. The
midpoint is computed with `arf_sqrt` in the relaxed mode `ARF_RND_FAST`
(Newton iteration from 100000 bits); as for `arb_div`, its inexact flag
is reliable, so exact input gives radius 0 exactly when the square root
is exact.
:::
:::

::: function
void arb_sqrtpos(arb_t z, const arb_t x, slong prec)

Sets `z` to a ball containing $\sqrt{t}$ for all nonnegative $t \in x$,
ignoring the negative part of `x`.
:::

::: function
void arb_hypot(arb_t z, const arb_t x, const arb_t y, slong prec)

Sets `z` to $\sqrt{x^2 + y^2}$.
:::

::: function
void arb_rsqrt(arb_t z, const arb_t x, slong prec)
void arb_rsqrt_ui(arb_t z, ulong x, slong prec)

Set `z` to $1/\sqrt{x}$. At high precision this is faster than computing
a square root.
:::

::: function
void arb_sqrt1pm1(arb_t r, const arb_t z, slong prec)

Sets `r` to $\sqrt{1 + z} - 1$, accurately also for $z$ near 0.
:::

::: function
void arb_root_ui(arb_t z, const arb_t x, ulong k, slong prec)
void arb_root(arb_t z, const arb_t x, ulong k, slong prec)

Set `z` to the principal `k`-th root of `x`. `arb_root` is an alias kept
for compatibility.

::: details
For small `k`, the root of the midpoint is computed with `arf_root` and
for $x = [m \pm r]$, $r < m$, the propagated error is bounded by

$$
m^{1/k} - (m - r)^{1/k} \le m^{1/k} \min\left(1, \frac{1}{k}
\log\left(1 + \frac{r}{m - r}\right)\right),
$$

evaluated with `mag_log1p`. For large `k`, the root is computed as
$\exp(\log(x)/k)$.
:::
:::

::: function
void arb_pow_ui(arb_t y, const arb_t b, ulong e, slong prec)
void arb_pow_fmpz(arb_t y, const arb_t b, const fmpz_t e, slong prec)
void arb_pow_fmpz_binexp(arb_t y, const arb_t b, const fmpz_t e, slong prec)
void arb_ui_pow_ui(arb_t y, ulong b, ulong e, slong prec)
void arb_si_pow_ui(arb_t y, slong b, ulong e, slong prec)

Set `y` to $b^e$ by binary exponentiation, with an initial inversion if
$e < 0$. With exact `b` and $e \ge 0$, `prec` may be `ARF_PREC_EXACT` to
get the exact power if it fits in memory.

::: details
The cost grows with the size of $e$; for huge exponents, `arb_pow` (which
uses $\exp(e \log b)$) can be much faster.
:::
:::

::: function
void arb_pow_fmpq(arb_t y, const arb_t x, const fmpq_t a, slong prec)

Sets `y` to $x^a$, computed as $(x^{1/q})^p$ if the denominator $q$ of
$a = p/q$ is small and as $\exp(a \log x)$ otherwise.
:::

::: function
void arb_pow(arb_t z, const arb_t x, const arb_t y, slong prec)

Sets `z` to $x^y$: by binary exponentiation if `y` is a small exact
integer, as $(x^{1/2})^{2y}$ if `y` is a small exact half-integer, and as
$\exp(y \log x)$ otherwise, with $x^y = \pm(-x)^y$ for negative $x$ and
exact integer $y$.

::: details
$x^0 = 1$ for every $x$, and $0^y$ is 0 for positive $y$ and
indeterminate otherwise. For a ball $x = [a \pm a]$ with $a > 0$, which
contains 0, and nonnegative $y$, the result is the finite ball $[u/2 \pm
u/2]$ with $u \ge (2a)^y$ rather than an indeterminate value. Otherwise a
base containing nonpositive numbers gives an indeterminate result unless
$y$ is an exact integer.
:::
:::

## Exponentials and logarithms {#exp-log}

The guide [](arb/elementary) explains the accuracy and the
special cases of the elementary functions; the algorithms are described in
[](algorithms/special-functions/elementary).

::: function
void arb_exp(arb_t z, const arb_t x, slong prec)

Sets `z` to $e^x$.

::: details
For $x = [m \pm r]$, the midpoint is $e^m$ computed with `arb_exp_arf`
and the propagated error is bounded by $e^m (e^r - 1)$, evaluated with
`mag_expm1`; wide balls are evaluated at the endpoints. If the input is
inexact, $e^m$ is computed only to about the accuracy of the input plus
30 bits. If $|m| > 2^{\max(128, 2\,\mathit{prec})}$, the result is
$[0 \pm \infty]$ for positive $m$ and a tiny ball near 0 for negative
$m$, so that the cost stays polynomial in `prec`.
:::
:::

::: function
void arb_expm1(arb_t z, const arb_t x, slong prec)

Sets `z` to $e^x - 1$, with small relative error also near $x = 0$.
:::

::: function
void arb_exp_invexp(arb_t z, arb_t w, const arb_t x, slong prec)

Sets `z` to $e^x$ and `w` to $e^{-x}$. The second is computed from the
first by a division, with its propagated error bounded separately.
:::

::: function
void arb_log(arb_t z, const arb_t x, slong prec)
void arb_log_arf(arb_t z, const arf_t x, slong prec)
void arb_log_ui(arb_t z, ulong x, slong prec)
void arb_log_fmpz(arb_t z, const fmpz_t x, slong prec)

Set `z` to $\log x$. If `x` contains nonpositive numbers, the result is
indeterminate.

::: details
`arb_log` computes the logarithm of the midpoint with `arb_log_arf` and
bounds the propagated error for $x = [m \pm r]$, $r < m$, by
$\log(1 + r/(m - r))$. Up to about 2800 bits, `arb_log_arf` reduces the
argument with precomputed tables and sums an arctangent series by
rectangular splitting; above, it uses a Newton-type iteration on
`arb_exp`, and logarithms of smooth integers use precomputed logarithms
of primes.
:::
:::

::: function
void arb_log_ui_from_prev(arb_t s, ulong k, arb_t log_prev, ulong prev, slong prec)

Sets `s` to $\log k$, given `log_prev` $= \log(\mathit{prev})$ with
$\mathit{prev} < k$. At high precision it uses $\log k = \log
\mathit{prev} + 2 \operatorname{atanh}((k - \mathit{prev}) / (k +
\mathit{prev}))$ with the atanh series evaluated by binary splitting,
which is fast when $\mathit{prev}$ is large and $k - \mathit{prev}$
small; otherwise it ignores `log_prev`.
:::

::: function
void arb_log1p(arb_t r, const arb_t z, slong prec)

Sets `r` to $\log(1 + z)$, with small relative error also near $z = 0$.
:::

::: function
void arb_log_base_ui(arb_t res, const arb_t x, ulong b, slong prec)

Sets `res` to $\log_b x$. The result is exact when $x$ is an exact power
of $b$.
:::

::: function
void arb_log_hypot(arb_t res, const arb_t x, const arb_t y, slong prec)

Sets `res` to $\log \sqrt{x^2 + y^2}$, accurately also when $x^2 + y^2$
is close to 1. This is the real part of the complex logarithm.
:::

## Trigonometric functions {#trig}

::: function
void arb_sin(arb_t s, const arb_t x, slong prec)
void arb_cos(arb_t c, const arb_t x, slong prec)
void arb_sin_cos(arb_t s, arb_t c, const arb_t x, slong prec)

Set `s` to $\sin x$ and `c` to $\cos x$. Computing both together costs
little more than one.

::: details
The propagated error for a radius $r$ is at most $\min(2, r)$. Balls with
a radius larger than about $2^{-24}$ are handled by `arb_sin_cos_wide`,
which computes a near-optimal enclosure of the range. If the midpoint
exceeds $2^{\max(65536, 4\,\mathit{prec})}$, the result is $[0 \pm 1]$:
reducing such an argument modulo $\pi$ would take time exponential in the
size of the input.
:::
:::

::: function
void arb_sin_pi(arb_t s, const arb_t x, slong prec)
void arb_cos_pi(arb_t c, const arb_t x, slong prec)
void arb_sin_cos_pi(arb_t s, arb_t c, const arb_t x, slong prec)

Set `s` to $\sin(\pi x)$ and `c` to $\cos(\pi x)$. These are exact at
exact integers and half-integers and accurate near the zeros, which
`arb_sin` applied to a multiple of an enclosure of $\pi$ is not.
:::

::: function
void arb_tan(arb_t y, const arb_t x, slong prec)
void arb_cot(arb_t y, const arb_t x, slong prec)

Set `y` to $\tan x = \sin x / \cos x$, respectively $\cot x = \cos x /
\sin x$.
:::

::: function
void arb_tan_pi(arb_t y, const arb_t x, slong prec)
void arb_cot_pi(arb_t y, const arb_t x, slong prec)

Set `y` to $\tan(\pi x)$, respectively $\cot(\pi x)$.
:::

::: function
void arb_sin_cos_pi_fmpq(arb_t s, arb_t c, const fmpq_t x, slong prec)
void arb_sin_pi_fmpq(arb_t s, const fmpq_t x, slong prec)
void arb_cos_pi_fmpq(arb_t c, const fmpq_t x, slong prec)

Set `s` to $\sin(\pi x)$ and `c` to $\cos(\pi x)$ for a rational $x$ in
lowest terms.

::: details
Symmetries first reduce $x$ to $[0, 1/4]$. Then either $\pi x$ is
approximated and the usual functions are called, or, for small
denominators, the value is computed as an algebraic number (a root of a
Chebyshev-type polynomial), whichever is estimated to be faster. Since the
reduced argument is in the first octant, the first method gives full
relative accuracy even near the zeros.
:::
:::

::: function
void arb_sec(arb_t res, const arb_t x, slong prec)
void arb_csc(arb_t res, const arb_t x, slong prec)
void arb_csc_pi(arb_t res, const arb_t x, slong prec)

Set `res` to $\sec x = 1/\cos x$, $\csc x = 1/\sin x$, respectively
$\csc(\pi x)$.
:::

::: function
void arb_sinc(arb_t z, const arb_t x, slong prec)
void arb_sinc_pi(arb_t z, const arb_t x, slong prec)

Set `z` to $\operatorname{sinc} x = \sin(x)/x$, respectively
$\operatorname{sinc}(\pi x) = \sin(\pi x)/(\pi x)$, with the value 1 at
$x = 0$ and accurate enclosures for balls containing 0.
:::

## Inverse trigonometric functions {#inverse-trig}

::: function
void arb_atan(arb_t z, const arb_t x, slong prec)
void arb_atan_arf(arb_t z, const arf_t x, slong prec)

Set `z` to $\operatorname{atan} x$.

::: details
`arb_atan` computes the arctangent of the midpoint with `arb_atan_arf`
and adds the propagated error, bounded by $r / (1 + (|m| - r)^2)$ for $r <
|m|$ and by $\min(\pi, r)$ in general; inaccurate input is evaluated at
the endpoints. `arb_atan_arf` uses tables and a Taylor series below 3400
bits and a Newton iteration on `arb_sin_cos` above; the bit-burst
algorithm `arb_atan_arf_bb` is the fallback when the tables are not
precise enough.
:::
:::

::: function
void arb_atan2(arb_t z, const arb_t b, const arb_t a, slong prec)

Sets `z` to the argument (phase) of the complex number $a + bi$, in $(-\pi,
\pi]$, with the branch cut on $(-\infty, 0]$. By definition
$\operatorname{atan2}(0, 0) = 0$ and $\operatorname{atan2}(0, a) = \pi$
for $a < 0$. A ball that straddles the branch cut gives an enclosure of
$[-\pi, \pi]$.
:::

::: function
void arb_asin(arb_t z, const arb_t x, slong prec)
void arb_acos(arb_t z, const arb_t x, slong prec)

Set `z` to $\operatorname{asin} x = \operatorname{atan}(x / \sqrt{1 -
x^2})$, respectively $\operatorname{acos} x = \pi/2 - \operatorname{asin}
x$. If `x` is not contained in $[-1, 1]$, the result is indeterminate.
:::

## Hyperbolic functions {#hyperbolic}

::: function
void arb_sinh(arb_t z, const arb_t x, slong prec)
void arb_cosh(arb_t z, const arb_t x, slong prec)
void arb_sinh_cosh(arb_t s, arb_t c, const arb_t x, slong prec)

Set `s` to $\sinh x$ and `c` to $\cosh x$. Near 0, $\sinh x$ is computed
as $(e^{2x} - 1)/(2e^x)$ via `arb_expm1` to avoid cancellation; otherwise
$(e^x \pm e^{-x})/2$ is used.
:::

::: function
void arb_tanh(arb_t y, const arb_t x, slong prec)
void arb_coth(arb_t y, const arb_t x, slong prec)

Set `y` to $\tanh x$, respectively $\coth x$, computed as $(e^{2x} -
1)/(e^{2x} + 1)$ via `arb_expm1` for small $|x|$ and as $\pm(1 - 2e^{\mp
2x}/(1 + e^{\mp 2x}))$ for large $|x|$ (and the reciprocals).
:::

::: function
void arb_sech(arb_t res, const arb_t x, slong prec)
void arb_csch(arb_t res, const arb_t x, slong prec)

Set `res` to $\operatorname{sech} x = 1/\cosh x$, respectively
$\operatorname{csch} x = 1/\sinh x$.
:::

## Inverse hyperbolic functions {#inverse-hyperbolic}

::: function
void arb_atanh(arb_t z, const arb_t x, slong prec)
void arb_asinh(arb_t z, const arb_t x, slong prec)

Set `z` to $\operatorname{atanh} x$, respectively $\operatorname{asinh}
x$. For `arb_atanh`, an input not contained in $(-1, 1)$ gives an
indeterminate result.
:::

::: function
void arb_acosh(arb_t z, const arb_t x, slong prec)

Sets `z` to $\operatorname{acosh} x$. If `x` contains numbers less than
1, the result is indeterminate.
:::

## Constants {#constants}

Each function returns a ball of radius about one unit in the last place
at precision `prec`. Values are cached per thread ($\pi$ and $\log 2$
are read from static tables below about 4600 bits and $\gamma$ below
about 3400 bits), so
repeated calls at the same or a lower precision cost a rounding;
`flint_cleanup` frees the caches. The formulas and their cost are
described in [](algorithms/special-functions/constants).

::: function
void arb_const_pi(arb_t z, slong prec)
void arb_const_sqrt_pi(arb_t z, slong prec)
void arb_const_log_sqrt2pi(arb_t z, slong prec)

Set `z` to $\pi$, $\sqrt{\pi}$, respectively $\log \sqrt{2\pi}$.
:::

::: function
void arb_const_e(arb_t z, slong prec)

Sets `z` to $e = \exp(1)$.
:::

::: function
void arb_const_log2(arb_t z, slong prec)
void arb_const_log10(arb_t z, slong prec)

Set `z` to $\log 2$, respectively $\log 10$.
:::

::: function
void arb_const_euler(arb_t z, slong prec)

Sets `z` to Euler's constant $\gamma = \lim_{k \to \infty} (H_k - \log
k)$, where $H_k = 1 + 1/2 + \cdots + 1/k$.
:::

::: function
void arb_const_catalan(arb_t z, slong prec)

Sets `z` to Catalan's constant $G = \sum_{n=0}^{\infty} (-1)^n /
(2n+1)^2$.
:::

::: function
void arb_const_apery(arb_t s, slong prec)
void arb_const_zeta5(arb_t s, slong prec)

Set `s` to Apéry's constant $\zeta(3)$, respectively $\zeta(5)$.
:::

::: function
void arb_const_khinchin(arb_t z, slong prec)

Sets `z` to Khinchin's constant $K_0$.
:::

::: function
void arb_const_glaisher(arb_t z, slong prec)

Sets `z` to the Glaisher-Kinkelin constant $A = \exp(1/12 -
\zeta'(-1))$, with the derivative from `acb_dirichlet_zeta_jet`: by
Euler-Maclaurin summation, or at very high precision by the approximate
functional equation (`acb_dirichlet_l_fmpq_afe_jet`) in time
$\tilde O(p^{1.5})$.
:::

::: function
void arb_const_reciprocal_fibonacci(arb_t res, slong prec)

Sets `res` to the reciprocal Fibonacci constant $\sum_{n=1}^{\infty}
1/F_n$.
:::

## Gamma function, factorials and binomial coefficients {#gamma}

Most of these functions are wrappers of `arb_hypgeom` functions; see
[](algorithms/special-functions/gamma) for the algorithms and
[](arb_hypgeom/gamma) for how the method is chosen.

::: function
void arb_gamma(arb_t z, const arb_t x, slong prec)
void arb_gamma_fmpq(arb_t z, const fmpq_t x, slong prec)
void arb_gamma_fmpz(arb_t z, const fmpz_t x, slong prec)

Set `z` to $\Gamma(x)$. Aliases of `arb_hypgeom_gamma`,
`arb_hypgeom_gamma_fmpq` and `arb_hypgeom_gamma_fmpz`.
:::

::: function
void arb_lgamma(arb_t z, const arb_t x, slong prec)

Sets `z` to $\log \Gamma(x)$. The branch structure of the complex log-gamma
function is assumed, so the result is indeterminate if `x` contains
nonpositive numbers. Alias of `arb_hypgeom_lgamma`.
:::

::: function
void arb_rgamma(arb_t z, const arb_t x, slong prec)

Sets `z` to $1/\Gamma(x)$, which is finite (zero) at the poles of
$\Gamma$. Alias of `arb_hypgeom_rgamma`.
:::

::: function
void arb_digamma(arb_t y, const arb_t x, slong prec)

Sets `y` to the digamma function $\psi(x) = \Gamma'(x)/\Gamma(x)$.
:::

::: function
void arb_rising_ui(arb_t z, const arb_t x, ulong n, slong prec)
void arb_rising(arb_t z, const arb_t x, const arb_t n, slong prec)

Set `z` to the rising factorial $x (x+1) \cdots (x+n-1) =
\Gamma(x+n)/\Gamma(x)$. Aliases of `arb_hypgeom_rising_ui` and
`arb_hypgeom_rising`.
:::

::: function
void arb_rising_fmpq_ui(arb_t y, const fmpq_t x, ulong n, slong prec)

Sets `y` to the rising factorial of the rational number `x`, by binary
splitting. If the numerator or denominator of `x` is large compared to
`prec`, `arb_rising_ui` on an approximation of `x` is faster.
:::

::: function
void arb_rising2_ui(arb_t u, arb_t v, const arb_t x, ulong n, slong prec)

Sets `u` to the rising factorial $u(x) = x (x+1) \cdots (x+n-1)$ and `v`
to its derivative $u'(x)$. Wrapper of `arb_hypgeom_rising_ui_jet`.
:::

::: function
void arb_fac_ui(arb_t z, ulong n, slong prec)
void arb_doublefac_ui(arb_t z, ulong n, slong prec)

Set `z` to $n!$, respectively the double factorial $n!!$.
:::

::: function
void arb_bin_ui(arb_t z, const arb_t n, ulong k, slong prec)
void arb_bin_uiui(arb_t z, ulong n, ulong k, slong prec)

Set `z` to the binomial coefficient $\binom{n}{k} = (n-k+1)_k / k!$.
:::

## Zeta function {#zeta}

The algorithms for $\zeta(s)$ at integers and in general are described in
[](algorithms/special-functions/zeta).

::: function
void arb_zeta(arb_t z, const arb_t s, slong prec)
void arb_hurwitz_zeta(arb_t z, const arb_t s, const arb_t a, slong prec)

Set `z` to the Riemann zeta function $\zeta(s)$, respectively the Hurwitz
zeta function $\zeta(s, a)$. For derivatives with respect to $s$, use
`arb_poly_zeta_series`.
:::

::: function
void arb_zeta_ui(arb_t x, ulong n, slong prec)

Sets `x` to $\zeta(n)$ for an integer $n \ge 0$, $n \ne 1$, choosing an
algorithm automatically. For many values, the vector functions below are
more efficient.
:::

::: function
void arb_zeta_ui_vec(arb_ptr x, ulong start, slong num, slong prec)
void arb_zeta_ui_vec_even(arb_ptr x, ulong start, slong num, slong prec)
void arb_zeta_ui_vec_odd(arb_ptr x, ulong start, slong num, slong prec)

Set `x` to $\zeta(s)$ for `num` consecutive integers (respectively even
or odd integers) $s$ starting at `start` $\ge 2$, choosing an algorithm
automatically.
:::

::: function
void arb_zeta_ui_vec_borwein(arb_ptr z, ulong start, slong num, ulong step, slong prec)

Sets `z` to $\zeta(s)$ for $s = \mathit{start}, \mathit{start} +
\mathit{step}, \ldots$ (`num` values), $\mathit{start} \ge 2$, by Borwein's
algorithm [@Bor2000; @GS2003], organised for multi-evaluation.

::: details
For efficiency the largest $s$ should be at most about `prec`; beyond
that, the Euler product is better. Arguments near `WORD_MAX` overflow.
For a single $s$ the algorithm is essentially that of MPFR [@MPFR2012]:
the sum is evaluated backwards to avoid storing more than one coefficient
$d_k$, in integer arithmetic since the terms are slightly larger than
$2^{\mathit{prec}}$. The only error in the main loop comes from the
division by $k^s$, at most one unit per term; for multi-evaluation the
terms are divided repeatedly by $k^{\mathit{step}}$, and by induction the
error per term stays below two units.
:::
:::

::: function
void arb_zeta_ui_borwein_bsplit(arb_t x, ulong s, slong prec)

Sets `x` to $\zeta(s)$, $s \ge 2$, by binary splitting applied to
Borwein's algorithm, in quasilinear time in `prec` for fixed $s$.
:::

::: function
void arb_zeta_ui_euler_product(arb_t z, ulong s, slong prec)
void arb_zeta_ui_asymp(arb_t x, ulong s, slong prec)

Set the output to $\zeta(s)$ using the Euler product, which is fast only
when $s$ is large compared to `prec`. Both are wrappers of
`_acb_dirichlet_euler_product_real_ui`.
:::

::: function
void arb_zeta_ui_bernoulli(arb_t x, ulong n, slong prec)

Sets `x` to $\zeta(n)$ for even $n$ from the Bernoulli number $B_n$.
:::

## Bernoulli numbers and polynomials {#bernoulli}

::: function
void arb_bernoulli_ui(arb_t b, ulong n, slong prec)
void arb_bernoulli_fmpz(arb_t b, const fmpz_t n, slong prec)

Set `b` to the Bernoulli number $B_n$ to `prec` bits of accuracy.

::: details
The working precision is increased automatically. If $B_n$ is in the
global cache (see `bernoulli_cache_compute`) or its exact numerator has
about `prec` bits or more, the exact fraction is divided once; otherwise
$B_n$ is computed from $\zeta(n)$ by `arb_bernoulli_ui_zeta`. The
functions read the cache but do not extend it. A huge `fmpz` index gives a
result with a huge exponent, which is slower. The methods are described
in [](algorithms/number-theory/bernoulli).
:::
:::

::: function
void arb_bernoulli_ui_zeta(arb_t b, ulong n, slong prec)

Sets `b` to $B_n$ using $B_{2n} = (-1)^{n+1} 2 (2n)!\, \zeta(2n) /
(2\pi)^{2n}$, with $\zeta$ evaluated by the Euler product. This is
accurate only if `prec` is small enough compared to $n$ for the product
to converge quickly.
:::

::: function
void arb_bernoulli_poly_ui(arb_t res, ulong n, const arb_t x, slong prec)

Sets `res` to the Bernoulli polynomial $B_n(x)$. This is fast only if $n$
or $x$ is a small integer. Bernoulli numbers are read from the cache but
the cache is not extended.
:::

::: function
void arb_power_sum_vec(arb_ptr res, const arb_t a, const arb_t b, slong len, slong prec)

Sets `res[n]` to $S_n(a, b) = (B_{n+1}(b) - B_{n+1}(a))/(n+1)$ for $0 \le
n < \mathit{len}$, which equals $\sum_{k=a}^{b-1} k^n$ for integers $b
\ge a$. The values are computed together from the generating function of
the Bernoulli polynomials.
:::

## Polylogarithms {#polylog}

::: function
void arb_polylog(arb_t w, const arb_t s, const arb_t z, slong prec)
void arb_polylog_si(arb_t w, slong s, const arb_t z, slong prec)

Set `w` to the polylogarithm $\operatorname{Li}_s(z)$, computed with
`acb_polylog`. The result is indeterminate if `z` contains 1, and if the
midpoint of `z` is at least 1 unless $s$ is a nonpositive integer (where
the function is not real). In particular $\operatorname{Li}_s(1) =
\zeta(s)$ is not returned; use `acb_polylog` or `arb_zeta` for it.
:::

## Lambert W function {#lambertw}

::: function
void arb_lambertw(arb_t res, const arb_t x, int flags, slong prec)

Sets `res` to the Lambert W function, the solution $w$ of $w e^w = x$:
the principal branch $W_0$ if `flags` is 0, and the branch $W_{-1}$ if
`flags` is 1.

::: details
$W_0$ is real for $x \ge -1/e$, with values in $[-1, +\infty)$, and
$W_{-1}$ is real for $-1/e \le x < 0$, with values in $(-\infty, -1]$.
Elsewhere the function is complex; use `acb_lambertw`. The algorithm
computes a floating-point approximation heuristically and then certifies
an enclosure around it, with special code for asymptotic cases; it is
described in [@Joh2017b] and follows [@CGHJK1996]. See also
[](algorithms/special-functions/lambert-w).
:::
:::

## Other special functions {#other}

::: function
void arb_fib_fmpz(arb_t z, const fmpz_t n, slong prec)
void arb_fib_ui(arb_t z, ulong n, slong prec)

Set `z` to the Fibonacci number $F_n$, by binary squaring.
:::

::: function
void arb_agm(arb_t z, const arb_t x, const arb_t y, slong prec)

Sets `z` to the arithmetic-geometric mean of `x` and `y`. The result is
indeterminate if `x` or `y` contains negative numbers.
:::

::: function
void arb_chebyshev_t_ui(arb_t a, ulong n, const arb_t x, slong prec)
void arb_chebyshev_u_ui(arb_t a, ulong n, const arb_t x, slong prec)

Set `a` to the Chebyshev polynomial $T_n(x)$ of the first kind,
respectively $U_n(x)$ of the second kind.
:::

::: function
void arb_chebyshev_t2_ui(arb_t a, arb_t b, ulong n, const arb_t x, slong prec)
void arb_chebyshev_u2_ui(arb_t a, arb_t b, ulong n, const arb_t x, slong prec)

Set $a = T_n(x)$, $b = T_{n-1}(x)$, respectively $a = U_n(x)$, $b =
U_{n-1}(x)$. The outputs `a`, `b` and the input `x` must be distinct
objects.
:::

::: function
void arb_bell_fmpz(arb_t res, const fmpz_t n, slong prec)
void arb_bell_ui(arb_t res, ulong n, slong prec)

Set `res` to the Bell number $B_n$, exactly if it fits in `prec` bits and
otherwise as a numerical approximation computed efficiently by the
algorithm of [@Joh2015].
:::

::: function
void arb_bell_sum_bsplit(arb_t res, const fmpz_t n, const fmpz_t a, const fmpz_t b, const fmpz_t mmag, slong prec)
void arb_bell_sum_taylor(arb_t res, const fmpz_t n, const fmpz_t a, const fmpz_t b, const fmpz_t mmag, slong prec)

Set `res` to $\sum_{k=a}^{b-1} k^n / k!$, the sums used by `arb_bell_fmpz`.
If `mmag` is not `NULL`, the target error is $2^{\mathit{mmag} -
\mathit{prec}}$.
:::

::: function
void arb_euler_number_fmpz(arb_t res, const fmpz_t n, slong prec)
void arb_euler_number_ui(arb_t res, ulong n, slong prec)

Set `res` to the Euler number $E_n$, defined by the exponential
generating function $1/\cosh x$. The result is exact if $E_n$ fits in
`prec` bits.
:::

::: function
void arb_fmpz_euler_number_ui(fmpz_t res, ulong n)
void arb_fmpz_euler_number_ui_multi_mod(fmpz_t res, ulong n, double alpha)

Set `res` to the Euler number $E_n$ exactly. The default algorithm uses a
table, the Dirichlet beta function or a hybrid multimodular algorithm
depending on $n$ ([](algorithms/number-theory/bernoulli#euler)); the `multi_mod` version takes a tuning parameter
`alpha`, with a negative value selecting the default.
:::

::: function
void arb_partitions_fmpz(arb_t res, const fmpz_t n, slong prec)
void arb_partitions_ui(arb_t res, ulong n, slong prec)

Set `res` to the partition function $p(n)$. If $\log_2 p(n)$ exceeds
twice `prec`, the leading term of the Hardy-Ramanujan-Rademacher series
with an error bound is used; otherwise the exact value is computed and
rounded. See [](algorithms/number-theory/partitions).
:::

::: function
void arb_primorial_ui(arb_t res, ulong n, slong prec)
void arb_primorial_nth_ui(arb_t res, ulong n, slong prec)

Set `res` to the product of the primes up to $n$, respectively of the
first $n$ primes, in time quasilinear in $n$.
:::

## Random generation {#random}

::: function
void arb_randtest(arb_t x, flint_rand_t state, slong prec, slong mag_bits)
void arb_randtest_exact(arb_t x, flint_rand_t state, slong prec, slong mag_bits)
void arb_randtest_precise(arb_t x, flint_rand_t state, slong prec, slong mag_bits)
void arb_randtest_positive(arb_t x, flint_rand_t state, slong prec, slong mag_bits)
void arb_randtest_wide(arb_t x, flint_rand_t state, slong prec, slong mag_bits)
void arb_randtest_special(arb_t x, flint_rand_t state, slong prec, slong mag_bits)

Set `x` to a random ball for testing, with a midpoint of at most `prec`
bits and exponents of at most `mag_bits` bits: a finite ball, an exact
number, a ball with radius about $2^{-\mathit{prec}}$ times the
midpoint, such a ball that is certainly positive, a ball with independent
midpoint and radius (possibly very wide), and a ball that may have an
infinite or NaN midpoint or an infinite radius.
:::

::: function
void arb_urandom(arb_t x, flint_rand_t state, slong prec)

Sets `x` to a uniformly distributed random number in $[0, 1]$, rounded to
`prec` bits, so the radius may be nonzero.
:::

::: function
void arb_get_rand_fmpq(fmpq_t q, flint_rand_t state, const arb_t x, slong bits)

Sets `q` to a random rational number in the interval of `x`, whose
denominator is the binary denominator of `x` times a random integer of up
to `bits` bits. The midpoint and radius must be finite and of reasonable
size.
:::

## Vector functions {#vectors}

::: function
arb_ptr _arb_vec_entry_ptr(arb_ptr vec, slong i)

Returns a pointer to entry `i` of `vec`.
:::

::: function
void _arb_vec_zero(arb_ptr A, slong n)
void _arb_vec_indeterminate(arb_ptr vec, slong len)

Set all entries to zero, respectively to indeterminate balls.
:::

::: function
int _arb_vec_is_zero(arb_srcptr vec, slong len)
int _arb_vec_is_finite(arb_srcptr x, slong len)

Return whether all entries are exactly zero, respectively finite.
:::

::: function
int _arb_vec_equal(arb_srcptr vec1, arb_srcptr vec2, slong len)
int _arb_vec_overlaps(arb_srcptr vec1, arb_srcptr vec2, slong len)
int _arb_vec_contains(arb_srcptr vec1, arb_srcptr vec2, slong len)

Apply `arb_equal`, `arb_overlaps`, respectively `arb_contains`
entrywise and return whether it holds for all entries.
:::

::: function
void _arb_vec_set(arb_ptr res, arb_srcptr vec, slong len)
void _arb_vec_set_round(arb_ptr res, arb_srcptr vec, slong len, slong prec)
void _arb_vec_swap(arb_ptr res, arb_ptr vec, slong len)
void _arb_vec_neg(arb_ptr B, arb_srcptr A, slong n)
void _arb_vec_trim(arb_ptr res, arb_srcptr vec, slong len)

Copy, copy with rounding, swap, negate, and apply `arb_trim`
entrywise.
:::

::: function
void _arb_vec_add(arb_ptr C, arb_srcptr A, arb_srcptr B, slong n, slong prec)
void _arb_vec_sub(arb_ptr C, arb_srcptr A, arb_srcptr B, slong n, slong prec)

Set $C_i = A_i \pm B_i$.
:::

::: function
void _arb_vec_scalar_mul(arb_ptr res, arb_srcptr vec, slong len, const arb_t c, slong prec)
void _arb_vec_scalar_mul_arf(arb_ptr res, arb_srcptr vec, slong len, const arf_t c, slong prec)
void _arb_vec_scalar_mul_ui(arb_ptr res, arb_srcptr vec, slong len, ulong c, slong prec)
void _arb_vec_scalar_mul_si(arb_ptr res, arb_srcptr vec, slong len, slong c, slong prec)
void _arb_vec_scalar_mul_fmpz(arb_ptr res, arb_srcptr vec, slong len, const fmpz_t c, slong prec)
void _arb_vec_scalar_mul_2exp_si(arb_ptr res, arb_srcptr src, slong len, slong c)
void _arb_vec_scalar_div(arb_ptr res, arb_srcptr vec, slong len, const arb_t c, slong prec)
void _arb_vec_scalar_addmul(arb_ptr res, arb_srcptr vec, slong len, const arb_t c, slong prec)

Multiply each entry of `vec` by $c$ (by $2^c$ for the `_2exp_si`
version) or divide it by $c$, writing the result to `res`;
`_arb_vec_scalar_addmul` instead adds $c \cdot \mathit{vec}_i$ to
$\mathit{res}_i$.
:::

::: function
void _arb_vec_get_mag(mag_t bound, arb_srcptr vec, slong len)

Sets `bound` to an upper bound for the absolute values of all entries.
:::

::: function
slong _arb_vec_bits(arb_srcptr x, slong len)

Returns the maximum of `arb_bits` over the entries.
:::

::: function
void _arb_vec_set_powers(arb_ptr xs, const arb_t x, slong len, slong prec)

Sets `xs` to $1, x, x^2, \ldots, x^{\mathit{len}-1}$, computing even
powers by squaring and odd powers by one multiplication, which keeps the
error growth logarithmic in the exponent.
:::

::: function
void _arb_vec_add_error_arf_vec(arb_ptr res, arf_srcptr err, slong len)
void _arb_vec_add_error_mag_vec(arb_ptr res, mag_srcptr err, slong len)

Add $|\mathit{err}_i|$ to the radius of $\mathit{res}_i$.
:::

::: function
int _arb_vec_get_unique_fmpz_vec(fmpz * res, arb_srcptr vec, slong len)

Applies `arb_get_unique_fmpz` entrywise and returns nonzero if every
entry contains a unique integer; otherwise returns zero, with the
contents of `res` unspecified.
:::

::: function
void _arb_vec_printn(arb_srcptr vec, slong len, slong ndigits, ulong flags)
void _arb_vec_printd(arb_srcptr vec, slong len, slong ndigits)

Print the entries with `arb_printn`, respectively `arb_printd`.
:::

## Internals of the elementary functions {#internals}

These functions and tables implement the elementary functions and are
described in [](algorithms/special-functions/elementary). They use
fixed-point numbers: arrays of $n$ limbs representing a number in $[0,
1)$ (or with an extra integer limb where stated) with an error measured in
units of $2^{-64n}$ (ulps).

::: function
void _arb_exp_taylor_naive(nn_ptr y, ulong * error, nn_srcptr x, slong xn, ulong N)
void _arb_exp_taylor_rs(nn_ptr y, ulong * error, nn_srcptr x, slong xn, ulong N)

Set $y \approx \sum_{k=0}^{N-1} x^k/k!$ and write an error bound in ulps
to `error`; `y` has `xn` fractional limbs and one integer limb. The naive
version uses the forward recurrence, the `rs` version rectangular
splitting with a table of reciprocal factorials, avoiding divisions.
Require $N \le 287$, $0 \le x \le 1/16$ and $\mathit{xn} > 0$.
:::

::: function
void _arb_sin_cos_taylor_naive(nn_ptr ysin, nn_ptr ycos, ulong * error, nn_srcptr x, slong xn, ulong N)
void _arb_sin_cos_taylor_rs(nn_ptr ysin, nn_ptr ycos, ulong * error, nn_srcptr x, slong xn, ulong N, int sinonly, int alternating)

Set `ysin` and `ycos` to approximations of $\sum_{k=0}^{N-1} (-1)^k
x^{2k+1}/(2k+1)!$ and $\sum_{k=0}^{N-1} (-1)^k x^{2k}/(2k)!$, with `xn`
fractional limbs, and write an error bound in ulps to `error`. Require $N
\le 143$, $0 \le x \le 1/16$ and $\mathit{xn} > 0$.

::: details
In the `rs` version, `sinonly` = 1 computes only the sine, and
`alternating` = 0 computes the series without signs (the hyperbolic
sine, intended for use with `sinonly`, as in the exponential at high
precision, where $\cosh$ is then obtained from $\sqrt{1 + \sinh^2}$).
:::
:::

::: function
void _arb_atan_taylor_naive(nn_ptr y, ulong * error, nn_srcptr x, slong xn, ulong N, int alternating)
void _arb_atan_taylor_rs(nn_ptr y, ulong * error, nn_srcptr x, slong xn, ulong N, int alternating)

Set $y \approx \sum_{k=0}^{N-1} (\pm 1)^k x^{2k+1}/(2k+1)$, alternating
(the arctangent series) if `alternating` is 1 and not (the atanh series,
used for logarithms) if it is 0, and write an error bound in ulps to
`error`. Require $N \le 255$, $0 \le x \le 1/16$ and $\mathit{xn} > 0$.
:::

::: function
int _arb_get_mpn_fixed_mod_log2(nn_ptr w, fmpz_t q, ulong * error, const arf_t x, slong wn)

Attempts to write $x = q \log 2 + w$ with $0 \le w < \log 2$, `w` a
fixed-point number of `wn` limbs with error `error` ulps. Returns 0 if
the precision of the internal table of $\log 2$ is insufficient.
:::

::: function
int _arb_get_mpn_fixed_mod_pi4(nn_ptr w, fmpz_t q, int * octant, ulong * error, const arf_t x, slong wn)

Attempts to write $|x| = q \pi/4 + w$ with $0 \le w < \pi/4$, `w` a
fixed-point number of `wn` limbs with error `error` ulps, and writes $q
\bmod 8$ to `octant`. `q` may be `NULL`. Returns 0 if the precision of the
internal table of $\pi$ is insufficient.
:::

::: function
slong _arb_exp_taylor_bound(slong mag, slong prec)

Returns an $N$ such that $|\sum_{k=N}^{\infty} x^k/k!| \le
2^{-\mathit{prec}}$ for $|x| \le 2^{\mathit{mag}} \le 1/4$.
:::

::: function
flint_bitcnt_t _arb_mpn_leading_zeros(nn_srcptr d, slong n)

Returns the number of leading zero bits of the `n`-limb number `d`.
:::

::: function
void arb_exp_arf(arb_t z, const arf_t x, slong prec, int minus_one, slong maglim)

Sets `z` to $e^x$, or $e^x - 1$ if `minus_one` is set, for an exact `x`,
choosing among all algorithms; inputs with $|x| \ge 2^{\mathit{maglim}}$
give the crude overflow results of `arb_exp`.

::: details
Tiny $x$ uses $1 + x$ with error $x^2$. Below 2240 bits the argument is
reduced modulo $\log 2$ and by one or two table lookups in
`arb_exp_tab1`, respectively `arb_exp_tab21` and `arb_exp_tab22`, and a
Taylor series is evaluated with `_arb_exp_taylor_rs`. From 2240 bits to 4
million bits, `arb_exp_arf_log_reduction` is used, and above that
`arb_exp_arf_generic`.
:::
:::

::: function
void arb_exp_arf_generic(arb_t z, const arf_t x, slong prec, int minus_one)

Computes $e^x$ (or $e^x - 1$) without tables: reduces by $\log 2$ if
$|x|$ is large and then calls `arb_exp_arf_rs_generic` or
`arb_exp_arf_bb` depending on the precision, the size of $x$, the length
of its mantissa and the number of available threads.
:::

::: function
void arb_exp_arf_rs_generic(arb_t res, const arf_t x, slong prec, int minus_one)

Computes $e^x$ (or $e^x - 1$) by argument halving and rectangular
splitting of the Taylor series in ball arithmetic, for intermediate
precision.
:::

::: function
void arb_exp_arf_bb(arb_t z, const arf_t x, slong prec, int minus_one)

Computes $e^x$ (or $e^x - 1$, accurately) by the bit-burst algorithm:
$x$ is split into chunks of doubling length and the series for each chunk
is evaluated by binary splitting.

::: details
Raises an exception for extremely small or large $x$, which must be
handled by the caller. The argument is first halved until $|x| <
2^{-d}$ with $d \approx 16$ to speed up convergence, so $|x| \approx
2^m$ costs $m + d$ squarings at the end.
:::
:::

::: function
void arb_exp_arf_log_reduction(arb_t res, const arf_t x, slong prec, int minus_one)

Computes $e^x$ (or $e^x - 1$) by reduction with logarithms of small
primes [@HJ2024]: finds integers $c_i$ with $r = x - \sum c_i \log p_i$
small, evaluates $e^r$ and multiplies by $\prod p_i^{c_i}$.
:::

::: function
void _arb_exp_sum_bs_simple(fmpz_t T, fmpz_t Q, flint_bitcnt_t * Qexp, const fmpz_t x, flint_bitcnt_t r, slong N)
void _arb_exp_sum_bs_powtab(fmpz_t T, fmpz_t Q, flint_bitcnt_t * Qexp, const fmpz_t x, flint_bitcnt_t r, slong N)

Compute $T$, $Q$ and $Q_{\exp}$ with $T/(Q \cdot 2^{Q_{\exp}}) =
\sum_{k=1}^{N} (x/2^r)^k/k!$ by binary splitting (the sum runs to $N$
inclusive and omits the constant term). The `powtab` version precomputes
powers of $x$, which uses more memory and is faster; $N$ should have
many trailing zero bits.
:::

::: function
void arb_log_newton(arb_t res, const arb_t x, slong prec)
void arb_log_arf_newton(arb_t res, const arf_t x, slong prec)

Compute $\log x$ from a lower-precision value $t \approx \log x$ by the
higher-order iteration $\log x = t + \log(1 + (x e^{-t} - 1))$, with the
logarithm of the small correction evaluated by a short series. Used for
`arb_log_arf` above 2800 bits.
:::

::: function
void arb_atan_arf_bb(arb_t z, const arf_t x, slong prec)

Computes $\operatorname{atan} x$ by the bit-burst algorithm.

::: details
The argument-halving formula $\operatorname{atan}(x) = 2
\operatorname{atan}(x / (1 + \sqrt{1 + x^2}))$ is applied up to 8 times,
and then $\operatorname{atan}(x) = \operatorname{atan}(p/q) +
\operatorname{atan}(w)$ with $w = (qx - p)/(px + q)$, $p = \lfloor qx
\rfloor$, is applied repeatedly with growing $q = 2^k$, the arctangents of
the rational numbers $p/q$ being computed by binary splitting.
:::
:::

::: function
void _arb_atan_sum_bs_simple(fmpz_t T, fmpz_t Q, flint_bitcnt_t * Qexp, const fmpz_t x, flint_bitcnt_t r, slong N)
void _arb_atan_sum_bs_powtab(fmpz_t T, fmpz_t Q, flint_bitcnt_t * Qexp, const fmpz_t x, flint_bitcnt_t r, slong N)

Compute $T$, $Q$ and $Q_{\exp}$ with $T/(Q \cdot 2^{Q_{\exp}}) =
\sum_{k=1}^{N} (-1)^k (x/2^r)^{2k}/(2k+1)$ by binary splitting; the sum
must be multiplied by $x/2^r$ and the linear term added to give the
arctangent.
:::

::: function
void arb_atan_frac_bsplit(arb_t s, const fmpz_t p, const fmpz_t q, int hyperbolic, slong prec)

Sets `s` to $\operatorname{atan}(p/q)$, or $\operatorname{atanh}(p/q)$ if
`hyperbolic` is set, by binary splitting of the Taylor series.
:::

::: function
void arb_atan_newton(arb_t res, const arb_t x, slong prec)
void arb_atan_arf_newton(arb_t res, const arf_t x, slong prec)

Compute $\operatorname{atan} x$ from a lower-precision value $t$ by
$\operatorname{atan} x = t + \operatorname{atan}((x \cos t - \sin
t)/(\cos t + x \sin t))$, with the correction evaluated by a short
Taylor series, so that the precision multiplies by 5 to 11 per step.
Arguments with $|x| > 16$ use $\operatorname{atan} x = \pm\pi/2 -
\operatorname{atan}(1/x)$. Used for `arb_atan_arf` from 3400 bits.
:::

::: function
void arb_sin_cos_arf_generic(arb_t res_sin, arb_t res_cos, const arf_t x, slong prec)

Computes $\sin x$ and $\cos x$ without tables. If $|x|$ is not certainly
below $\pi/2$, the argument is reduced to the first octant using an
enclosure of $\pi$ and the function is called recursively; otherwise
`arb_sin_cos_arf_rs_generic` or `arb_sin_cos_arf_bb` is used depending on
the precision. Either output may be `NULL`.
:::

::: function
void arb_sin_cos_arf_rs_generic(arb_t res_sin, arb_t res_cos, const arf_t x, slong prec)

Computes $\sin x$ and $\cos x$ for $|x| \le 1$ by argument reduction
$x \to x/2^q$, rectangular splitting of the Taylor series and the double-angle
formula, in ball arithmetic, for $|x| < 1$. Either output may be
`NULL`.
:::

::: function
void arb_sin_cos_arf_bb(arb_t zsin, arb_t zcos, const arf_t x, slong prec)

Computes $\sin x$ and $\cos x$ by the bit-burst algorithm. Requires $|x|
< \pi/2$, which is not checked.
:::

::: function
void arb_sin_cos_arf_atan_reduction(arb_t res1, arb_t res2, const arf_t x, slong prec)

Computes $\sin x$ and $\cos x$ by reduction with the arguments of small
Gaussian primes [@HJ2024]: $x - \sum c_i \theta_i$ is made small with
$\theta_i = \arg(a_i + b_i i)$, and the result is multiplied by the
Gaussian integer $\prod (a_i + b_i i)^{c_i}$ normalised to modulus 1. Used
for `arb_sin_cos` from 2600 bits to 4 million bits.
:::

::: function
void arb_sin_cos_wide(arb_t s, arb_t c, const arb_t x, slong prec)
void _arb_sin_cos_wide(arb_t s, arb_t c, const arf_t x, const mag_t r, slong prec)

Compute enclosures of the ranges of $\sin$ and $\cos$ on a wide ball $x =
[m \pm r]$, with both endpoints optimal to within about $2^{-30}$.

::: details
The sine and cosine are evaluated in `double` arithmetic at the endpoints,
and the extrema $\pm 1$ are included when a peak lies between them; the
result is converted back to a ball. The accumulated error of the simple
floating-point algorithm is below $2^{-30}$, the dominant part coming from
a single approximate reduction by $\pi/2$, which is safe in `double`
arithmetic up to about $2^{20}$; larger arguments are first reduced in
`arb_t` arithmetic.
:::
:::

::: function
void arb_sin_cos_generic(arb_t s, arb_t c, const arb_t x, slong prec)
void _arb_sin_cos_generic(arb_t s, arb_t c, const arf_t x, const mag_t xrad, slong prec)

Compute $\sin$ and $\cos$ of the ball $[x \pm \mathit{xrad}]$, handling
special values and the propagated error and calling
`arb_sin_cos_arf_generic`. This is the fallback of `arb_sin_cos` for the
cases without a fast path.
:::

::: function
void _arb_sin_pi_fmpq_algebraic(arb_t s, ulong p, ulong q, slong prec)
void _arb_cos_pi_fmpq_algebraic(arb_t c, ulong p, ulong q, slong prec)
void _arb_sin_cos_pi_fmpq_algebraic(arb_t s, arb_t c, ulong p, ulong q, slong prec)

Compute $\sin(\pi p/q)$ and $\cos(\pi p/q)$ as algebraic numbers: an
isolated root of the minimal polynomial of the cosine
(`fmpz_poly_cos_minpoly`) is refined by Newton iteration. This is used for
small denominators, where it is faster than evaluating the functions.
:::

::: function
void arb_log_primes_vec_bsplit(arb_ptr res, slong n, slong prec)

Sets `res` to $\log p$ for the first `n` primes, by binary splitting
applied to simultaneous Machin-type formulas. Not optimised for large `n`
or small `prec`.
:::

::: function
void arb_atan_gauss_primes_vec_bsplit(arb_ptr res, slong n, slong prec)

Sets `res` to the arguments of the first `n` Gaussian primes $a + bi$
with $a > b > 0$ (one per rational prime $p \equiv 1 \bmod 4$, and $1 +
i$), by binary splitting applied to simultaneous Machin-type formulas.
:::

::: function
void _arb_log_p_ensure_cached(slong prec)
arb_srcptr _arb_log_p_cache_vec(void)
void _arb_atan_gauss_p_ensure_cached(slong prec)
arb_srcptr _arb_atan_gauss_p_cache_vec(void)

Make sure that the thread-local caches of the logarithms of the first
`ARB_LOG_PRIME_CACHE_NUM` primes, respectively of the arguments of the
first `ARB_ATAN_GAUSS_PRIME_CACHE_NUM` Gaussian primes, hold at least
`prec` bits, and return the cached vectors.
:::

::: function
void _arb_log_precompute_reductions(short * rel, double * eps, arb_srcptr alpha, slong n, slong max_rel, double C)

Precomputes, by lattice reduction, short integer relations
$\sum_j c_j \alpha_j = \varepsilon$ with small $|\varepsilon|$ between
the constants $\alpha_0, \ldots, \alpha_{n-1}$ (logarithms of primes or
arguments of Gaussian primes), as used by the reduction algorithm of
[@HJ2024]: `rel` receives the coefficient vectors and `eps` the values
$\varepsilon$, for successively smaller $\varepsilon$ decreasing roughly
by the factor `C`, up to `max_rel` relations. The tables actually used
are precomputed in `mp_real`.
:::

::: macro
ARB_LOG_PRIME_CACHE_NUM
ARB_ATAN_GAUSS_PRIME_CACHE_NUM
ARB_LOG_REDUCTION_DEFAULT_MAX_PREC
ARB_SIN_COS_ATAN_REDUCTION_DEFAULT_MAX_PREC
ARB_EXP_LOG_REDUCTION_PREC
ARB_SIN_COS_ATAN_REDUCTION_PREC
ARB_LOG_NEWTON_PREC
ARB_ATAN_NEWTON_PREC

The number of cached logarithms of primes and arguments of Gaussian
primes (13 each), the precision range in which the reduction algorithms
are used for the exponential and for sine and cosine (from 2240,
respectively 2600 bits, to 4 million bits), and the precisions from which
the logarithm (2800 bits) and the arctangent (3400 bits) use Newton
iteration.
:::

::: variable
arb_exp_tab1
arb_exp_tab21
arb_exp_tab22
arb_sin_cos_tab1
arb_sin_cos_tab21
arb_sin_cos_tab22
arb_log_tab11
arb_log_tab12
arb_log_tab21
arb_log_tab22
arb_atan_tab1
arb_atan_tab21
arb_atan_tab22

Precomputed fixed-point tables for the table-based argument reduction.
The `tab1` tables (512 bits) serve one-step reduction with 8 index bits;
the `tab21` and `tab22` tables (4608 bits) serve two-step reduction with
5 + 5 index bits. They hold $e^{p/2^8}$, $\sin$ and $\cos$ of $p/2^8$,
$\log(1 + p/2^7)$ and $\operatorname{atan}(p/2^8)$ and the analogous
values for the two-step tables.
:::

::: variable
arb_pi4_tab
arb_log_log2_tab
arb_atan_pi2_minus_one
arb_log_p_tab
arb_atan_gauss_tab

Fixed-point values of $\pi/4$, $\log 2$ and $\pi/2 - 1$ to 4608 bits, and
of the logarithms of the first 13 primes and the arguments of the first
13 Gaussian primes, used by the argument reductions and as the initial
cache of the constants.
:::

## Other internal functions {#internal}

::: function
void _arb_digits_round_inplace(char * s, flint_bitcnt_t * shift, fmpz_t error, slong n, arf_rnd_t rnd)

Rounds the decimal digit string `s` to `n` digits in the direction `rnd`,
in place, writing the number of removed digits to `shift` and the
rounding error in units of the last kept digit to `error`. Used by the
decimal conversion.
:::
