---
title: mag.h
header: mag.h
description: Unsigned floating-point numbers with a 30-bit mantissa for fast upper and lower bounds
category: Real and complex numbers
order: 1
---

A `mag_t` holds a nonnegative real number of the form $m \cdot 2^{e - 30}$
with a 30-bit mantissa $m$ and an arbitrary-precision exponent $e$, or
positive infinity. It is the type of the radius of an `arb_t` ball and,
more generally, the type FLINT uses for error bounds and other magnitudes
that only need to be known to a few digits. Every rounding operation has
a direction fixed by its name: the plain functions return an upper bound,
the functions with suffix `_lower` a lower bound. No attempt is made to
compute the best possible bound; results can be off by a few units in the
last place, always in the safe direction.

Because the mantissa fits in a word and the exponent is almost always a
small `fmpz` stored inline, `mag_t` arithmetic is several times faster
than `arf_t` arithmetic at any precision, which is why ball arithmetic
pays almost nothing for its radius. The algorithm page
[](algorithms/ball-arithmetic/representation) describes the format and
the rounding rules in detail; the guide [](mag/bounds) explains how to
combine upper and lower bounds correctly.

For floating-point arithmetic with more precision, correct rounding or
negative numbers, use `arf_t`.

## Representation {#representation}

::: type
mag_struct
mag_t
mag_ptr
mag_srcptr

A `mag_struct` has two fields, an `fmpz` exponent `exp` and a `ulong`
mantissa `man`. A finite nonzero value is $\mathit{man} \cdot
2^{\mathit{exp} - 30}$ with $2^{29} \le \mathit{man} < 2^{30}$, that is,
a number in $[2^{\mathit{exp}-1}, 2^{\mathit{exp}})$. Zero is stored as
$\mathit{man} = 0$, $\mathit{exp} = 0$, and $+\infty$ as $\mathit{man} = 0$
with a nonzero exponent. There is no NaN. A `mag_t` is an array of one
`mag_struct`, so that it is passed by reference; `mag_ptr` and
`mag_srcptr` point to arrays of `mag_struct`.
:::

::: macro
MAG_BITS
MAG_EXP(x)
MAG_EXPREF(x)
MAG_MAN(x)

`MAG_BITS` is the number of mantissa bits, 30. `MAG_EXP(x)` and
`MAG_MAN(x)` access the exponent (as an `fmpz` value) and the mantissa of
`x`, and `MAG_EXPREF(x)` is a pointer to the exponent, for use with `fmpz`
functions.
:::

Conventions for all functions below:

- Every function allows aliasing between its inputs and its output.
- A function named `mag_f` sets the output to an upper bound for the
  exact result; `mag_f_lower` sets it to a lower bound. Functions that are
  always exact, such as `mag_set` and `mag_mul_2exp_si`, have no lower
  version.
- Infinity is handled as the limit value: an upper bound may always be
  $+\infty$, and lower bounds of quantities that are not bounded below by a
  positive number are 0. Operations with no meaningful value, such as
  $0 \cdot \infty$, give the safe answer for the direction ($+\infty$ for
  upper bounds, 0 for lower bounds) unless stated otherwise.
- Only the `mag_fast_` functions require exponents to be small; all other
  functions work for exponents of any size.

## Memory management {#memory}

::: function
void mag_init(mag_t x)
void mag_init_set(mag_t x, const mag_t y)

Initialise `x` to zero, respectively to a copy of `y`.
:::

::: function
void mag_clear(mag_t x)

Frees the memory used by `x`. Only an exponent too large for an inline
`fmpz` uses heap memory.
:::

::: function
void mag_swap(mag_t x, mag_t y)

Swaps `x` and `y` without copying.
:::

::: function
mag_ptr _mag_vec_init(slong n)
void _mag_vec_clear(mag_ptr v, slong n)

Allocate a vector of `n` entries initialised to zero, and free such a
vector.
:::

::: function
slong mag_allocated_bytes(const mag_t x)

Returns the number of bytes allocated on the heap by `x`, excluding
`sizeof(mag_struct)`.
:::

## Special values {#special}

::: function
void mag_zero(mag_t res)
void mag_one(mag_t res)
void mag_inf(mag_t res)

Set `res` to 0, 1 and $+\infty$ respectively.
:::

::: function
int mag_is_special(const mag_t x)
int mag_is_zero(const mag_t x)
int mag_is_inf(const mag_t x)
int mag_is_finite(const mag_t x)

Return whether `x` is zero or infinite, zero, infinite, and finite
respectively. Since there is no NaN, `mag_is_finite` is the negation of
`mag_is_inf`.
:::

## Assignment and conversion {#conversion}

::: function
void mag_set(mag_t res, const mag_t x)

Sets `res` to `x`. Exact.
:::

::: function
void mag_set_d(mag_t res, double x)
void mag_set_ui(mag_t res, ulong x)
void mag_set_fmpz(mag_t res, const fmpz_t x)

Set `res` to an upper bound for $|x|$. The result can be larger than
$|x|$ even when $|x|$ is representable: `mag_set_d` always adds one unit
in the 30th bit. A NaN `double` gives $+\infty$.
:::

::: function
void mag_set_d_lower(mag_t res, double x)
void mag_set_ui_lower(mag_t res, ulong x)
void mag_set_fmpz_lower(mag_t res, const fmpz_t x)

Set `res` to a lower bound for $|x|$. A NaN `double` gives 0.
:::

::: function
void mag_set_d_2exp_fmpz(mag_t res, double x, const fmpz_t y)
void mag_set_fmpz_2exp_fmpz(mag_t res, const fmpz_t x, const fmpz_t y)
void mag_set_ui_2exp_si(mag_t res, ulong x, slong y)

Set `res` to an upper bound for $|x| \cdot 2^y$.
:::

::: function
void mag_set_d_2exp_fmpz_lower(mag_t res, double x, const fmpz_t y)
void mag_set_fmpz_2exp_fmpz_lower(mag_t res, const fmpz_t x, const fmpz_t y)

Set `res` to a lower bound for $|x| \cdot 2^y$.
:::

::: function
double mag_get_d(const mag_t x)

Returns a `double` upper bound for `x`. The result is exact when the
exponent of `x` lies in $[-1000, 1000]$; outside this range it is
$2^{-1000}$ for small `x` and $+\infty$ for large `x`.
:::

::: function
double mag_get_d_log2_approx(const mag_t x)

Returns an approximation of $\log_2 x$, clamped to
$[\mathrm{COEFF\_MIN}, \mathrm{COEFF\_MAX}]$, for estimating magnitudes.
It is not a rigorous bound.
:::

::: function
void mag_get_fmpq(fmpq_t res, const mag_t x)
void mag_get_fmpz(fmpz_t res, const mag_t x)
void mag_get_fmpz_lower(fmpz_t res, const mag_t x)

Set `res` to the exact value of `x`, to $\lceil x \rceil$, and to $\lfloor
x \rfloor$ respectively.

::: details
The caller must make sure that `x` has a reasonable size. An infinite `x`
or one with a huge exponent raises an exception, and an exponent that is
merely large can exhaust memory.
:::
:::

## Comparison {#comparison}

::: function
int mag_equal(const mag_t x, const mag_t y)

Returns whether `x` and `y` are equal.
:::

::: function
int mag_cmp(const mag_t x, const mag_t y)
int mag_cmp_2exp_si(const mag_t x, slong e)

Return a negative number, zero or a positive number according as `x` is
less than, equal to or greater than `y`, respectively $2^e$. Both
comparisons are exact and handle infinity.
:::

::: function
void mag_min(mag_t res, const mag_t x, const mag_t y)
void mag_max(mag_t res, const mag_t x, const mag_t y)

Set `res` to the smaller, respectively the larger, of `x` and `y`. Exact.
:::

## Arithmetic {#arithmetic}

::: function
void mag_add(mag_t res, const mag_t x, const mag_t y)
void mag_add_ui(mag_t res, const mag_t x, ulong y)

Set `res` to an upper bound for $x + y$.
:::

::: function
void mag_add_lower(mag_t res, const mag_t x, const mag_t y)
void mag_add_ui_lower(mag_t res, const mag_t x, ulong y)

Set `res` to a lower bound for $x + y$.
:::

::: function
void mag_add_2exp_fmpz(mag_t res, const mag_t x, const fmpz_t e)
void mag_add_ui_2exp_si(mag_t res, const mag_t x, ulong y, slong e)

Set `res` to an upper bound for $x + 2^e$, respectively $x + y 2^e$.
:::

::: function
void mag_sub(mag_t res, const mag_t x, const mag_t y)
void mag_sub_lower(mag_t res, const mag_t x, const mag_t y)

Set `res` to an upper, respectively lower, bound for $\max(x - y, 0)$.

::: details
To bound $x - y$ from above, `x` should itself be an upper bound and `y`
a lower bound; `mag_sub_lower` is typically used to bound a denominator
such as $1 - |z|$ away from zero.
:::
:::

::: function
void mag_mul_2exp_si(mag_t res, const mag_t x, slong y)
void mag_mul_2exp_fmpz(mag_t res, const mag_t x, const fmpz_t y)

Set `res` to $x \cdot 2^y$. Exact.
:::

::: function
void mag_mul(mag_t res, const mag_t x, const mag_t y)
void mag_mul_ui(mag_t res, const mag_t x, ulong y)
void mag_mul_fmpz(mag_t res, const mag_t x, const fmpz_t y)

Set `res` to an upper bound for $xy$ (with $|y|$ for the integer
versions). $0 \cdot \infty$ gives $+\infty$.
:::

::: function
void mag_mul_lower(mag_t res, const mag_t x, const mag_t y)
void mag_mul_ui_lower(mag_t res, const mag_t x, ulong y)
void mag_mul_fmpz_lower(mag_t res, const mag_t x, const fmpz_t y)

Set `res` to a lower bound for $xy$. $0 \cdot \infty$ gives 0.
:::

::: function
void mag_addmul(mag_t z, const mag_t x, const mag_t y)

Sets `z` to an upper bound for $z + xy$.
:::

::: function
void mag_div(mag_t res, const mag_t x, const mag_t y)
void mag_div_ui(mag_t res, const mag_t x, ulong y)
void mag_div_fmpz(mag_t res, const mag_t x, const fmpz_t y)

Set `res` to an upper bound for $x / y$. Division by zero, and
$\infty/\infty$, give $+\infty$.
:::

::: function
void mag_div_lower(mag_t res, const mag_t x, const mag_t y)

Sets `res` to a lower bound for $x / y$. Division by zero, and
$\infty/\infty$, give 0.
:::

::: function
void mag_inv(mag_t res, const mag_t x)
void mag_inv_lower(mag_t res, const mag_t x)

Set `res` to an upper, respectively lower, bound for $1 / x$. The upper
bound for $1/0$ is $+\infty$; the lower bound for $1/0$ is 0, which is
valid but weak.
:::

## Powers, roots, exponentials and logarithms {#elementary}

::: function
void mag_pow_ui(mag_t res, const mag_t x, ulong e)
void mag_pow_fmpz(mag_t res, const mag_t x, const fmpz_t e)

Set `res` to an upper bound for $x^e$. A negative exponent in
`mag_pow_fmpz` is handled as $(1/x)^{-e}$.
:::

::: function
void mag_pow_ui_lower(mag_t res, const mag_t x, ulong e)
void mag_pow_fmpz_lower(mag_t res, const mag_t x, const fmpz_t e)

Set `res` to a lower bound for $x^e$.
:::

::: function
void mag_sqrt(mag_t res, const mag_t x)
void mag_sqrt_lower(mag_t res, const mag_t x)

Set `res` to an upper, respectively lower, bound for $\sqrt{x}$.
:::

::: function
void mag_rsqrt(mag_t res, const mag_t x)
void mag_rsqrt_lower(mag_t res, const mag_t x)

Set `res` to an upper, respectively lower, bound for $1/\sqrt{x}$, with
$1/\sqrt{0} = +\infty$.
:::

::: function
void mag_hypot(mag_t res, const mag_t x, const mag_t y)

Sets `res` to an upper bound for $\sqrt{x^2 + y^2}$.
:::

::: function
void mag_root(mag_t res, const mag_t x, ulong n)

Sets `res` to an upper bound for $x^{1/n}$. For $n = 0$ the result is
$+\infty$.
:::

::: function
void mag_log(mag_t res, const mag_t x)
void mag_log_lower(mag_t res, const mag_t x)

Set `res` to an upper, respectively lower, bound for $\log(\max(1, x))$.
:::

::: function
void mag_neg_log(mag_t res, const mag_t x)
void mag_neg_log_lower(mag_t res, const mag_t x)

Set `res` to an upper, respectively lower, bound for $-\log(\min(1, x))$,
which is $|\log x|$ for $x \le 1$. The value for $x = 0$ is $+\infty$.
:::

::: function
void mag_log_ui(mag_t res, ulong n)

Sets `res` to an upper bound for $\log n$. For $n = 0$ the result is
$+\infty$.
:::

::: function
void mag_log1p(mag_t res, const mag_t x)

Sets `res` to an upper bound for $\log(1 + x)$, with small relative error
also when `x` is small.
:::

::: function
void mag_exp(mag_t res, const mag_t x)
void mag_exp_lower(mag_t res, const mag_t x)

Set `res` to an upper, respectively lower, bound for $e^x$.
:::

::: function
void mag_expinv(mag_t res, const mag_t x)
void mag_expinv_lower(mag_t res, const mag_t x)

Set `res` to an upper, respectively lower, bound for $e^{-x}$.
:::

::: function
void mag_expm1(mag_t res, const mag_t x)

Sets `res` to an upper bound for $e^x - 1$, with small relative error also
when `x` is small.
:::

::: function
void mag_cosh(mag_t res, const mag_t x)
void mag_cosh_lower(mag_t res, const mag_t x)
void mag_sinh(mag_t res, const mag_t x)
void mag_sinh_lower(mag_t res, const mag_t x)

Set `res` to an upper or lower bound for $\cosh x$ or $\sinh x$.
:::

::: function
void mag_atan(mag_t res, const mag_t x)
void mag_atan_lower(mag_t res, const mag_t x)

Set `res` to an upper, respectively lower, bound for $\operatorname{atan}
x$.
:::

::: function
void mag_const_pi(mag_t res)
void mag_const_pi_lower(mag_t res)

Set `res` to an upper, respectively lower, bound for $\pi$.
:::

::: function
double mag_d_log_upper_bound(double x)
double mag_d_log_lower_bound(double x)

Return an upper, respectively lower, bound for $\log x$ as a `double`,
for $x > 0$ (with $\log 0 = -\infty$).

::: details
These are the building blocks of the logarithm and exponential bounds.
They do not call the `libm` logarithm, whose error is not specified by
any standard, but evaluate $\log x = n \log 2 + \log(m/32) + \log(1 +
v)$ with a 32-entry table and a short series, which has an error of at
most 10 units in the last place, and then multiply by $1 \pm 10^{-14}$.
The bounds are valid even if the floating-point operations are not
correctly rounded, as long as each has error at most about one unit in the
last place.
:::
:::

## Bounds for factorials, series and special functions {#series}

::: function
void mag_fac_ui(mag_t res, ulong n)
void mag_rfac_ui(mag_t res, ulong n)

Set `res` to an upper bound for $n!$, respectively $1/n!$.
:::

::: function
void mag_bin_uiui(mag_t res, ulong n, ulong k)

Sets `res` to an upper bound for the binomial coefficient $\binom{n}{k}$.
:::

::: function
void mag_binpow_uiui(mag_t res, ulong m, ulong n)

Sets `res` to an upper bound for $(1 + 1/m)^n$.
:::

::: function
void mag_bernoulli_div_fac_ui(mag_t res, ulong n)

Sets `res` to an upper bound for $|B_n| / n!$, where $B_n$ is a Bernoulli
number.
:::

::: function
void mag_geom_series(mag_t res, const mag_t x, ulong N)

Sets `res` to an upper bound for $\sum_{k=N}^{\infty} x^k = x^N / (1 -
x)$, which is $+\infty$ if $x \ge 1$.
:::

::: function
void mag_exp_tail(mag_t res, const mag_t x, ulong N)

Sets `res` to an upper bound for $\sum_{k=N}^{\infty} x^k / k!$.

::: details
If $N \ge 2x$, the terms decrease at least geometrically with ratio
$1/2$ and the bound $2 x^N / N!$ is used. Otherwise the crude bound $e^x$
is returned.
:::
:::

::: function
void mag_polylog_tail(mag_t res, const mag_t z, slong s, ulong d, ulong N)

Sets `res` to an upper bound for

$$
\sum_{k=N}^{\infty} \frac{z^k \log^d(k)}{k^s}.
$$

::: details
This is the tail bound used when summing polylogarithms and their
derivatives. When $s$ is real or complex, any integer $s' \le
\operatorname{Re}(s)$ can be substituted. The derivation is described in
[](algorithms/special-functions/zeta).
:::
:::

::: function
void mag_hurwitz_zeta_uiui(mag_t res, ulong s, ulong a)

Sets `res` to an upper bound for $\zeta(s, a) = \sum_{k=0}^{\infty}
(k+a)^{-s}$.

::: details
The bound is

$$
\zeta(s,a) \le \frac{1}{a^s} + \frac{1}{(s-1) a^{s-1}},
$$

obtained by comparing the sum with an integral. It is $+\infty$ if $s \le
1$ or $a = 0$.
:::
:::

## Input and output {#io}

::: function
void mag_print(const mag_t x)
void mag_fprint(FILE * file, const mag_t x)

Print the exact value of `x` as `(man * 2^exp)`, with the mantissa as an
integer, to standard output or to `file`. Zero and infinity are printed as
`(0)` and `(inf)`.
:::

::: function
void mag_printd(const mag_t x, slong d)
void mag_fprintd(FILE * file, const mag_t x, slong d)

Print `x` in decimal with `d` significant digits, to standard output or
to `file`. The decimal conversion is not rigorous.
:::

::: function
char * mag_dump_str(const mag_t x)
int mag_load_str(mag_t x, const char * str)

`mag_dump_str` returns a string, to be freed with `flint_free`, holding a
binary representation of `x` that `mag_load_str` reads back exactly.
`mag_load_str` returns nonzero if `str` is not in this format.
:::

::: function
int mag_dump_file(FILE * stream, const mag_t x)
int mag_load_file(mag_t x, FILE * stream)

Write `x` to `stream` in the format of `mag_dump_str`, and read it back.
Both return nonzero on failure. The data is delimited by whitespace or end
of file, so several values written in sequence must be separated by
whitespace.
:::

## Random generation {#random}

::: function
void mag_randtest(mag_t res, flint_rand_t state, slong expbits)
void mag_randtest_special(mag_t res, flint_rand_t state, slong expbits)

Set `res` to a random finite value with an exponent of up to `expbits`
bits. `mag_randtest_special` also produces 0 and $+\infty$ sometimes.
:::

## Fast arithmetic with small exponents {#fast}

The following functions are used in inner loops of ball arithmetic, for
example in `arb_dot`. They assume that all inputs are finite and that all
exponents, in the inputs and in the result, fit in an inline `fmpz`
(below $2^{60}$ in absolute value on 64-bit systems). They overwrite the
exponent of the output directly, so the output must not hold a heap
allocated exponent, which would be leaked. Nothing is checked.

::: function
void mag_fast_init_set(mag_t x, const mag_t y)
void mag_fast_zero(mag_t res)
int mag_fast_is_zero(const mag_t x)

Initialise `x` to a copy of `y`, set `res` to zero, and test whether `x`
is zero.
:::

::: function
void mag_fast_mul(mag_t res, const mag_t x, const mag_t y)
void mag_fast_addmul(mag_t z, const mag_t x, const mag_t y)
void mag_fast_add_2exp_si(mag_t res, const mag_t x, slong e)

Set `res` to an upper bound for $xy$, `z` to an upper bound for $z + xy$,
and `res` to an upper bound for $x + 2^e$.
:::

::: function
void mag_fast_mul_2exp_si(mag_t res, const mag_t x, slong e)

Sets `res` to $x \cdot 2^e$. Exact.
:::
