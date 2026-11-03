---
title: arf.h
header: arf.h
description: Arbitrary-precision binary floating-point numbers with correct rounding and unbounded exponents
category: Real and complex numbers
order: 2
---

An `arf_t` holds a binary floating-point number $m \cdot 2^e$ with an
integer mantissa $m$ of any length and an exponent $e$ of any size, or one
of the special values $0$, $+\infty$, $-\infty$ and NaN. It is the type of
the midpoint of an `arb_t` ball, and it can be used on its own as a
correctly rounded floating-point type in the style of MPFR. Unlike MPFR,
an `arf_t` has no fixed precision: the precision is a parameter of each
operation, the mantissa grows and shrinks as needed, and numbers whose
mantissa is short (an integer, a power of two, a `double`) stay short.
Exponents are `fmpz` integers, so overflow and underflow cannot occur.
There is no negative zero, no unsigned infinity and no NaN payload.

The guide [](arf/rounding) explains the rounding modes, the inexact flag
and the precision argument with examples. The representation is described
in [](algorithms/ball-arithmetic/representation).

## Representation and conventions {#conventions}

::: type
arf_struct
arf_t
arf_ptr
arf_srcptr

An `arf_struct` has four words: an `fmpz` exponent, a `size` field
holding the number of limbs of the mantissa and, in its lowest bit, the
sign, and two words that store a mantissa of at most two limbs inline or
otherwise the allocation size and a pointer to the limbs. Numbers with at
most 128 bits (64 bits on a 32-bit system) therefore use no heap memory.
A nonzero finite value is stored as $\pm 0.d_{n-1} \ldots d_0 \cdot
2^{\mathit{exp}}$ with $n$ limbs, the top bit of $d_{n-1}$ set and $d_0
\ne 0$, so the stored exponent `exp` is the $e$ with $2^{e-1} \le |x| <
2^e$. A special value has `size` zero and the exponent encodes which one
it is. `arf_t` is an array of one `arf_struct`; `arf_ptr` and
`arf_srcptr` point to arrays.
:::

::: type
arf_rnd_t

A rounding mode: one of `ARF_RND_DOWN`, `ARF_RND_UP`, `ARF_RND_FLOOR`,
`ARF_RND_CEIL`, `ARF_RND_NEAR`, or one of the relaxed modes
`ARF_RND_FAST` and `ARF_RND_ACCURATE`.
:::

::: macro
ARF_RND_DOWN
ARF_RND_UP
ARF_RND_FLOOR
ARF_RND_CEIL
ARF_RND_NEAR

Round towards zero, away from zero, towards $-\infty$, towards $+\infty$,
and to the nearest representable number with ties to even.
:::

::: macro
ARF_RND_FAST
ARF_RND_ACCURATE

Relaxed modes, accepted only by division, square root and reciprocal
square root (`arf_div`, `arf_sqrt`, `arf_rsqrt` and the functions built on
them). The result differs from the exact value by at most one unit in the
last place in either direction (`ARF_RND_FAST`), respectively by at most
0.51 units (`ARF_RND_ACCURATE`). The return value is the inexact flag, as
with the strict modes: it is 0 if and only if the result equals the exact
value.

::: details
A unit in the last place here is $2^{e - \mathit{prec}}$ where $2^{e-1}
\le |f| < 2^e$ for the exact result $f$. The functions treat the
relaxed modes as the strict modes of `arf_rnd_relaxed_to_strict`
(`ARF_RND_DOWN` for `ARF_RND_FAST`, `ARF_RND_NEAR` for
`ARF_RND_ACCURATE`), which the limb-array kernels (see [](#mpn))
perform correctly at every precision, so the inexact flag is as
reliable as with the strict modes.
:::
:::

::: macro
ARF_PREC_EXACT

A precision (`WORD_MAX`) meaning that no rounding takes place. It may only
be passed to operations whose exact result is a dyadic number, such as
addition, multiplication and conversion from integers, and only when that
result fits in memory: adding $1$ and $2^{-10^{12}}$ exactly needs
$10^{12}$ bits. If in doubt, use a large finite precision and check the
inexact flag instead.
:::

Unless stated otherwise:

- A function with `prec` and `rnd` arguments that returns `int` rounds its
  result to at most `prec` bits of mantissa in the direction `rnd`, and
  the rounding is *correct*: the result is the exact result of the
  operation on the exact inputs, rounded once. The return value is 0 if
  the result is exact and 1 if it was rounded (the *inexact flag*). The
  error is less than one unit in the last place for directed rounding and
  at most half a unit for `ARF_RND_NEAR`.
- A function without a precision argument is exact.
- The precision must be at least 2 bits.
- Every function allows aliasing between its inputs and its outputs.
- Special values follow IEEE 754 conventions: $0 \cdot \infty$, $\infty
  - \infty$, $0/0$ and $\sqrt{-1}$ give NaN, and so does division of a
  nonzero number by zero, since there is no signed zero to choose the sign
  of the infinity.

## Memory management {#memory}

::: function
void arf_init(arf_t x)

Initialises `x` to zero.
:::

::: function
void arf_clear(arf_t x)

Frees the memory used by `x`.
:::

::: function
arf_ptr _arf_vec_init(slong n)
void _arf_vec_clear(arf_ptr v, slong n)

Allocate a vector of `n` entries initialised to zero, and free such a
vector.
:::

::: function
slong arf_allocated_bytes(const arf_t x)

Returns the number of bytes allocated on the heap by `x`, excluding
`sizeof(arf_struct)`.
:::

## Special values {#special}

::: function
void arf_zero(arf_t x)
void arf_one(arf_t x)
void arf_pos_inf(arf_t x)
void arf_neg_inf(arf_t x)
void arf_nan(arf_t x)

Set `x` to $0$, $1$, $+\infty$, $-\infty$ and NaN respectively.
:::

::: function
int arf_is_zero(const arf_t x)
int arf_is_one(const arf_t x)
int arf_is_pos_inf(const arf_t x)
int arf_is_neg_inf(const arf_t x)
int arf_is_nan(const arf_t x)
int arf_is_inf(const arf_t x)

Return whether `x` is $0$, $1$, $+\infty$, $-\infty$, NaN, and
$\pm\infty$ respectively.
:::

::: function
int arf_is_normal(const arf_t x)
int arf_is_special(const arf_t x)
int arf_is_finite(const arf_t x)

Return whether `x` is a finite nonzero number, whether it is one of the
special values $0$, $\pm\infty$, NaN, and whether it is finite (zero
included). `arf_is_finite` is not the negation of `arf_is_inf`,
because of NaN.
:::

## Assignment and rounding {#assignment}

::: function
void arf_set(arf_t y, const arf_t x)
void arf_set_ui(arf_t x, ulong v)
void arf_set_si(arf_t x, slong v)
void arf_set_mpz(arf_t y, const mpz_t x)
void arf_set_fmpz(arf_t y, const fmpz_t x)
void arf_set_d(arf_t x, double v)
void arf_set_mpfr(arf_t x, const mpfr_t y)

Set the output to the exact value of the input. A `double` NaN gives NaN
and $-0.0$ gives $0$.
:::

::: function
void arf_init_set_ui(arf_t x, ulong v)
void arf_init_set_si(arf_t x, slong v)

Initialise `x` and set it to `v`.
:::

::: function
void arf_set_mpn(arf_t y, nn_srcptr x, slong xn, int sgnbit)

Sets `y` to the integer stored in the `xn` limbs at `x`, negated if
`sgnbit` is 1. Requires $\mathit{xn} > 0$ and `x[xn - 1]` nonzero.
:::

::: function
void arf_set_si_2exp_si(arf_t x, slong man, slong exp)
void arf_set_ui_2exp_si(arf_t x, ulong man, slong exp)
void arf_set_fmpz_2exp(arf_t x, const fmpz_t man, const fmpz_t exp)

Set `x` to $\mathit{man} \cdot 2^{\mathit{exp}}$ exactly.
:::

::: function
void arf_swap(arf_t y, arf_t x)

Swaps `x` and `y` without copying.
:::

::: function
int arf_set_round(arf_t y, const arf_t x, slong prec, arf_rnd_t rnd)
int arf_set_round_si(arf_t x, slong v, slong prec, arf_rnd_t rnd)
int arf_set_round_ui(arf_t x, ulong v, slong prec, arf_rnd_t rnd)
int arf_set_round_mpz(arf_t y, const mpz_t x, slong prec, arf_rnd_t rnd)
int arf_set_round_fmpz(arf_t y, const fmpz_t x, slong prec, arf_rnd_t rnd)

Set the output to the input rounded to `prec` bits in the direction
`rnd`, and return the inexact flag.
:::

::: function
int arf_set_round_fmpz_2exp(arf_t y, const fmpz_t x, const fmpz_t exp, slong prec, arf_rnd_t rnd)

Sets `y` to $x \cdot 2^{\mathit{exp}}$ rounded to `prec` bits in the
direction `rnd`, and returns the inexact flag.
:::

::: function
int arf_set_fmpq(arf_t y, const fmpq_t x, slong prec, arf_rnd_t rnd)

Sets `y` to the rational number `x` rounded to `prec` bits in the
direction `rnd`, and returns the inexact flag. This is
`arf_fmpz_div_fmpz` applied to the numerator and denominator.
:::

::: function
void arf_neg(arf_t y, const arf_t x)
void arf_abs(arf_t y, const arf_t x)

Set `y` to $-x$, respectively $|x|$. Exact.
:::

::: function
int arf_neg_round(arf_t y, const arf_t x, slong prec, arf_rnd_t rnd)

Sets `y` to $-x$ rounded to `prec` bits in the direction `rnd`, and
returns the inexact flag.
:::

::: function
void arf_mul_2exp_si(arf_t y, const arf_t x, slong e)
void arf_mul_2exp_fmpz(arf_t y, const arf_t x, const fmpz_t e)

Set `y` to $x \cdot 2^e$. Exact. This is the inverse of `arf_frexp`.
:::

::: function
void arf_min(arf_t z, const arf_t a, const arf_t b)
void arf_max(arf_t z, const arf_t a, const arf_t b)

Set `z` to the smaller, respectively larger, of `a` and `b`.
:::

## Conversion {#conversion}

::: function
void arf_get_fmpz_2exp(fmpz_t man, fmpz_t exp, const arf_t x)

Sets `man` and `exp` to the integers with $x = \mathit{man} \cdot
2^{\mathit{exp}}$ and `man` odd, for finite nonzero `x`. For $x = 0$ both
are set to 0. The result is undefined for infinite or NaN `x`.
:::

::: function
void arf_frexp(arf_t man, fmpz_t exp, const arf_t x)

Writes $x = m \cdot 2^e$ with $1/2 \le |m| < 1$, setting `man` to $m$ and
`exp` to $e$. A special value is copied to `man`, and `exp` is set to 0.
:::

::: function
double arf_get_d(const arf_t x, arf_rnd_t rnd)

Returns `x` rounded to a `double` in the direction `rnd`, with correct
rounding also on overflow (to the largest finite `double` or to
$\pm\infty$, depending on the direction) and on underflow (to zero or
the smallest subnormal).
:::

::: function
int arf_get_mpfr(mpfr_t x, const arf_t y, mpfr_rnd_t rnd)

Sets `x` to `y` rounded to the precision of `x` with the MPFR rounding
mode `rnd`, and returns the MPFR ternary value ($-1$, $0$ or $1$). An
exponent outside the current MPFR exponent range gives an infinity or a
zero and sets the MPFR overflow or underflow flag.
:::

::: function
int arf_get_fmpz(fmpz_t z, const arf_t x, arf_rnd_t rnd)

Sets `z` to `x` rounded to an integer in the direction `rnd` (ties to
even for `ARF_RND_NEAR`), and returns 1 if `x` is not an integer, 0
otherwise.

::: details
An infinite or NaN `x` raises an exception. The integer is computed
exactly, so a huge exponent can use all available memory or swap space
before the allocation fails; check the size of `x` (for example with
`arf_cmpabs_2exp_si`) before converting values of unknown magnitude.
:::
:::

::: function
slong arf_get_si(const arf_t x, arf_rnd_t rnd)

Returns `x` rounded to an integer in the direction `rnd` (ties to even
for `ARF_RND_NEAR`). Raises an exception if `x` is infinite, NaN, or if
the rounded value does not fit in an `slong`.
:::

::: function
int arf_get_fmpz_fixed_fmpz(fmpz_t y, const arf_t x, const fmpz_t e)
int arf_get_fmpz_fixed_si(fmpz_t y, const arf_t x, slong e)

Set `y` to the integer $\operatorname{trunc}(x / 2^e)$, so that $y \cdot
2^e$ is `x` truncated towards zero to a multiple of $2^e$, and return 1 if
truncation occurred. The size caveats of `arf_get_fmpz` apply.
:::

::: function
void arf_floor(arf_t z, const arf_t x)
void arf_ceil(arf_t z, const arf_t x)
void arf_nint(arf_t z, const arf_t x)

Set `z` to $\lfloor x \rfloor$, $\lceil x \rceil$ and the nearest integer
to `x` with ties to even. The result is exact and never needs more bits
than `x`; round it with `arf_set_round` if needed.
:::

::: function
void arf_get_fmpq(fmpq_t y, const arf_t x)

Sets `y` to the exact rational value of `x`. Raises an exception for
infinite or NaN `x`; the size caveats of `arf_get_fmpz` apply.
:::

## Comparison and bounds {#comparison}

::: function
int arf_equal(const arf_t x, const arf_t y)
int arf_equal_si(const arf_t x, slong y)
int arf_equal_ui(const arf_t x, ulong y)
int arf_equal_d(const arf_t x, double y)

Return whether `x` and `y` are the same number. Unlike IEEE 754
comparison, NaN is equal to itself. For a `double` `y`, both zeros
compare equal to 0 and every NaN to NaN.
:::

::: function
int arf_cmp(const arf_t x, const arf_t y)
int arf_cmp_si(const arf_t x, slong y)
int arf_cmp_ui(const arf_t x, ulong y)
int arf_cmp_d(const arf_t x, double y)

Return a negative number, zero or a positive number according as $x <
y$, $x = y$ or $x > y$. If either argument is NaN, the result is 0.
:::

::: function
int arf_cmpabs(const arf_t x, const arf_t y)
int arf_cmpabs_ui(const arf_t x, ulong y)
int arf_cmpabs_d(const arf_t x, double y)
int arf_cmpabs_mag(const arf_t x, const mag_t y)
int arf_mag_cmpabs(const mag_t x, const arf_t y)

Compare $|x|$ with $|y|$, with the same return convention as `arf_cmp`.
:::

::: function
int arf_cmp_2exp_si(const arf_t x, slong e)
int arf_cmpabs_2exp_si(const arf_t x, slong e)

Compare $x$, respectively $|x|$, with $2^e$.
:::

::: function
int arf_sgn(const arf_t x)

Returns $-1$, $0$ or $1$ according to the sign of `x`, and 0 for NaN.
:::

::: function
slong arf_bits(const arf_t x)

Returns the number of bits of the odd mantissa of `x`, that is, the
smallest precision at which `x` is exactly representable, and 0 for a
special value.
:::

::: function
void arf_bot(fmpz_t e, const arf_t x)

Sets `e` to the exponent of the lowest nonzero bit of `x`, that is, to
the exponent of `arf_get_fmpz_2exp`, so that $x$ is an odd multiple of
$2^e$. Sets `e` to 0 for a special value.
:::

::: function
int arf_is_int(const arf_t x)
int arf_is_int_2exp_si(const arf_t x, slong e)

Return whether `x` is an integer, respectively an integer multiple of
$2^e$.
:::

::: function
void arf_abs_bound_lt_2exp_fmpz(fmpz_t b, const arf_t x)
void arf_abs_bound_le_2exp_fmpz(fmpz_t b, const arf_t x)

Set `b` to the smallest integer with $|x| < 2^b$, respectively $|x| \le
2^b$, for finite nonzero `x`. The result is undefined for special values.
:::

::: function
slong arf_abs_bound_lt_2exp_si(const arf_t x)

Returns the smallest integer $b$ with $|x| < 2^b$, clamped to
$[-\mathrm{ARF\_PREC\_EXACT}, \mathrm{ARF\_PREC\_EXACT}]$. Returns
$-\mathrm{ARF\_PREC\_EXACT}$ for zero and $\mathrm{ARF\_PREC\_EXACT}$ for
infinities and NaN.
:::

## Conversion to and from mag {#mag}

::: function
void arf_get_mag(mag_t y, const arf_t x)
void arf_get_mag_lower(mag_t y, const arf_t x)

Set `y` to an upper, respectively lower, bound for $|x|$. An infinite or
NaN `x` gives $+\infty$ for the upper bound.
:::

::: function
void arf_set_mag(arf_t y, const mag_t x)

Sets `y` to the value of `x`. Exact.
:::

::: function
void mag_init_set_arf(mag_t y, const arf_t x)

Initialises `y` and sets it to an upper bound for $|x|$.
:::

::: function
void mag_fast_init_set_arf(mag_t y, const arf_t x)

Initialises `y` to an upper bound for $|x|$, for finite `x` whose exponent
is small enough for the `mag_fast_` functions. Nothing is checked.
:::

::: function
void arf_mag_set_ulp(mag_t z, const arf_t y, slong prec)

Sets `z` to $2^{e - \mathit{prec}}$, where $e$ is the exponent of `y`
($2^{e-1} \le |y| < 2^e$): one unit in the last place of a `prec`-bit
number of the size of `y`. Raises an exception for a special `y`.
:::

::: function
void arf_mag_add_ulp(mag_t z, const mag_t x, const arf_t y, slong prec)
void arf_mag_fast_add_ulp(mag_t z, const mag_t x, const arf_t y, slong prec)

Set `z` to an upper bound for $x + 2^{e - \mathit{prec}}$, with $e$ the
exponent of `y`. This adds the rounding error of a `prec`-bit operation
with result `y` to the radius `x`. The fast version requires small
exponents; the other raises an exception for a special `y`.
:::

## Shallow copies {#shallow}

::: function
void arf_init_set_shallow(arf_t z, const arf_t x)
void arf_init_set_mag_shallow(arf_t y, const mag_t x)
void arf_init_neg_shallow(arf_t z, const arf_t x)
void arf_init_neg_mag_shallow(arf_t z, const mag_t x)

Initialise the output as a shallow copy of `x` or of $-x$, sharing its
limbs, without allocating memory.

::: details
The copy may only be used as a constant input, must not be cleared, and
is invalid once `x` is modified or cleared. The `mag_t` versions make it
possible to call `arf` functions on radii, as in
`arf_cmpabs_mag`.
:::
:::

## Arithmetic {#arithmetic}

::: function
int arf_add(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd)
int arf_add_si(arf_ptr z, arf_srcptr x, slong y, slong prec, arf_rnd_t rnd)
int arf_add_ui(arf_ptr z, arf_srcptr x, ulong y, slong prec, arf_rnd_t rnd)
int arf_add_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, slong prec, arf_rnd_t rnd)

Set `z` to $x + y$, correctly rounded.

::: details
The cost depends on the precision and not on the distance between the
exponents: when `y` is far below the last bit of `x` at precision `prec`,
only its sign and position matter for the rounding.
:::
:::

::: function
int arf_add_fmpz_2exp(arf_ptr z, arf_srcptr x, const fmpz_t y, const fmpz_t exp, slong prec, arf_rnd_t rnd)

Sets `z` to $x + y \cdot 2^{\mathit{exp}}$, correctly rounded.
:::

::: function
int arf_sub(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd)
int arf_sub_si(arf_ptr z, arf_srcptr x, slong y, slong prec, arf_rnd_t rnd)
int arf_sub_ui(arf_ptr z, arf_srcptr x, ulong y, slong prec, arf_rnd_t rnd)
int arf_sub_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, slong prec, arf_rnd_t rnd)

Set `z` to $x - y$, correctly rounded.
:::

::: function
int arf_mul(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd)
int arf_mul_ui(arf_ptr z, arf_srcptr x, ulong y, slong prec, arf_rnd_t rnd)
int arf_mul_si(arf_ptr z, arf_srcptr x, slong y, slong prec, arf_rnd_t rnd)
int arf_mul_mpz(arf_ptr z, arf_srcptr x, const mpz_t y, slong prec, arf_rnd_t rnd)
int arf_mul_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, slong prec, arf_rnd_t rnd)

Set `z` to $xy$, correctly rounded.

::: details
`arf_mul` is a macro that calls `arf_mul_rnd_down` for `ARF_RND_DOWN`
and `arf_mul_rnd_any` otherwise. Rounding down is the mode used by
ball arithmetic, and `arf_mul_rnd_down` has dedicated code for one- and
two-limb operands (`arf_mul_rnd_any` rounds those in registers too).
Longer products are formed in full and rounded when the full product
is cheap or not much longer than the output (up to 16 limbs in total,
and the shapes of `_arf_mul_prefer_full`); otherwise only their top
limbs are computed, with a guard limb (the high products of
[mpn_extras](mpn_extras#short), `flint_mpn_mulhigh_normalised` and the
balanced kernels), and the result is rounded from the guard limb when
it certifies the rounding (`_arf_mpn_round_guarded`), the full product
deciding otherwise (probability about $2^{-55}$ for random operands,
and always for a product that is exact at `prec` bits). MPFR is no
longer used.
:::
:::

::: function
int arf_neg_mul(arf_t z, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)

Sets `z` to $-xy$, correctly rounded.
:::

::: function
int arf_addmul(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd)
int arf_addmul_ui(arf_ptr z, arf_srcptr x, ulong y, slong prec, arf_rnd_t rnd)
int arf_addmul_si(arf_ptr z, arf_srcptr x, slong y, slong prec, arf_rnd_t rnd)
int arf_addmul_mpz(arf_ptr z, arf_srcptr x, const mpz_t y, slong prec, arf_rnd_t rnd)
int arf_addmul_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, slong prec, arf_rnd_t rnd)
int arf_submul(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd)
int arf_submul_ui(arf_ptr z, arf_srcptr x, ulong y, slong prec, arf_rnd_t rnd)
int arf_submul_si(arf_ptr z, arf_srcptr x, slong y, slong prec, arf_rnd_t rnd)
int arf_submul_mpz(arf_ptr z, arf_srcptr x, const mpz_t y, slong prec, arf_rnd_t rnd)
int arf_submul_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, slong prec, arf_rnd_t rnd)

Set `z` to $z + xy$, respectively $z - xy$, with a single rounding (fused
multiply-add).
:::

::: function
int arf_fma(arf_ptr res, arf_srcptr x, arf_srcptr y, arf_srcptr z, slong prec, arf_rnd_t rnd)

Sets `res` to $xy + z$ with a single rounding.

::: details
`arf_addmul`, `arf_submul` and `arf_fma` share one core: the exact
product and one exact addition, rounded once. From 24 limbs of
precision and of each factor, the product is truncated instead
(`flint_mpn_mulhigh`: its top limbs above the rounding point with two
guard limbs, a lower bound of the product within a known number of
units of the lowest limb), the exact sum of `z` and that bound is
formed, and it is rounded when its bits between the rounding point and
the error band certify the rounding (`_arf_rounding_certain`); the
exact product decides otherwise.
:::
:::

::: function
int arf_sosq(arf_t z, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)

Sets `z` to $x^2 + y^2$ with a single rounding.
:::

::: function
int arf_div(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd)
int arf_div_ui(arf_ptr z, arf_srcptr x, ulong y, slong prec, arf_rnd_t rnd)
int arf_ui_div(arf_ptr z, ulong x, arf_srcptr y, slong prec, arf_rnd_t rnd)
int arf_div_si(arf_ptr z, arf_srcptr x, slong y, slong prec, arf_rnd_t rnd)
int arf_si_div(arf_ptr z, slong x, arf_srcptr y, slong prec, arf_rnd_t rnd)
int arf_div_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, slong prec, arf_rnd_t rnd)
int arf_fmpz_div(arf_ptr z, const fmpz_t x, arf_srcptr y, slong prec, arf_rnd_t rnd)
int arf_fmpz_div_fmpz(arf_ptr z, const fmpz_t x, const fmpz_t y, slong prec, arf_rnd_t rnd)

Set `z` to $x / y$, correctly rounded for the strict modes. Division by
zero gives NaN. The relaxed modes `ARF_RND_FAST` and `ARF_RND_ACCURATE`
are accepted.

::: details
The quotient kernel `_arf_mpn_div` (see [](#mpn)): exact register
divisions at up to two limbs and for one-limb divisors (such as the
integers of `arf_div_ui`), the floor with the remainder deciding the
rounding; otherwise `flint_mpn_divhigh` (see [](mpn_extras#fphigh)):
the rounding is certified by the guard limb unless it lies within the
error bound of a rounding boundary (probability about $2^{-55}$, but
always for exact quotients), in which case the exact truncated quotient
and remainder (`flint_mpn_divhigh_exact`) decide the rounding and the
inexact flag. Division by a power of two is a rounding of the
numerator.
:::
:::

::: function
int arf_sqrt(arf_ptr z, arf_srcptr x, slong prec, arf_rnd_t rnd)
int arf_sqrt_ui(arf_t z, ulong x, slong prec, arf_rnd_t rnd)
int arf_sqrt_fmpz(arf_t z, const fmpz_t x, slong prec, arf_rnd_t rnd)

Set `z` to $\sqrt{x}$, correctly rounded for the strict modes. A negative
`x` gives NaN. The relaxed modes are accepted.

::: details
As `arf_div`, with `_arf_mpn_sqrt`: the exact truncated root of the
mantissa padded to $2 \, \mathit{rn}$ limbs (`flint_mpn_sqrthigh_exact`,
with the half bit for rounding to nearest at a limb boundary) up to the
Newton range, which is cheaper than the guard limb; the one-limb root
in registers (`_flint_mpn_isqrt_2_rem`), here inlined for one-limb
operands at up to one limb of precision. A power of two with an odd
exponent is exact.
:::
:::

::: function
int arf_rsqrt(arf_ptr z, arf_srcptr x, slong prec, arf_rnd_t rnd)

Sets `z` to $1/\sqrt{x}$, correctly rounded for the strict modes. A
negative `x` gives NaN and $x = 0$ gives $+\infty$. The relaxed modes
are accepted.

::: details
As `arf_div`, with `_arf_mpn_rsqrt` (`flint_mpn_rsqrthigh` and
`flint_mpn_rsqrthigh_exact`). The result is dyadic only for a power of
four, handled exactly; all other results are inexact.
:::
:::

::: function
int arf_root(arf_ptr z, arf_srcptr x, ulong k, slong prec, arf_rnd_t rnd)

Sets `z` to $x^{1/k}$, correctly rounded. A negative `x` gives NaN.

::: details
This is a wrapper of the MPFR root function, which becomes slow and uses
much memory for large $k$; `arb_root_ui` is usually a better choice.
:::
:::

## Sums and dot products {#sums}

::: function
int arf_sum(arf_t s, arf_srcptr terms, slong len, slong prec, arf_rnd_t rnd)

Sets `s` to the sum of the `len` entries of `terms` with a single
rounding, and returns the inexact flag.

::: details
Unlike repeated calls to `arf_add` with `ARF_PREC_EXACT`, this does not
build exact intermediate sums whose size depends on the spread of the
exponents: terms that are far below the current precision are replaced
by tiny perturbations that only affect the rounding. The implementation
merges terms pairwise and takes quadratic time in `len` in the worst case.
:::
:::

::: function
int arf_dot(arf_t res, const arf_t initial, int subtract, arf_srcptr x, slong xstep, arf_srcptr y, slong ystep, slong len, slong prec, arf_rnd_t rnd)

Sets `res` to $s + (-1)^{\mathit{subtract}} \sum_{i=0}^{\mathit{len}-1}
x_i y_i$ with a single rounding, where $s$ is `initial` (0 if `initial` is
`NULL`) and the entries are read with strides `xstep` and `ystep`. Returns
the inexact flag.

::: details
The products are formed exactly and summed with `arf_sum`. This is a
reference implementation, not optimised for speed; `arf_approx_dot` and
`arb_dot` are much faster.
:::
:::

::: function
void arf_approx_dot(arf_t res, const arf_t initial, int subtract, arf_srcptr x, slong xstep, arf_srcptr y, slong ystep, slong len, slong prec, arf_rnd_t rnd)

Computes the same dot product as `arf_dot` approximately: the final
rounding is in the direction `rnd`, but the intermediate roundings are
not specified, so neither correct rounding nor an error bound is
guaranteed. It uses the fixed-point accumulation of `arb_dot` without the
error bound, and is the midpoint computation of
`arb_approx_dot`.
:::

## Complex arithmetic {#complex}

::: function
int arf_complex_mul(arf_t e, arf_t f, const arf_t a, const arf_t b, const arf_t c, const arf_t d, slong prec, arf_rnd_t rnd)
int arf_complex_mul_fallback(arf_t e, arf_t f, const arf_t a, const arf_t b, const arf_t c, const arf_t d, slong prec, arf_rnd_t rnd)

Set $e + fi = (a + bi)(c + di)$ with both parts correctly rounded. Bit 0
of the return value is the inexact flag of $e$ and bit 1 that of $f$.

::: details
If one of $a, b, c, d$ is zero, the product is computed with two real
multiplications and no additions, even if another component is infinite
or NaN; special values follow from this rule. From 13 limbs of
precision (up to 2000 limbs), for balanced parts of comparable
magnitude, the complex high product (`flint_mpn_mulhigh_n_complex`)
gives both parts within a few units of a guard limb, and the rounding
is certified as in `arf_mul`. Otherwise `arf_complex_mul` forms the
exact products and rounds $ac - bd$ and $ad + bc$ once each: when all
four mantissas have at least 20 limbs and the exponents of the real
and imaginary parts are within 64 bits of each other, with the exact
complex product `flint_mpn_mul_complex` (Karatsuba's three products,
or shared `fft_small` transforms), else with four real products.
`arf_complex_mul_fallback` is the plain version with four products,
kept for testing.
:::
:::

::: function
int arf_complex_sqr(arf_t e, arf_t f, const arf_t a, const arf_t b, slong prec, arf_rnd_t rnd)

Sets $e + fi = (a + bi)^2$, with the same conventions as
`arf_complex_mul`, but faster.

::: details
As `arf_complex_mul` with the square kernels: the complex high square
(`flint_mpn_sqrhigh_n_complex`), and for the exact square
`flint_mpn_sqr_complex` from 24 limbs, below that two real squares and
one product.
:::
:::

## Input and output {#io}

::: function
void arf_print(const arf_t x)
void arf_fprint(FILE * file, const arf_t x)

Print the exact value of `x` as `(m * 2^e)` with an odd integer mantissa,
or as `(0)`, `(+inf)`, `(-inf)` or `(nan)`, to standard output or to
`file`.
:::

::: function
void arf_printd(const arf_t y, slong d)
void arf_fprintd(FILE * file, const arf_t y, slong d)
char * arf_get_str(const arf_t x, slong d)

Print `x` in decimal with `d` significant digits, or return such a
string (to be freed with `flint_free`). The decimal conversion is
faithful: the error is at most one unit in the last digit.
:::

::: function
char * arf_get_str_rnd(const arf_t x, slong d, arf_rnd_t rnd)

Return `x` in decimal with `d` significant digits, correctly rounded in
the direction `rnd`, as a string to be freed with `flint_free`.
:::

::: function
void arf_debug(const arf_t x)

Prints the internal fields of `x` (exponent, size, sign and limbs).
:::

::: function
char * arf_dump_str(const arf_t x)
int arf_load_str(arf_t res, const char * data)

`arf_dump_str` returns a string, to be freed with `flint_free`, holding
the exact value of `x` as two hexadecimal integers (mantissa and
exponent), which `arf_load_str` reads back exactly. `arf_load_str`
returns nonzero if the string is not in this format.
:::

::: function
int arf_dump_file(FILE* stream, const arf_t x)
int arf_load_file(arf_t res, FILE *stream)

Write `x` to `stream` in the format of `arf_dump_str`, and read it back.
Both return nonzero on failure. Several values written in sequence must
be separated by whitespace.
:::

## Random generation {#random}

::: function
void arf_randtest(arf_t x, flint_rand_t state, slong bits, slong mag_bits)
void arf_randtest_not_zero(arf_t x, flint_rand_t state, slong bits, slong mag_bits)
void arf_randtest_special(arf_t x, flint_rand_t state, slong bits, slong mag_bits)

Set `x` to a random finite number with a mantissa of at most `bits` bits
and an exponent of at most `mag_bits` bits, with a non-uniform
distribution that favours corner cases. The second version never returns
zero; the third sometimes returns $\pm\infty$ or NaN.
:::

::: function
void arf_urandom(arf_t x, flint_rand_t state, slong bits, arf_rnd_t rnd)

Sets `x` to a uniformly distributed random number in $[0, 1]$ with `bits`
bits, obtained by rounding a random integer in the direction `rnd`.
:::

## Rounding helpers {#rounding-helpers}

::: function
int arf_rounds_down(arf_rnd_t rnd, int sgnbit)
int arf_rounds_up(arf_rnd_t rnd, int sgnbit)

Return whether the strict mode `rnd` rounds a number with sign bit
`sgnbit` (1 for negative) towards zero, respectively away from zero.
Both return 0 for `ARF_RND_NEAR`.
:::

::: function
mpfr_rnd_t arf_rnd_to_mpfr(arf_rnd_t rnd)

Returns the MPFR rounding mode corresponding to `rnd` (`MPFR_RNDZ`,
`MPFR_RNDA`, `MPFR_RNDD`, `MPFR_RNDU`, and `MPFR_RNDN` for any other
mode). Available only if `mpfr.h` is included before `arf.h`.
:::

::: function
arf_rnd_t arf_rnd_relaxed_to_strict(arf_rnd_t rnd)

Returns the strict mode by which the standard algorithms implement a
relaxed mode: `ARF_RND_DOWN` for `ARF_RND_FAST` and `ARF_RND_NEAR` for
`ARF_RND_ACCURATE`. Other modes are returned unchanged.
:::

## Internal functions {#internal}

These functions are used to implement the module and ball arithmetic.
Their interfaces may change.

::: function
void arf_mul_special(arf_t z, const arf_t x, const arf_t y)

Sets `z` to $xy$ when at least one of `x` and `y` is a special value.
:::

::: function
int arf_mul_rnd_any(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd)
int arf_mul_rnd_down(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec)

The implementations of `arf_mul`: for any rounding mode, and for
`ARF_RND_DOWN`.
:::

::: function
int _arf_mul_prefer_full(slong xn, slong yn, slong n)

Whether the product of mantissas of `xn` and `yn` limbs of which the
top `n` limbs are wanted is computed in full rather than as a high
product: small products (up to `ARF_MUL_FULL_SMALL` = 64 limb
products), and the shapes for which the general high product would
form the full product anyway (`_flint_mpn_mulhigh_prefer_full`).
:::

::: function
int _arf_set_round_ui(arf_t x, ulong v, int sgnbit, slong prec, arf_rnd_t rnd)

Sets `x` to $\pm v$ (negative if `sgnbit` is 1) rounded to `prec` bits.
:::

::: function
int _arf_set_round_uiui(arf_t z, slong * fix, ulong hi, ulong lo, int sgnbit, slong prec, arf_rnd_t rnd)

Sets the mantissa of `z` to the two-limb value $(\mathit{hi},
\mathit{lo})$, negated if `sgnbit` is 1, rounded to `prec` bits. Requires
a nonzero value. Writes the exponent shift to `fix` without setting the
exponent of `z`.
:::

::: function
int _arf_set_round_mpn(arf_t y, slong * exp_shift, nn_srcptr x, slong xn, int sgnbit, slong prec, arf_rnd_t rnd)

Sets the mantissa of `y` to the `xn` limbs at `x` (top limb nonzero),
negated if `sgnbit` is 1, rounded to `prec` bits, and returns the inexact
flag. Writes the exponent shift caused by leading zero bits and by
rounding up to a power of two to `exp_shift` without setting the exponent
of `y`. `x` may point into the limbs of `y`.
:::

::: function
int _arf_set_mpn_fixed(arf_t z, nn_srcptr xp, slong xn, slong fixn, int negative, slong prec, arf_rnd_t rnd)

Sets `z` to the fixed-point number with `xn` limbs of which `fixn` are
fractional, negated if `negative` is set, rounded to `prec` bits.
:::

::: function
int _arf_get_integer_mpn(nn_ptr y, nn_srcptr x, slong xn, slong exp)

Writes the integer part of the number with mantissa $(x, \mathit{xn})$
(top bit set, `x[0]` nonzero) and exponent `exp` to `y`, using exactly as
many limbs as needed, and returns whether a fractional part was
discarded.
:::

::: function
int _arf_add_mpn(arf_t z, nn_srcptr xp, slong xn, int xsgnbit, const fmpz_t xexp, nn_srcptr yp, slong yn, int ysgnbit, flint_bitcnt_t shift, slong prec, arf_rnd_t rnd)

Adds two mantissas, the second shifted right by `shift` bits relative to
the first, which has exponent `xexp`, and rounds the result into `z`.
This is the kernel of `arf_add`.
:::

::: function
int _arf_call_mpfr_func(arf_ptr r1, arf_ptr r2, int (*func)(void), arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd)

Evaluates the MPFR function `func` (cast to a generic function pointer)
on `x` (and `y`, if not `NULL`) at precision `prec`, writing one result
or, if `r2` is not `NULL`, two results. Returns the MPFR ternary value.
:::

::: function
void _arf_promote(arf_t x, slong n)
void _arf_demote(arf_t x)

Switch the mantissa storage of `x` to a heap array of `n` limbs, and back
to inline storage.
:::

::: function
void _arf_set_mpn_canonical(arf_ptr z, slong n, int sgnbit)
void _arf_set_mpn_canonical_copy(arf_ptr z, nn_srcptr d, slong n, int sgnbit)

Set the mantissa of `z` from the `n` normalised limbs a kernel wrote at
the start of its buffer (allocated with `ARF_GET_MPN_WRITE` for at
least `n` limbs), respectively from the limbs `d` in scratch: zero low
limbs are stripped and a result of two or fewer limbs moves into the
inline limbs, so that `z` is in canonical form. The exponent is left
for the caller. `ARF_KERNEL_DECL`, `ARF_KERNEL_ALLOC(zd, alloc, z, x, y)`
and `ARF_KERNEL_FINISH(z, zd, n, sgnbit)` wrap the buffer management
for a kernel writing `alloc` limbs into `z`: the buffer of `z` itself,
unless `z` is one of the operands `x`, `y` or the result may be short
enough for the inline limbs, in which case the kernel works in scratch
and the result is copied.
:::

## Limb-array kernels {#mpn}

The arithmetic of long mantissas is done by kernels on limb arrays,
shared with [nfloat](nfloat) and [acf](acf). A mantissa is a
normalised limb array $\{x, \mathit{xn}\}$ (top bit set), standing for
the value
$$x \cdot 2^{e - \mathit{xn} \cdot \bits}$$
for an exponent $e$ which stays with the caller, as does the sign. A
kernel writes its result into $\{z, \mathit{rn}\}$ with $\mathit{rn}
= \lceil \mathit{prec} / \bits \rceil$ limbs, normalised with the bits
below `prec` zero, and sets `*zexp` to the *correction* of the natural
exponent of the operation ($e_x - e_y$ for a quotient, and so on); the
caller adds it to that exponent. The rounded kernels take a sign-free
mode (`ARF_RND_DOWN`, `ARF_RND_UP` or `ARF_RND_NEAR`; `_arf_rnd_resolve`
gives it) and return the inexact flag: inside, they form an
approximation with a guard limb, round it when the guard limb certifies
that the exact value rounds the same way, and compute the exact value
otherwise. The `_approx` kernels write $\mathit{rn}$ limbs within an
ulp of the magnitude of the result and nothing else: the form for
`nfloat`, which rounds nothing. No kernel allocates memory visibly:
scratch beyond a small stack buffer (`ARF_MPN_STACK_LIMBS` limbs,
`ARF_MPN_TMP_DECL`, `ARF_MPN_TMP_ALLOC`) comes from a thread-local
buffer (`_arf_mpn_tmp`, one buffer, not for nested use within the
kernels).

::: function
arf_rnd_t _arf_rnd_resolve(arf_rnd_t rnd, int sgnbit)

The sign-free strict mode rounding a value of sign bit `sgnbit` in the
mode `rnd`: `ARF_RND_FLOOR` and `ARF_RND_CEIL` become `ARF_RND_DOWN` or
`ARF_RND_UP`, the relaxed modes become the strict modes of
`arf_rnd_relaxed_to_strict`, and the others are returned as they are.
:::

::: function
int _arf_mpn_set_round(nn_ptr res, slong * exp, nn_srcptr x, slong xn, slong prec, arf_rnd_t rnd)

Rounds the exact value $(x / B^{\mathit{xn}}) 2^E$, where `x` has a
nonzero top limb (leading zero bits allowed), to `prec` bits in the
sign-free mode `rnd`, writing the $\mathit{rn}$ normalised limbs `res`
with the bits below `prec` zero: the result is $(\mathit{res} /
B^{\mathit{rn}}) 2^{E + \mathit{exp}}$. Returns the inexact flag. `res`
may be `x` or lie below it in the same buffer (in-place rounding; then
$\mathit{xn} \ge \mathit{rn}$).
:::

::: function
int _arf_mpn_round_guarded(nn_ptr res, slong * exp, nn_srcptr a, slong n, ulong g, int e, slong prec, arf_rnd_t rnd, int side)

Given the $n$ normalised limbs `a` and the limb `g` below them, which
together approximate a value to within $2^e$ ulps of `g` (in either
direction for `side` = 0, strictly above the approximation for `side` =
1, strictly below it for `side` = $-1$), rounds the value to
$\mathit{prec} \le n \cdot \bits$ bits when the guard limb certifies
the rounding: writes the $\lceil \mathit{prec} / \bits \rceil$ limbs to
`res` (which may be `a`, or lie below it in the same buffer), sets `exp`
to the change of the exponent ($-1$, 0 or 1) and returns 1, the value
being inexact. Returns $-1$ when the window does not certify the
rounding; the exact value must then be computed.

::: details
With the approximation $A$ within $2^e$ ulps of the exact value $P$ and
$k \ge \bits$ discarded bits, the bits $e, \ldots, \bits - 2$ of `g`
being neither all zero nor all ones gives $2^e \le A \bmod 2^{k-1} <
2^{k-1} - 2^e$, so that $P$ lies strictly between the same rounding
boundaries as $A$ (the truncation point and the half point) and is not
exact; the rounding of $A$ is then the rounding of $P$ in every mode
[@AJ2025]. The check fails with probability about $2^{e + 1 - \bits}$.
With the direction of the error known, only the boundary in that
direction matters: all zeros will do for a value above the
approximation (which then rounds like the approximation with sticky
bits), all ones for one below; and a value just below an approximation
that is exact at `prec` bits rounds to it, or to the value below it when
rounding toward zero.
:::
:::

::: function
int _arf_rounding_certain(nn_srcptr a, slong an, slong e, slong prec)

Whether rounding $\{a, \mathit{an}\}$ (an approximation of a value $P$
within $2^e$ units of $a_0$, on either side) to `prec` bits gives the
same result as rounding $P$, in every strict rounding mode, with $P$
inexact: the bits of `a` from bit $e$ of its lowest limb up to the bit
below the rounding point (at most one limb of them) are neither all
zero nor all ones. The same argument as above, for an approximation
whose error band is given by a bit position rather than a guard limb.
:::

::: function
int _arf_mpn_increment(nn_ptr res, slong rn, slong sb)
int _arf_mpn_decrement(nn_ptr res, slong rn, slong sb)

The $\mathit{rn}$ normalised limbs `res`, whose bits below `sb` are
zero, moved up or down by one unit of $2^{\mathit{sb}}$. The increment
returns the carry (1 when the result is a power of two, with the
exponent increased); the decrement returns $-1$ when `res` was a power
of two, the result being the $\mathit{rn} \cdot \bits - \mathit{sb}$
ones just below it, normalised (the exponent decreases), else 0.
:::

::: function
void _arf_mpn_mul_approx(nn_ptr z, slong * zexp, nn_srcptr x, slong xn, nn_srcptr y, slong yn, slong prec)

The product $xy$ truncated to $\mathit{rn}$ limbs, within an ulp of
$\mathit{rn}$ limbs: the exact product when it fits or the shape prefers
it (`_arf_mul_prefer_full`), else the normalised high product with its
guard limb dropped. `z` must not alias an operand.
:::

::: function
int _arf_mpn_div(nn_ptr z, slong * zexp, nn_srcptr x, slong xn, nn_srcptr y, slong yn, slong prec, arf_rnd_t rnd)
int _arf_mpn_inv(nn_ptr z, slong * zexp, nn_srcptr y, slong yn, slong prec, arf_rnd_t rnd)
int _arf_mpn_sqrt(nn_ptr z, slong * zexp, nn_srcptr x, slong xn, int t, slong prec, arf_rnd_t rnd)
int _arf_mpn_rsqrt(nn_ptr z, slong * zexp, nn_srcptr x, slong xn, int t, slong prec, arf_rnd_t rnd)

The quotient $x / y$, the reciprocal $1 / y$, the square root of $x$
and its reciprocal, correctly rounded to `prec` bits in the sign-free
mode `rnd` into $\{z, \mathit{rn}\}$, returning the inexact flag, with
`*zexp` the correction to the natural exponent: $e_x - e_y$ for the
quotient (the correction 0 or 1), $-e_y$ for the reciprocal (1 or 2),
and for the roots of $x = A \cdot 2^e$ with $t = e \bmod 2$ given,
$(e + t) / 2$ for $\sqrt{x}$ and $(1 - t) - (e - t) / 2$ for
$1 / \sqrt{x}$ (the corrections 0 or 1).

The output may alias an operand. Zero low limbs of the mantissas
(nfloat pads its mantissas to a fixed length) are stripped first, so
that integers and other short operands take the fast paths.

::: details
Small sizes are exact and in registers (`mpn_extras/fp_small.h`, with
`arf/mpn_small.h` for callers with fixed small sizes): a quotient of
one limb from one-limb operands by one $2/1$ division, a quotient of
up to two limbs from operands of up to two (divisor: three) limbs by
two or three $3/2$ divisions, each giving the floor with the remainder,
from which the rounding is decided (the remainder encoded as the limb
below the floor: its top bit the half bit, a low bit the sticky bit); a
divisor of a single limb at any size by a chain of $2/1$ divisions
(`flint_mpn_divrem_1_preinv`, GMP's `mpn_divrem_1` from eight limbs).
The one-limb square root is the exact integer root of the two-limb
radicand (`_flint_mpn_isqrt_2_rem`, the remainder deciding the half
bit); the two-limb roots and the one-limb reciprocal root come from the
three- resp. two-limb register kernels, certified when their guard limb
lies away from the rounding boundaries.

Longer quotients and reciprocal roots use the guard-limb kernels of
`mpn_extras` (`flint_mpn_divhigh`, `flint_mpn_rsqrthigh`, see
[](mpn_extras#fphigh)): $\mathit{rn} + 1$ limbs with the top bit set,
the lowest a guard limb within a small bound of the exact value in
either direction, certified by `_arf_mpn_round_guarded`, with the exact
floor of the same quantity (the `_exact` variants, with an inexact
flag) as the fallback. Longer square roots are exact integer roots of
the mantissa padded to $2 \, \mathit{rn}$ limbs
(`flint_mpn_sqrthigh_exact`; for rounding to nearest at a limb boundary
the variant with the half bit up to 24 limbs, the guard-limb root
beyond), which is cheaper than the guard limb below the Newton range.
A divisor that is a power of two shifts; a root of a power of two with
an odd exponent is a power of two.
:::
:::

::: function
void _arf_mpn_div_approx(nn_ptr z, slong * zexp, nn_srcptr x, slong xn, nn_srcptr y, slong yn, slong prec)
void _arf_mpn_inv_approx(nn_ptr z, slong * zexp, nn_srcptr y, slong yn, slong prec)
void _arf_mpn_sqrt_approx(nn_ptr z, slong * zexp, nn_srcptr x, slong xn, int t, slong prec)
void _arf_mpn_rsqrt_approx(nn_ptr z, slong * zexp, nn_srcptr x, slong xn, int t, slong prec)

As above, but the $\mathit{rn}$ limbs are written without rounding and
without a flag (the kernels of nfloat without directed rounding): the
quotients are the truncation or one ulp more (an error below an ulp,
up to the truncation of operands beyond $\mathit{rn} + 2$ limbs; the
truncation itself for a divisor with a single nonzero limb), the square
root is the truncation, and the reciprocal square root is within
$1 + 2^{-37}$ ulp.

::: details
The quotient uses the register and chain paths above without the
rounding step, and otherwise the schoolbook approximate division at
$\mathit{rn}$ limbs (`_flint_mpn_divapprox_basecase_preinv1`, the floor
or one more) up to 80 limbs, with the divisor truncated to
$\mathit{rn} + 2$ limbs and the dividend to $\mathit{rn} + \mathit{yn}$;
beyond that, and for two-limb divisors, the guard-limb quotient
truncated. The square root is the exact root up to the Newton range,
and the reciprocal square root the one- or two-limb register kernel
without certification, or the guard-limb root truncated.
:::
:::

::: function
void _arf_mpn_cinv_approx(nn_ptr re, slong * reexp, nn_ptr im, slong * imexp, nn_srcptr a, slong an, slong ea, nn_srcptr b, slong bn, slong eb, slong prec)
void _arf_mpn_cdiv_approx(nn_ptr re, slong * reexp, int * reneg, nn_ptr im, slong * imexp, int * imneg, nn_srcptr a, slong an, slong ea, int sa, nn_srcptr b, slong bn, slong eb, int sb, nn_srcptr c, slong cn, slong ec, int sc, nn_srcptr d, slong dn, slong ed, int sd, slong prec)
void _arf_mpn_csqrt_approx(nn_ptr re, slong * reexp, nn_ptr im, slong * imexp, nn_srcptr a, slong an, slong ea, int sa, nn_srcptr b, slong bn, slong eb, slong prec)
void _arf_mpn_crsqrt_approx(nn_ptr re, slong * reexp, nn_ptr im, slong * imexp, nn_srcptr a, slong an, slong ea, int sa, nn_srcptr b, slong bn, slong eb, slong prec)

The reciprocal, the quotient, the square root and the reciprocal
square root of complex numbers given as mantissa/exponent pairs: a part
is $\{a, \mathit{an}\} \, 2^{e_a - \mathit{an} B}$ with a normalised
mantissa, a nonzero magnitude and an `slong` exponent far from overflow
(the callers scale extreme operands). The `_approx` flavour writes
$\mathit{rn}$ limbs per part, each within about an ulp of the magnitude
of the result (the parts truncated), with the exponents of
the parts (`WORD_MIN` for an exactly cancelled part of the quotient);
the outputs are the magnitudes of the parts, the signs being the
caller's unless stated. The reciprocal is $(a - b i) t$ with
$t = 1/(a^2 + b^2)$, written as $a t$ and $b t$. The quotient
$(a + b i)/(c + d i) = (u + v i) t$ with $u = a c + b d$,
$v = b c - a d$ and $t = 1/(c^2 + d^2)$ takes the signs of the operands
and returns those of the parts. The square root is $(q, b/(2q))$ with
$q = \sqrt{(|z| + |a|)/2}$ for $a \ge 0$ and $(|b|/(2q), q)$ for $a < 0$
(the imaginary part has the sign of $b$); the reciprocal square root is
$(T v, b v)$ resp. $(b v, T v)$ with $T = |z| + |a|$ and
$v = 1/\sqrt{2 |z|^2 T}$ (the imaginary part has the sign opposite to
$b$). The outputs may alias the operands' limbs.

::: details
The compositions of `nfloat_complex` on the `_approx` kernels, every
intermediate held at $\mathit{rn} + 1$ limbs (the guard limb keeps the
few ulps of the compositions below the ulp of the output, as the
$2^{2 - \mathit{prec}}$ of `acf`'s callers requires) and each kernel
within an ulp: $|z|^2$ from two approximate squares and a truncated
sum, one reciprocal shared by the two parts of a quotient (a division
costs several products), $|z|$ by the exact root, the larger part of a
root first so that nothing cancels.
:::
:::

::: function
void _arf_mul_tmp_cleanup(void)
void _arf_add_tmp_cleanup(void)

Free the thread-local scratch space of multiplication and addition;
registered with `flint_register_cleanup_function` and called by
`flint_cleanup`.
:::
