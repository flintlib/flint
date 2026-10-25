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
\le |f| < 2^e$ for the exact result $f$. Below a cutoff the functions use
the standard algorithm with `ARF_RND_DOWN`, respectively
`ARF_RND_NEAR` (see `arf_rnd_relaxed_to_strict`). Above it, they use the
Newton iterations of [mp_real](mp_real): inversion from 5000 bits (with
a divisor of more than $\mathit{prec}/2$ bits), division from 12000 bits,
reciprocal square root from 4000 bits and square root from 100000 bits.
The Newton result carries 8 (fast) or 11 (accurate) guard bits and is
rounded to nearest, which gives errors of at most 0.5625 and 0.508 units;
its exponent is that of $f$, or it is exactly $2^e$, so adding one unit in
the last place of the result is a valid error bound. The inexact flag
is not that of the final rounding of the Newton approximation, which
says nothing about $f$, but is determined separately: $1/x$ and
$1/\sqrt{x}$ are dyadic only for powers of two (resp. powers of four),
which are handled without Newton iteration, and an exact quotient or
square root that fits in `prec` bits is found exactly by the iteration,
so it is verified by one exact multiplication, made only when the bit
lengths and the low bits of the mantissas allow an exact result.
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
two-limb operands; for operands between 25 limbs and a few hundred limbs
it calls MPFR, whose short product computes only the top half of the
mantissa.
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
The standard algorithm computes a truncated quotient with at least 32
extra bits (`flint_mpn_tdiv_q`) and decides exactness by multiplying back
only when the low bits of the quotient are zero. Division by a power of
two is a rounding of the numerator. With a relaxed mode and a large
precision, Newton iteration from [mp_real](mp_real) is used instead (see
`ARF_RND_FAST`).
:::
:::

::: function
int arf_sqrt(arf_ptr z, arf_srcptr x, slong prec, arf_rnd_t rnd)
int arf_sqrt_ui(arf_t z, ulong x, slong prec, arf_rnd_t rnd)
int arf_sqrt_fmpz(arf_t z, const fmpz_t x, slong prec, arf_rnd_t rnd)

Set `z` to $\sqrt{x}$, correctly rounded for the strict modes. A negative
`x` gives NaN. The relaxed modes are accepted.
:::

::: function
int arf_rsqrt(arf_ptr z, arf_srcptr x, slong prec, arf_rnd_t rnd)

Sets `z` to $1/\sqrt{x}$, correctly rounded for the strict modes. A
negative `x` gives NaN and $x = 0$ gives $+\infty$. The relaxed modes
are accepted.
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
or NaN; special values follow from this rule. `arf_complex_mul` forms
the exact products of the mantissas and rounds $ac - bd$ and $ad + bc$
once each; when all four mantissas have at least 20 limbs and the
exponents of the real and imaginary parts are within 64 bits of each
other, it uses Karatsuba's three-product formula $ad + bc = (a+b)(c+d) -
ac - bd$. `arf_complex_mul_fallback` is the plain version with four
products, kept for testing.
:::
:::

::: function
int arf_complex_sqr(arf_t e, arf_t f, const arf_t a, const arf_t b, slong prec, arf_rnd_t rnd)

Sets $e + fi = (a + bi)^2$, with the same conventions as
`arf_complex_mul`, but faster.
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
int arf_mul_via_mpfr(arf_t z, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)

The implementations of `arf_mul`: for any rounding mode, for
`ARF_RND_DOWN`, and through `mpfr_mul`.
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
of `y`. `x` may point into the limbs of `y`, as the Newton division code
does.
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
int _arf_want_newton_inv(const arf_t x, slong prec)
int _arf_want_newton_div(const arf_t x, const arf_t y, slong prec)
int _arf_want_newton_sqrt(const arf_t x, slong prec)
int _arf_want_newton_rsqrt(const arf_t x, slong prec)
int _arf_want_newton_inv_large(const arf_t x, slong prec)
int _arf_want_newton_div_large(const arf_t x, const arf_t y, slong prec)

Decide whether a relaxed-mode inversion, division, square root or
reciprocal square root at precision `prec` uses Newton iteration; the
`_large` versions contain the tests that are made only above the
precision cutoffs.
:::

::: function
int _arf_inv_newton(arf_t res, const arf_t x, slong prec, arf_rnd_t rnd)
int _arf_div_newton(arf_t res, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)
int _arf_sqrt_newton(arf_t res, const arf_t x, slong prec, arf_rnd_t rnd)
int _arf_rsqrt_newton(arf_t res, const arf_t x, slong prec, arf_rnd_t rnd)

Compute $1/x$, $x/y$, $\sqrt{x}$ and $1/\sqrt{x}$ for a relaxed mode
`rnd` by Newton iteration, with the error bounds stated for
`ARF_RND_FAST`. Require finite nonzero (positive, for the roots) input.
Return the inexact flag of the result (0 if and only if `res` equals
the exact value; see `ARF_RND_FAST` for how it is determined).
:::

::: function
void _arf_mul_tmp_cleanup(void)
void _arf_add_tmp_cleanup(void)

Free the thread-local scratch space of multiplication and addition;
registered with `flint_register_cleanup_function` and called by
`flint_cleanup`.
:::
