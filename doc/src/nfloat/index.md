---
title: nfloat.h
header: nfloat.h
description: Packed floating-point numbers with precision in whole words
category: Floating-point arithmetic
order: 1
---

This module provides binary floating-point numbers with a precision of
$n$ full words, $1 \le n \le 66$ on 64-bit machines (64 to 4224 bits),
stored in a flat array of $n + 2$ words with no heap allocation, and
complex numbers built from them. Because an element has a fixed size
determined by the context, vectors and matrices of `nfloat` numbers are
packed contiguously, temporaries can live on the stack, and every
operation can switch to a kernel specialised for its number of words.
The module is designed for the generic rings interface: a context
created with `nfloat_ctx_init` is a `gr_ctx_t`, the functions return
status codes, and `gr_poly`, `gr_mat` and the other generic modules
work over it. For approximate polynomial root-finding and linear
algebra this is about twice as fast as the same algorithms over `arf_t`
or MPFR [@Joh2025gr].

The focus is on speed, not on bitwise-specified results. Without
directed rounding, the basic operations are correctly rounded toward
zero or at most a few ulps off, and the elementary functions are
accurate to about one ulp; with the directed rounding flags, most
basic operations are correctly rounded and the other supported
operations return valid bounds (see [](nfloat/rounding)). Results that
are not correctly rounded may differ between 64-bit and 32-bit
machines. For variable precision, precision above 4224 bits, or
rigorous error bounds, use `arf_t` and `arb_t`; for 53 to 212 bits with
SIMD-friendly double arithmetic, use `dfloat`; see
[](algorithms/floating-point/overview) for a comparison.

The guide [](nfloat/getting-started) shows how to use the module. The
fixed-size kernels, the dot product and the matrix multiplication
algorithms are described in [](algorithms/floating-point/nfloat).

## Representation {#representation}

An element with $n$-word precision consists of $n + 2$ words: the
exponent (an `slong`), the sign (a whole word, 0 or 1, to avoid bit
manipulation), and the mantissa $d_0, \ldots, d_{n-1}$, least
significant word first. A normal (nonzero finite) value is
$$x = (-1)^{\mathit{sign}} \cdot 0.d_{n-1} d_{n-2} \ldots d_0 \cdot 2^e,$$
where the mantissa is a fraction in $[1/2, 1)$, so that the top bit of
$d_{n-1}$ is set and $|x| \in [2^{e-1}, 2^e)$. The exponent is limited
to about a quarter of the word range (`NFLOAT_MIN_EXP` to
`NFLOAT_MAX_EXP`) so that sums of exponents cannot overflow. Zero,
infinities and NaN are encoded by special exponent values, with
arbitrary data in the mantissa.

::: macro
NFLOAT_MIN_LIMBS
NFLOAT_MAX_LIMBS

The allowed numbers of words $n$: $1 \le n \le 66$ on 64-bit machines
and $1 \le n \le 132$ on 32-bit machines, so at most 4224 bits. The
upper limit keeps elements and temporary buffers small enough for the
stack and makes copying and swapping cheap compared with a
pointer-and-size representation.
:::

::: macro
NFLOAT_HEADER_LIMBS
NFLOAT_EXP(x)
NFLOAT_SGNBIT(x)
NFLOAT_D(x)
NFLOAT_DATA(x)

The number of header words, 2, and access to the exponent, the sign
word, the mantissa (a pointer to $d_0$) and the whole element as an
array of words.
:::

::: macro
NFLOAT_MAX_ALLOC

The number of words of a temporary that can hold an element of any
precision, `NFLOAT_HEADER_LIMBS + NFLOAT_MAX_LIMBS`.
:::

::: macro
NFLOAT_MIN_EXP
NFLOAT_MAX_EXP

The exponent range of normal values, $\pm$`WORD_MAX / 4`. Results
outside this range underflow or overflow.
:::

::: macro
NFLOAT_EXP_ZERO
NFLOAT_EXP_POS_INF
NFLOAT_EXP_NEG_INF
NFLOAT_EXP_NAN
NFLOAT_IS_SPECIAL(x)
NFLOAT_IS_ZERO(x)
NFLOAT_IS_POS_INF(x)
NFLOAT_IS_NEG_INF(x)
NFLOAT_IS_INF(x)
NFLOAT_IS_NAN(x)

The exponent values encoding zero, $+\infty$, $-\infty$ and NaN (all
below `NFLOAT_MIN_EXP`), and tests for them. `NFLOAT_IS_SPECIAL` is
true for zero, infinities and NaN.
:::

::: macro
LIMB_MSB_IS_SET(n)

Whether the most significant bit of the word $n$ is set.
:::

## Types {#types}

::: type
nfloat_ptr
nfloat_srcptr

Pointers to an element or a vector of elements of any precision, as
`void *` and `const void *`. Pointer arithmetic needs the element size
from the context, for example through `GR_ENTRY`.
:::

::: type
nfloat64_struct
nfloat128_struct
nfloat192_struct
nfloat256_struct
nfloat384_struct
nfloat512_struct
nfloat1024_struct
nfloat2048_struct
nfloat4096_struct
nfloat64_t
nfloat128_t
nfloat192_t
nfloat256_t
nfloat384_t
nfloat512_t
nfloat1024_t
nfloat2048_t
nfloat4096_t

Structs of the right size for elements of some common precisions, for
declaring elements on the stack; `nfloatX_t` is a length-one array of
`nfloatX_struct`, so it is passed by reference. The precision of the
context must match the type.
:::

::: type
_nfloat_ctx_struct

The context data: the number of words `nlimbs`, the `flags`, an unused
rounding field `rnd`, and the function precision `func_prec`.
:::

::: macro
NFLOAT_CTX(ctx)
NFLOAT_CTX_NLIMBS(ctx)
NFLOAT_CTX_PREC(ctx)
NFLOAT_CTX_FLAGS(ctx)
NFLOAT_CTX_RND(ctx)
NFLOAT_CTX_FUNC_PREC(ctx)
NFLOAT_CTX_DATA_NLIMBS(ctx)
NFLOAT_CTX_HAS_INF_NAN(ctx)
NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx)

Access to the context data: the number of words $n$, the precision
$n \cdot$ `FLINT_BITS`, the flags, the function precision, the element
size $n + 2$ in words, and whether infinities or NaN, or directed
rounding, are enabled.
:::

## Contexts {#contexts}

::: function
int nfloat_ctx_init(gr_ctx_t ctx, slong prec, int flags)

Initialises `ctx` to the ring of floating-point numbers with precision
`prec` bits rounded up to a whole number of words (for example,
`prec = 53` gives 64-bit precision), with the given `flags` (a bitwise
or of the flags below, or 0). Returns `GR_UNABLE` without initialising
`ctx` if `prec` is not positive or exceeds the maximum
`NFLOAT_MAX_LIMBS * FLINT_BITS`, and `GR_SUCCESS` otherwise. The
context needs no memory and is cleared with `gr_ctx_clear`.
:::

::: macro
NFLOAT_ALLOW_UNDERFLOW

By default, a result that underflows the exponent range is set to zero
and the operation returns `GR_UNABLE`. With this flag, it returns
`GR_SUCCESS` with the zero result.
:::

::: macro
NFLOAT_ALLOW_INF

Allows infinities. By default, a result that overflows is set to zero
(a garbage value) and the operation returns `GR_UNABLE`, and
`nfloat_pos_inf` returns `GR_DOMAIN`; with this flag, overflow gives a
signed infinity and `GR_SUCCESS`.
:::

::: macro
NFLOAT_ALLOW_NAN

Allows NaN. By default, an operation without a meaningful result (such
as $0/0$ or $\log(-1)$) returns `GR_UNABLE` or `GR_DOMAIN`; with this
flag, it outputs NaN and returns `GR_SUCCESS`.
:::

::: macro
NFLOAT_RND_FLOOR
NFLOAT_RND_CEIL

Enable directed rounding toward $-\infty$ or $+\infty$ (experimental).
At most one of them may be set. Without either, results are rounded in
an unspecified direction, usually toward zero. Directed rounding may
not work correctly together with `NFLOAT_ALLOW_UNDERFLOW`. Which
operations respect it is listed in [](nfloat/rounding).
:::

Infinities and NaN are disabled by default because this lets many
functions skip checks for them.

::: function
int nfloat_ctx_set_func_prec(gr_ctx_t ctx, slong prec)
slong nfloat_ctx_get_func_prec(gr_ctx_t ctx)

Set and get the *function precision* of `ctx`: the precision in bits
to which the elementary functions (exponentials, logarithms,
trigonometric and hyperbolic functions, their inverses and powers) are
computed. It defaults to the full precision of the context; a value
larger than that is clamped to it, and a value less than 1 makes
`nfloat_ctx_set_func_prec` return `GR_UNABLE`. A smaller function
precision makes these functions cheaper (for example, at 64 bits with a
function precision of 50 bits, most functions evaluate their kernels
with one word instead of two); the arithmetic is unaffected.
:::

::: function
int nfloat_ctx_write(gr_stream_t out, gr_ctx_t ctx)

Writes a description of `ctx` to `out`, such as
`Floating-point numbers with prec = 128 (nfloat)`.
:::

::: function
int nfloat_ctx_set_real_prec(gr_ctx_t ctx, slong prec)
int nfloat_ctx_get_real_prec(slong * res, gr_ctx_t ctx)

The precision of a context cannot be changed:
`nfloat_ctx_set_real_prec` does nothing and returns `GR_UNABLE`.
`nfloat_ctx_get_real_prec` sets `res` to the precision in bits and
returns `GR_SUCCESS`.
:::

## Memory management and constants {#memory}

The functions in this and the following sections implement the methods
of the `gr` interface and can be used interchangeably with their `gr`
counterparts (`gr_add` and so on), which dispatch to them.

::: function
void nfloat_init(nfloat_ptr res, gr_ctx_t ctx)
void nfloat_clear(nfloat_ptr res, gr_ctx_t ctx)

Initialise `res` to zero, and clear it. Since elements own no memory,
clearing does nothing.
:::

::: function
void nfloat_swap(nfloat_ptr x, nfloat_ptr y, gr_ctx_t ctx)
int nfloat_set(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

Swap two elements word by word, and copy `x` to `res`.
:::

::: function
int nfloat_zero(nfloat_ptr res, gr_ctx_t ctx)
int nfloat_one(nfloat_ptr res, gr_ctx_t ctx)
int nfloat_neg_one(nfloat_ptr res, gr_ctx_t ctx)

Set `res` to 0, 1 or $-1$.
:::

::: function
int nfloat_pos_inf(nfloat_ptr res, gr_ctx_t ctx)
int nfloat_neg_inf(nfloat_ptr res, gr_ctx_t ctx)
int nfloat_nan(nfloat_ptr res, gr_ctx_t ctx)

Set `res` to $+\infty$, $-\infty$ or NaN if the context allows it.
Otherwise the infinities return `GR_DOMAIN` and NaN returns
`GR_UNABLE`, leaving `res` unchanged.
:::

::: function
int nfloat_pi(nfloat_ptr res, gr_ctx_t ctx)

Sets `res` to $\pi$ rounded to the precision of the context.
:::

## Conversions {#conversions}

::: function
int nfloat_set_ui(nfloat_ptr res, ulong x, gr_ctx_t ctx)
int nfloat_set_si(nfloat_ptr res, slong x, gr_ctx_t ctx)
int nfloat_set_fmpz(nfloat_ptr res, const fmpz_t x, gr_ctx_t ctx)
int nfloat_set_fmpq(nfloat_ptr res, const fmpq_t v, gr_ctx_t ctx)
int nfloat_set_d(nfloat_ptr res, double x, gr_ctx_t ctx)

Set `res` to an integer, a rational number or a double, rounded to the
precision of the context (toward zero, or in the direction of the
context's directed rounding).
:::

::: function
int nfloat_set_arf(nfloat_ptr res, const arf_t x, gr_ctx_t ctx)
int nfloat_get_arf(arf_t res, nfloat_srcptr x, gr_ctx_t ctx)

Convert from and to `arf_t`; the conversion to `arf_t` is exact. These
are declared only if `arf.h` is included before `nfloat.h`.
:::

::: function
int nfloat_get_fmpz(fmpz_t res, nfloat_srcptr x, gr_ctx_t ctx)

Sets `res` to $x$ if $x$ is an integer and returns `GR_SUCCESS`;
returns `GR_DOMAIN` if $x$ is not an integer and `GR_UNABLE` if it is
infinite or NaN.
:::

::: function
int nfloat_get_d_2exp_si(double * mant, slong * expo, nfloat_srcptr x, gr_ctx_t ctx)

Sets `mant` and `expo` such that $x \approx \mathit{mant} \cdot 2^{\mathit{expo}}$
with $|\mathit{mant}| \in [1/2, 1)$ (or both zero for $x = 0$), from the
leading word of the mantissa, so with a relative error of about
$2^{-53}$. This avoids the overflow of a direct conversion to `double`.
Returns `GR_UNABLE` for infinities and NaN.
:::

::: function
int nfloat_set_mpn_2exp(nfloat_ptr res, nn_srcptr x, slong xn, slong exp, int xsgnbit, gr_ctx_t ctx)
int _nfloat_set_mpn_2exp(nfloat_ptr res, nn_srcptr x, slong xn, slong exp, int xsgnbit, gr_ctx_t ctx)

Set `res` to $(-1)^{\mathit{xsgnbit}} \cdot f \cdot 2^{\mathit{exp}}$,
where $f = \{x, \mathit{xn}\} / 2^{\mathit{xn} \cdot \mathtt{FLINT\_BITS}}$
is the fraction in $[0, 1)$ given by the `xn` words of `x`, rounded to
the precision of the context. The version with an underscore requires
that the top word `x[xn - 1]` is nonzero; the other removes leading zero
words and handles $f = 0$.
:::

::: function
int nfloat_1_set_2_2exp(nfloat_ptr res, ulong x1, ulong x0, slong exp, int xsgnbit, gr_ctx_t ctx)
int nfloat_1_set_3_2exp(nfloat_ptr res, ulong x2, ulong x1, ulong x0, slong exp, int xsgnbit, gr_ctx_t ctx)
int nfloat_2_set_3_2exp(nfloat_ptr res, ulong x2, ulong x1, ulong x0, slong exp, int xsgnbit, gr_ctx_t ctx)
int nfloat_2_set_4_2exp(nfloat_ptr res, ulong x3, ulong x2, ulong x1, ulong x0, slong exp, int xsgnbit, gr_ctx_t ctx)

Inline versions of `nfloat_set_mpn_2exp` for a two-, three- or
four-word fraction (most significant word first in the argument list)
and a context of exactly one word (`nfloat_1_`) or two words
(`nfloat_2_`). The fraction is normalised and truncated (not rounded)
to the precision; zero, underflow and overflow are handled as by the
other functions.
:::

::: function
int nfloat_set_str(nfloat_ptr res, const char * x, gr_ctx_t ctx)

Parses a decimal number with `arb_set_str` and rounds it to nearest
(this ignores directed rounding); other expressions are parsed with
the generic string parser of the `gr` interface.
:::

::: function
int nfloat_set_other(nfloat_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)

Sets `res` to the element `x` of the ring `x_ctx`, rounded. Elements
of other `nfloat` contexts, `fmpz`, `fmpq`, `arf_t` and `arb_t` (whose
midpoint is used) are converted directly; other rings are converted
through `arf_t` with `gr_set_other`, whose status is returned (for
example `GR_DOMAIN` for a non-real number).
:::

::: function
int nfloat_write(gr_stream_t out, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_randtest(nfloat_ptr res, flint_rand_t state, gr_ctx_t ctx)

Write $x$ in decimal as the `arf_t` ring of the same precision does,
and set `res` to a random element (from `arf_randtest`) for testing.
:::

## Comparison {#comparison}

::: function
truth_t nfloat_equal(nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)

Returns `T_TRUE` if $x$ and $y$ have the same representation, which for
normal values means equal values, and `T_FALSE` otherwise. Unlike IEEE
comparison, a NaN equals a NaN.
:::

::: function
truth_t nfloat_is_zero(nfloat_srcptr x, gr_ctx_t ctx)
truth_t nfloat_is_one(nfloat_srcptr x, gr_ctx_t ctx)
truth_t nfloat_is_neg_one(nfloat_srcptr x, gr_ctx_t ctx)

Return whether $x$ is 0, 1 or $-1$; `T_UNKNOWN` for NaN.
:::

::: function
int nfloat_cmp(int * res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
int nfloat_cmpabs(int * res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)

Set `res` to $-1$, 0 or 1 according to the comparison of $x$ and $y$,
or of $|x|$ and $|y|$, and return `GR_SUCCESS`; return `GR_UNABLE` if
an operand is NaN.
:::

::: function
int nfloat_sgn(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_im(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

Set `res` to the sign of $x$ ($-1$, 0 or 1, and NaN for NaN), and to
the imaginary part of $x$, which is zero.
:::

## Arithmetic {#arithmetic}

::: function
int nfloat_neg(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_abs(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

Set `res` to $-x$ or $|x|$, exactly.
:::

::: function
int nfloat_add(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
int nfloat_sub(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
int nfloat_mul(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
int nfloat_sqr(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

Sum, difference, product and square. Without directed rounding the
error is less than 2 ulp for the sum and difference (the smaller
operand is truncated when it is aligned) and at most 1 ulp for the
products (truncation, except that the approximate high product used
from three words on can rarely give one ulp less in magnitude); with
directed rounding the results are correctly rounded.

::: details
The functions dispatch on the number of words to fixed-size kernels:
for one and two words the operations are done in registers, the
addition through `_nfloat_add_1`, `_nfloat_add_2`, ..., `_nfloat_add_4`
and their subtraction counterparts, and the product with a full
$1 \times 1$ or $2 \times 2$ word product. From three words the product
uses the high product `flint_mpn_mulhigh_n`, whose error is at most a
few units of the last word. See [](algorithms/floating-point/nfloat).
:::
:::

::: function
int nfloat_addmul(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
int nfloat_submul(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)

Set `res` to $\mathit{res} + xy$ or $\mathit{res} - xy$, with two
roundings, so the result can be far from the exact value (measured in
ulp of the result) when the terms cancel.
:::

::: function
int nfloat_mul_2exp_si(nfloat_ptr res, nfloat_srcptr x, slong y, gr_ctx_t ctx)

Sets `res` to $x \cdot 2^y$, exactly unless the result overflows or
underflows.
:::

::: function
int nfloat_inv(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_div(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
int nfloat_div_ui(nfloat_ptr res, nfloat_srcptr x, ulong y, gr_ctx_t ctx)
int nfloat_div_si(nfloat_ptr res, nfloat_srcptr x, slong y, gr_ctx_t ctx)

Inverse and quotients. Without directed rounding, the result is the
quotient truncated toward zero or, for divisors of three or more
nonzero words, possibly one ulp larger in magnitude (an error of less
than 1 ulp); `div_ui` and `div_si` are always truncated. With directed
rounding, the result is the correctly rounded floor or ceiling of the
exact quotient. Division by zero gives NaN, so it returns `GR_UNABLE`
unless NaN is allowed.

::: details
For one and two words, the mantissas are divided by an inlined $2/1$
division or two $3/2$ divisions, with the dividend shifted right by one
bit beforehand when it is larger than the divisor so that the quotient
comes out normalised; the remainder gives the exact rounding
information. These use hardware division where it is fast
(`FLINT_PREINVERT_LIMB_USE_NATIVE`) and a precomputed inverse otherwise
[@GraMol2010]. Divisors with a single nonzero word (such as small
integers) use a chain of $2/1$ divisions, $O(n)$. Otherwise a
schoolbook approximate division (`_flint_mpn_divapprox_basecase_preinv1`)
computes the quotient or the quotient plus one at the target precision.
With directed rounding, one guard word is computed and the truncation is
certified unless the guard word is 0 or 1, in which case (for example
for exact quotients) the exact quotient and remainder are computed.
:::
:::

::: function
int nfloat_sqrt(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_rsqrt(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

Square root and reciprocal square root. The square root is correctly
rounded in every mode (truncated toward zero by default); the
reciprocal square root has an error of less than $1 + 2^{-37}$ ulp by
default and is correctly rounded with directed rounding. Negative
arguments give NaN.

::: details
The square root is the truncated integer square root of the mantissa
padded to $2n$ words, with an exactness flag (`flint_mpn_sqrtrem`). At
one word on 64-bit machines it is computed in registers: a
double-precision root (within $2^{12}$ of the exact root), one
correction from the exact residual evaluated in double precision, and a
final exact residual, which gives the floor and the exactness flag. At
two words the root is first approximated at three words from the
two-word reciprocal square root below and one Newton step with the
exact residual, which certifies the truncation unless the guard word is
within 3 units of 0 or $2^{64}$.

The reciprocal square root refines a double-precision approximation by
one Newton step at 128 bits (one word), or a two-word starting value
by one Newton step with the residual evaluated exactly from the
six-word product of the mantissa and the square of the starting value
(two words). For $n \ge 3$ words, a fixed-point Newton iteration with
an added third-order term yields the result with one guard word and an
error of about one unit in the guard word. With directed rounding the
truncation is certified from these error bounds, and in the rare
uncertain cases the exact result
$\lfloor \sqrt{\lfloor z \rfloor} \rfloor$ is computed for
$z = 2^{3nw}/(4a)$ or $2^{3nw}/(2a)$, in terms of the $n$-word mantissa
$a$ and the word size $w$.
:::
:::

::: function
int nfloat_hypot(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)

Sets `res` to $\sqrt{x^2 + y^2}$, correctly rounded in every mode
(truncated by default), without intermediate overflow or underflow.
This also implements `nfloat_complex_abs`.

::: details
With $|x| \ge |y|$ and $d$ the difference of the exponents, the result
is $|x|(1 + \delta)$ with $0 < \delta < 2^{-p-3}$ ($p$ the precision)
when $d$ exceeds about half the precision. Otherwise
$x^2 + y^2 = S \cdot 2^{2e}$ for an integer $S$ of at most $3n + 2$
words (the squares of the mantissas, aligned exactly), whose truncated
square root with an inexactness flag gives the result. At one word on
64-bit machines, $S$ is formed in registers and the truncated root is
computed from the integer square root of its top two words, since the
discarded low bits of $S$ do not change the floor of the square root. At
two words, $S$ is formed in registers (the smaller square truncated to
the four words of the larger) and its root is approximated with a guard
word, which certifies the truncation in all but rare cases. From 12
words (or from 5 words when the exponents differ by at least $3/8$ of
the precision), the exact squares are replaced by high squarings
(`flint_mpn_sqrhigh`), whose sum is a lower bound within about $4n$
units of its last word, and the square root with a guard word
certifies the truncation likewise, with the exact computation as the
fallback.
:::
:::

::: function
int nfloat_floor(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_ceil(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_trunc(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_nint(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

Round $x$ to an integer toward $-\infty$, toward $+\infty$, toward zero,
or to the nearest integer (ties to even). The results are exact.
:::

## Elementary functions {#elementary}

The elementary functions are computed natively, on the mantissas, with
the fixed-point kernels of the `mp_real` module (table-based argument
reduction and rectangular splitting of the Taylor series [@Joh2014c];
see [](algorithms/floating-point/elementary-functions)); `nfloat`
precisions never need bit-burst evaluation. The
results are not correctly rounded: without directed rounding the error
is at most about one ulp at the function precision of the context (see
`nfloat_ctx_set_func_prec`), and with directed rounding the result is a
valid lower or upper bound within a few ulps. Exact special cases give
exact results (for example $\exp(0)$, $\log 1$, $\sin(\pi x)$ for
multiples of $1/2$, $\log_2$ of powers of two). Results that overflow or
underflow follow the flags of the context, and arguments outside the
real domain (such as $\log$ of a negative number or
$\operatorname{asin}$ of $|x| > 1$) give NaN, so they return
`GR_UNABLE` or `GR_DOMAIN` unless NaN is allowed.

::: details
On 64-bit machines, `nfloat_exp`, `nfloat_expm1`, `nfloat_log`,
`nfloat_log1p`, `nfloat_pow`, `nfloat_sin_cos`, `nfloat_sin`,
`nfloat_cos`, `nfloat_tan`, `nfloat_cot`, their versions at $\pi x$,
`nfloat_sinh_cosh`, `nfloat_sinh`, `nfloat_cosh`, `nfloat_tanh`,
`nfloat_coth`, `nfloat_atan` and `nfloat_atan2` have in-register fast
paths at one and two words for common arguments (for example every
argument below the overflow threshold for the exponential, and
$2^{-33} \le |x| < 2^{63}$ for the sine and cosine at one word): the
argument reduction on the mantissa, the table-driven kernels of
`mp_real` at two or three words (shared with its 128- and 192-bit
functions), and the final normalisation, without temporary arrays. The
hyperbolic functions use a paired kernel for $e^t$ and $e^{-t}$ (and
one division for `tanh` and `coth`), the tangent and cotangent divide
the sine and the cosine in registers, `nfloat_pow` multiplies a
logarithm with one more word (or two, for large $|y|$) by $y$ exactly
before the exponential, and the arctangent of $|x| > 1$ and
`nfloat_atan2` evaluate the arctangent of a ratio of mantissas with a
single division. Arguments where the reduced argument would lose
relative accuracy (close to 1 for the logarithm, close to a multiple of
$\pi/2$ for the sine, cosine and tangent, very small arguments of
`expm1`, `log1p` and `sinh`) and rarer ranges take the general path.

From three to nine words the general path runs on the dense-table
kernels of `mp_real` at four to ten words, and from ten to thirty-one
words on the dense-table kernels with tables computed at run time
(shared by all threads and built on first use).

Otherwise the functions wrap the ball functions of `mp_real` (for
example `mp_real_sin_cos_bits`), evaluated on the exact argument at
increasing precision until the ball determines the result. This handles
trigonometric functions of huge arguments (exponent above 65536, up to
$2^{22}$, beyond which `GR_UNABLE` is returned), the less common
functions for arguments whose intermediate results overflow or
underflow, and the maximum precision `NFLOAT_MAX_LIMBS`, which leaves
no room for a guard word.
:::

::: function
int nfloat_exp(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_expm1(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_exp2(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_log(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_log1p(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_log2(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

Exponentials $e^x$, $e^x - 1$, $2^x$ and logarithms $\log x$,
$\log(1 + x)$, $\log_2 x$. The exponentials reduce $|x|$ modulo
$\log 2$ in fixed point (with one integral word); `expm1` and `log1p`
have relative accuracy near zero, and `exp2` and `log2` are exact for
integers and powers of two respectively.
:::

::: function
int nfloat_sin(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_cos(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_sin_cos(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_tan(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_sin_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_cos_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_sin_cos_pi(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_tan_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

Sine, cosine and tangent of $x$ and of $\pi x$; the `sin_cos` versions
compute both. The argument is reduced modulo $\pi/2$ in fixed point,
repeated with more bits when the reduced argument is close to a zero of
the result, so that the results have relative accuracy; the functions
of $\pi x$ reduce $x$ modulo 2 exactly.
:::

::: function
int nfloat_atan(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_atan2(nfloat_ptr res, nfloat_srcptr y, nfloat_srcptr x, gr_ctx_t ctx)

The arctangent, and the two-argument arctangent
$\operatorname{atan2}(y, x) \in (-\pi, \pi]$. For $|x| > 1$, `atan` uses
$\pi/2 - \operatorname{atan}(1/|x|)$; `atan2` takes the arctangent of
the ratio of the smaller and the larger magnitude and adds the octant.
$\operatorname{atan2}(0, 0) = 0$.
:::

::: function
int nfloat_sinh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_cosh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_sinh_cosh(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_tanh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

Hyperbolic functions, from $E = e^{|x|}$ and $1/E$, with the difference
$E - 1/E$ computed in fixed point for $|x| < 1$.
:::

::: function
int nfloat_pow(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)

Sets `res` to $x^y = \exp(y \log|x|)$, with the sign $(-1)^y$ for
$x < 0$ and integer $y$ (NaN for $x < 0$ and noninteger $y$). Simple
cases ($y = \pm 1, 2, \pm 1/2$, $|x| = 1$, integer powers of powers of
two) are handled directly.
:::

::: function
int nfloat_exp10(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_log10(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_cot(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_sec(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_csc(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_sinc(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_cot_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_sec_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_csc_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_sinc_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_coth(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_sech(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_csch(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_asin(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_acos(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_asin_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_acos_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_atan_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_acot(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_asec(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_acsc(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_acot_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_asec_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_acsc_pi(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_asinh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_acosh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_atanh(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_acoth(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_asech(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_acsch(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

The less common elementary functions. After the exact special values
(zeros, poles, $|x| = 1$ for the inverse functions, multiples of $1/2$
for the functions of $\pi x$), they are evaluated as short compositions
of the functions above in an internal context with a guard word (none
when the function precision is sufficiently below the precision), using
formulas without cancellation, such as
$\operatorname{asin}(x) = \operatorname{atan2}(x, \sqrt{(1 - x)(1 + x)})$
and $\operatorname{asinh}(x) = \operatorname{log1p}(|x| + x^2/(1 + \sqrt{1 + x^2}))$,
and a rigorous bound for the composed error, which is added to obtain
valid bounds with directed rounding.

::: details
On 64-bit machines, at one and two words, the reciprocal functions
(`sec`, `csc`, `sech`, `csch`, also of $\pi x$) and `sinc` (also of
$\pi x$) instead divide the fixed-point outputs of the in-register fast
paths of the sine, the cosine and the paired exponentials, and the
inverse functions are evaluated in fixed point at one more word:
$y^2 \pm 1$ exactly, square roots by a refined reciprocal square root
(as for `nfloat_sqrt` at two words), the arctangent of a ratio (`asin`,
`acos`, `asec`, `acsc`, `acot` and the functions divided by $\pi$, which
multiply by $2/\pi$) or the logarithm kernels (the inverse hyperbolic
functions, as sums of logarithms of the same sign, for example
$\operatorname{atanh}(y) = (\operatorname{log1p}(y) - \log(1 - y))/2$).
Arguments where an intermediate result would lose relative accuracy at
that precision (results below about $2^{-32}$ or $2^{-40}$ from small
differences, $|x| \ge 2^{62}$) take the general path.
:::
:::

::: function
int nfloat_gamma(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
int nfloat_zeta(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)

The gamma function and the Riemann zeta function, computed by
converting to `arf_t` and calling the generic ball implementation
(`gr_gamma`, `gr_zeta` over `arb_t`), then rounding.
:::

## Vector functions {#vectors}

These override the generic `gr` vector operations with inlined or
partially inlined code, which avoids the per-element dispatch.

::: function
void _nfloat_vec_init(nfloat_ptr res, slong len, gr_ctx_t ctx)
void _nfloat_vec_clear(nfloat_ptr res, slong len, gr_ctx_t ctx)
int _nfloat_vec_set(nfloat_ptr res, nfloat_srcptr x, slong len, gr_ctx_t ctx)
int _nfloat_vec_zero(nfloat_ptr res, slong len, gr_ctx_t ctx)

Initialise a vector of `len` elements to zero, clear it (no-op), copy
and zero.
:::

::: function
int _nfloat_vec_add(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx)
int _nfloat_vec_sub(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx)
int _nfloat_vec_mul(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx)
int _nfloat_vec_mul_scalar(nfloat_ptr res, nfloat_srcptr x, slong len, nfloat_srcptr y, gr_ctx_t ctx)
int _nfloat_vec_addmul_scalar(nfloat_ptr res, nfloat_srcptr x, slong len, nfloat_srcptr y, gr_ctx_t ctx)
int _nfloat_vec_submul_scalar(nfloat_ptr res, nfloat_srcptr x, slong len, nfloat_srcptr y, gr_ctx_t ctx)

Elementwise sum, difference and product, product with the scalar $y$,
and $\mathit{res}_i \pm x_i y$, with the same results as the scalar
functions.
:::

::: function
int _nfloat_vec_div(nfloat_ptr res, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx)
int _nfloat_vec_div_scalar(nfloat_ptr res, nfloat_srcptr x, slong len, nfloat_srcptr c, gr_ctx_t ctx)
int _nfloat_vec_div_scalar_ui(nfloat_ptr res, nfloat_srcptr x, slong len, ulong c, gr_ctx_t ctx)
int _nfloat_vec_div_scalar_si(nfloat_ptr res, nfloat_srcptr x, slong len, slong c, gr_ctx_t ctx)

Elementwise quotients and quotients by a scalar, with the same results
as the scalar division, except that `_nfloat_vec_div_scalar` for three
or more words and a divisor of two or more nonzero words may differ by
up to 1 ulp without directed rounding.

::: details
The elementwise division inlines the one- and two-word kernels.
Division by a scalar precomputes the one-word inverse, the $3/2$
inverse, or for three or more words an approximate reciprocal with one
guard word, so that each quotient costs one high product of $n + 1$
words; with directed rounding, the truncation is certified using the
guard word and otherwise recomputed by exact division.
:::
:::

::: function
int _nfloat_vec_dot(nfloat_ptr res, nfloat_srcptr initial, int subtract, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx)
int _nfloat_vec_dot_rev(nfloat_ptr res, nfloat_srcptr initial, int subtract, nfloat_srcptr x, nfloat_srcptr y, slong len, gr_ctx_t ctx)

Set `res` to $s \pm \sum_{i} x_i y_i$ (or with $y$ reversed), where $s$
is `initial` (zero if `NULL`) and the sign is negative if `subtract` is
nonzero. The sum is accumulated in fixed point with one guard word
after a pass that finds the largest exponent, so the result has an
error of a few ulps relative to the largest term rather than one
rounding per term; with directed rounding it is a valid bound. See
[](algorithms/floating-point/nfloat#dot).
:::

## Matrix functions {#matrices}

These implement the `gr_mat` methods of the `nfloat` contexts; the
generic algorithms (solving, inversion, determinants, ...) use them
through the method table.

::: function
int nfloat_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)

Sets $C = AB$, choosing between classical multiplication with the dot
product, `nfloat_mat_mul_fixed` and `nfloat_mat_mul_block` according
to the dimensions and the precision (tuned cutoffs). Contexts with
infinities, NaN or directed rounding, and matrices with a dimension of
at most 2, always use classical multiplication. The methods are
described in [](algorithms/floating-point/nfloat#matmul).
:::

::: function
int nfloat_mat_mul_fixed(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, slong max_extra_prec, gr_ctx_t ctx)

Sets $C = AB$ by converting $A$ and $B$ to fixed-point matrices with
common exponents and multiplying them with `_nfixed_mat_mul`. The
fixed-point precision is the target precision plus the spread of the
exponents in $A$ and in $B$ plus a few guard bits; returns `GR_UNABLE`
if this exceeds the target by `max_extra_prec` bits or more, if a
matrix has zero entries (but is not zero), or if infinities or NaN are
enabled.
:::

::: function
int nfloat_mat_mul_block(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, slong min_block_size, gr_ctx_t ctx)

Sets $C = AB$ by splitting the inner dimension into blocks in which
the entries of each row of $A$ and each column of $B$ have similar
exponents, scaling the blocks to integer matrices and multiplying them
exactly with `fmpz_mat_mul` (multimodular or FFT methods), as
`arb_mat_mul_block` [@Joh2019dot]. Blocks shorter than
`min_block_size` are handled with dot products.
:::

::: function
int nfloat_mat_nonsingular_solve_tril(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
int nfloat_mat_nonsingular_solve_triu(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
int nfloat_mat_lu(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)
int nfloat_mat_lq(gr_mat_t L, gr_mat_t Q, const gr_mat_t A, gr_ctx_t ctx)

Triangular solving, LU factorisation and LQ factorisation, choosing
between the classical (Gram–Schmidt for LQ) and the recursive versions
of the generic algorithms (`gr_mat_nonsingular_solve_tril_classical`,
`gr_mat_lu_recursive`, `gr_mat_lq_recursive`, ...) by precision-dependent
cutoffs, so that large problems run on fast matrix multiplication.
:::

## Complex numbers {#complex}

Complex numbers are pairs of `nfloat` numbers, the real part followed by
the imaginary part, so an element has $2(n + 2)$ words. The context is
created with `nfloat_complex_ctx_init` and has the same data as the
real context, so that the real functions can be applied to the parts.
The arithmetic and elementary functions are approximate: each part
typically has an error of a few ulps relative to the magnitude of the
result, and directed rounding is not supported except by
`nfloat_complex_abs`.

::: type
nfloat_complex_ptr
nfloat_complex_srcptr

Pointers to complex elements (the same as `nfloat_ptr` and
`nfloat_srcptr`).
:::

::: function
int nfloat_complex_ctx_init(gr_ctx_t ctx, slong prec, int flags)

Initialises `ctx` to the ring of complex numbers whose parts have
precision `prec` bits rounded up to a whole number of words, with the
`flags` of `nfloat_ctx_init`. Returns `GR_UNABLE` if the precision is
out of range.
:::

::: macro
NFLOAT_COMPLEX_CTX_DATA_NLIMBS(ctx)
NFLOAT_COMPLEX_RE(ptr, ctx)
NFLOAT_COMPLEX_IM(ptr, ctx)
NFLOAT_COMPLEX_IS_SPECIAL(x, ctx)
NFLOAT_COMPLEX_IS_ZERO(x, ctx)

The size of a complex element in words, pointers to the real and
imaginary parts, and tests whether a part is special or both parts are
zero.
:::

::: function
void nfloat_complex_init(nfloat_complex_ptr res, gr_ctx_t ctx)
void nfloat_complex_clear(nfloat_complex_ptr res, gr_ctx_t ctx)
void nfloat_complex_swap(nfloat_complex_ptr x, nfloat_complex_ptr y, gr_ctx_t ctx)
int nfloat_complex_set(nfloat_complex_ptr res, nfloat_complex_ptr x, gr_ctx_t ctx)
int nfloat_complex_zero(nfloat_complex_ptr res, gr_ctx_t ctx)
int nfloat_complex_one(nfloat_complex_ptr res, gr_ctx_t ctx)
int nfloat_complex_neg_one(nfloat_complex_ptr res, gr_ctx_t ctx)
int nfloat_complex_i(nfloat_complex_ptr res, gr_ctx_t ctx)
int nfloat_complex_pi(nfloat_complex_ptr res, gr_ctx_t ctx)

Memory management and constants, as for the real numbers; `i` sets
`res` to the imaginary unit.
:::

::: function
int nfloat_complex_get_acf(acf_t res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_set_acf(nfloat_complex_ptr res, const acf_t x, gr_ctx_t ctx)
int nfloat_complex_get_acb(acb_t res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_set_acb(nfloat_complex_ptr res, const acb_t x, gr_ctx_t ctx)

Conversions from and to `acf_t` (exact for `get`) and `acb_t` (the
midpoint for `set`). They are declared only if `acf.h` or `acb.h` is
included before `nfloat.h`.
:::

::: function
int nfloat_complex_set_si(nfloat_complex_ptr res, slong x, gr_ctx_t ctx)
int nfloat_complex_set_ui(nfloat_complex_ptr res, ulong x, gr_ctx_t ctx)
int nfloat_complex_set_fmpz(nfloat_complex_ptr res, const fmpz_t x, gr_ctx_t ctx)
int nfloat_complex_set_fmpq(nfloat_complex_ptr res, const fmpq_t x, gr_ctx_t ctx)
int nfloat_complex_set_d(nfloat_complex_ptr res, double x, gr_ctx_t ctx)
int nfloat_complex_set_other(nfloat_complex_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
int nfloat_complex_write(gr_stream_t out, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_randtest(nfloat_complex_ptr res, flint_rand_t state, gr_ctx_t ctx)

Conversions from real values (with zero imaginary part) and from other
rings, output, and random elements.
:::

::: function
truth_t nfloat_complex_is_zero(nfloat_complex_srcptr x, gr_ctx_t ctx)
truth_t nfloat_complex_is_one(nfloat_complex_srcptr x, gr_ctx_t ctx)
truth_t nfloat_complex_is_neg_one(nfloat_complex_srcptr x, gr_ctx_t ctx)
truth_t nfloat_complex_equal(nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)
int nfloat_complex_cmp(int * res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)
int nfloat_complex_cmpabs(int * res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)

Predicates and comparisons. `cmp` compares real numbers and returns
`GR_DOMAIN` if an imaginary part is nonzero; `cmpabs` compares absolute
values.
:::

::: function
int nfloat_complex_conj(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_re(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_im(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_neg(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_add(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)
int nfloat_complex_sub(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)
int nfloat_complex_mul(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)
int nfloat_complex_sqr(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_mul_2exp_si(nfloat_complex_ptr res, nfloat_complex_srcptr x, slong y, gr_ctx_t ctx)
int nfloat_complex_abs(nfloat_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)

Conjugate, real and imaginary parts (as complex numbers), negation,
sum, difference, product, square, multiplication by $2^y$ and absolute
value (a real number, correctly rounded as by `nfloat_hypot`).
:::

::: function
int nfloat_complex_inv(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_div(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)
int nfloat_complex_sqrt(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_rsqrt(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)

Inverse, quotient, square root and reciprocal square root (principal
branches), with errors of a few ulps relative to the magnitude of the
result. Purely real or imaginary operands are handled with the real
functions.

::: details
The inverse and the quotient are $\bar{y}/|y|^2$ and $x\bar{y}/|y|^2$
with a single reciprocal of $|y|^2$ shared by both parts (operands whose
parts differ greatly in magnitude, or have extreme exponents, are
rescaled first). At one word, $|y|^2$ is formed in registers and
truncated to one word, and both parts of the numerator are divided by
it with one shared inverse.

The square root of $z = a + bi$ is computed from $r = |z|$ (by
`nfloat_hypot` at one and two words and for extreme exponents,
otherwise as $\sqrt{a^2 + b^2}$) as $q = \sqrt{(r + |a|)/2}$ and
$b/(2q)$, assigned to the real and imaginary parts according to the
sign of $a$, so that nothing cancels; the reciprocal square root uses
$v = 1/\sqrt{2r^2(r + |a|)}$ and the products $(r + |a|)v$ and $bv$.
At one word on 64-bit machines, both are computed in registers.
:::
:::

::: function
int nfloat_complex_exp(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_expm1(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_exp_pi_i(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_exp2(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_exp10(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_log(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_log1p(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_log_pi_i(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_log2(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_log10(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_sin(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_cos(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_sin_cos(nfloat_complex_ptr res1, nfloat_complex_ptr res2, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_tan(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_cot(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_sec(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_csc(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_sin_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_cos_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_sin_cos_pi(nfloat_complex_ptr res1, nfloat_complex_ptr res2, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_tan_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_cot_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_sec_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_csc_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_sinc(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_sinc_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_sinh(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_cosh(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_sinh_cosh(nfloat_complex_ptr res1, nfloat_complex_ptr res2, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_tanh(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_coth(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_sech(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_csch(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_asin(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_acos(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_atan(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_acot(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_asec(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_acsc(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_asinh(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_acosh(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_atanh(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_acoth(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_asech(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_acsch(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_asin_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_acos_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_atan_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_acot_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_asec_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_acsc_pi(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx)
int nfloat_complex_pow(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)

Complex elementary functions, as compositions of the real functions,
with the branch cuts of `acb_t`. Arguments in the real domain of a
function use the real function.

::: details
The formulas avoid cancellation except where the function itself is
ill-conditioned. For example,
$\tan(a + bi) = (\sin a \cos a + i \sinh b \cosh b)/(\cos^2 a + \sinh^2 b)$
for $|b| < 1/2$ and, with $E = e^{-2|b|}$,
$\tan(a + bi) = (2E\sin 2a + i\operatorname{sgn}(b)(1 - E^2))/(1 + 2E\cos 2a + E^2)$
for $|b| \ge 1/2$ (similarly for `cot`; `sec` and `csc` use the analog
of the first formula); $\log|z| = \log(a^2 + b^2)/2$, or
$\operatorname{log1p}(d)/2$ with $d = a^2 + b^2 - 1$ computed exactly in
fixed point when $\max(|a|, |b|) \in [1/2, 2)$ (also used for
$1 - |z|^2$ in `atan`; in `log1p`, $|1 + z|^2 - 1 = a(2 + a) + b^2$ is
computed in double precision when the terms cancel); `atan` by
Kahan-type formulas with `atan2` and `log1p`; and `asin`, `acos` and
`acosh` by Kahan's formulas using square roots of $1 \pm z$. The
hyperbolic functions and their inverses are evaluated as trigonometric
functions of $iz$, and `acot`, `asec`, `acsc` (and the hyperbolic
versions) as functions of $1/z$.

There are no internal guard bits: a small component of a large result
can have a large relative error, and the functions of $\pi z$ can have
an error of order $|\pi \operatorname{Im}(z)|$ ulp from rounding
$\pi \operatorname{Im}(z)$. Intermediate squares can overflow or
underflow (giving `GR_UNABLE`) for arguments near the limits of the
exponent range. Integer powers with $|y| < 2^{16}$ are computed by
binary exponentiation, other powers as $\exp(y \log x)$.
:::
:::

::: function
void _nfloat_complex_vec_init(nfloat_complex_ptr res, slong len, gr_ctx_t ctx)
void _nfloat_complex_vec_clear(nfloat_complex_ptr res, slong len, gr_ctx_t ctx)
int _nfloat_complex_vec_zero(nfloat_complex_ptr res, slong len, gr_ctx_t ctx)
int _nfloat_complex_vec_set(nfloat_complex_ptr res, nfloat_complex_srcptr x, slong len, gr_ctx_t ctx)
int _nfloat_complex_vec_add(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, slong len, gr_ctx_t ctx)
int _nfloat_complex_vec_sub(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, slong len, gr_ctx_t ctx)
int _nfloat_complex_vec_dot(nfloat_complex_ptr res, nfloat_complex_srcptr initial, int subtract, nfloat_complex_srcptr x, nfloat_complex_srcptr y, slong len, gr_ctx_t ctx)
int _nfloat_complex_vec_dot_rev(nfloat_complex_ptr res, nfloat_complex_srcptr initial, int subtract, nfloat_complex_srcptr x, nfloat_complex_srcptr y, slong len, gr_ctx_t ctx)

Vector operations on complex numbers, as for the real vectors; the dot
products accumulate the real and imaginary parts in fixed point.
:::

::: function
int nfloat_complex_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
int nfloat_complex_mat_mul_fixed(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, slong max_extra_prec, gr_ctx_t ctx)
int nfloat_complex_mat_mul_block(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, slong min_block_size, gr_ctx_t ctx)
int nfloat_complex_mat_mul_reorder(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)

Complex matrix multiplication, as for the real matrices. The fixed-point
and block versions use three real products (the Karatsuba formula
$(AC - BD) + ((A + B)(C + D) - AC - BD)i$);
`nfloat_complex_mat_mul_reorder` separates the real and imaginary parts
into real matrices and multiplies them with `nfloat_mat_mul` (with fewer
products when $A$ or $B$ is real). `nfloat_complex_mat_mul` uses the
reordering for real factors and for parts that are not well scaled.
:::

::: function
int nfloat_complex_mat_nonsingular_solve_tril(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
int nfloat_complex_mat_nonsingular_solve_triu(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
int nfloat_complex_mat_lu(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)

Triangular solving and LU factorisation, choosing between the classical
and recursive generic algorithms by precision-dependent cutoffs.
:::

## Fixed-point arithmetic {#fixed-point}

The matrix multiplication of `nfloat_mat_mul_fixed` works on packed
fixed-point numbers in $(-1, 1)$ with $n$-word precision, stored as
$n + 1$ words: a sign word followed by the $n$ words of the absolute
value, least significant first. In these functions, `nlimbs` is $n$
while the storage is `nlimbs + 1` words. There is no overflow handling:
inputs must be scaled so that all intermediate results, including
rounding errors, stay in $(-1, 1)$, which `_nfixed_mat_mul_bound`
checks.

::: function
void _nfixed_print(nn_srcptr x, slong nlimbs, slong exp)

Prints the fixed-point number `x` times $2^{\mathit{exp}}$ in decimal.
:::

::: function
void _nfixed_vec_add(nn_ptr res, nn_srcptr a, nn_srcptr b, slong len, slong nlimbs)
void _nfixed_vec_sub(nn_ptr res, nn_srcptr a, nn_srcptr b, slong len, slong nlimbs)

Elementwise sum and difference of vectors of `len` fixed-point numbers
(exact).
:::

::: function
void _nfixed_dot_2(nn_ptr res, nn_srcptr x, slong xstride, nn_srcptr y, slong ystride, slong len)
void _nfixed_dot_3(nn_ptr res, nn_srcptr x, slong xstride, nn_srcptr y, slong ystride, slong len)
void _nfixed_dot_4(nn_ptr res, nn_srcptr x, slong xstride, nn_srcptr y, slong ystride, slong len)
void _nfixed_dot_5(nn_ptr res, nn_srcptr x, slong xstride, nn_srcptr y, slong ystride, slong len)
void _nfixed_dot_6(nn_ptr res, nn_srcptr x, slong xstride, nn_srcptr y, slong ystride, slong len)
void _nfixed_dot_7(nn_ptr res, nn_srcptr x, slong xstride, nn_srcptr y, slong ystride, slong len)
void _nfixed_dot_8(nn_ptr res, nn_srcptr x, slong xstride, nn_srcptr y, slong ystride, slong len)

Dot products of fixed-point vectors with 2 to 8 words, accumulating
truncated high products in separate positive and negative sums. The
strides are the offsets in words between consecutive entries and may be
negative.
:::

::: function
void _nfixed_mat_mul(nn_ptr C, nn_srcptr A, nn_srcptr B, slong m, slong n, slong p, slong nlimbs)
void _nfixed_mat_mul_classical(nn_ptr C, nn_srcptr A, nn_srcptr B, slong m, slong n, slong p, slong nlimbs)
void _nfixed_mat_mul_classical_precise(nn_ptr C, nn_srcptr A, nn_srcptr B, slong m, slong n, slong p, slong nlimbs)
void _nfixed_mat_mul_waksman(nn_ptr C, nn_srcptr A, nn_srcptr B, slong m, slong n, slong p, slong nlimbs)
void _nfixed_mat_mul_strassen(nn_ptr C, nn_srcptr A, nn_srcptr B, slong m, slong n, slong p, slong cutoff, slong nlimbs)

Set the $m \times p$ matrix $C$ to the product of the $m \times n$
matrix $A$ and the $n \times p$ matrix $B$ (dense, row-major). The
classical version uses the dot products; the `_precise` version computes
with one extra word and is intended for testing. Waksman's algorithm
[@Waksman1970] halves the number of multiplications at the cost of more
additions, which pays off from about 9 words. Strassen's algorithm
[@Str1969] recurses down to `cutoff` (a negative value selects the
default). `_nfixed_mat_mul` chooses: Strassen for dimensions above
about 26 to 57 (depending on the precision and parity), otherwise Waksman
or classical.
:::

::: function
void _nfixed_mat_mul_bound(double * bound, double * error, slong m, slong n, slong p, double A, double B, slong nlimbs)
void _nfixed_mat_mul_bound_classical(double * bound, double * error, slong m, slong n, slong p, double A, double B, slong nlimbs)
void _nfixed_mat_mul_bound_waksman(double * bound, double * error, slong m, slong n, slong p, double A, double B, slong nlimbs)
void _nfixed_mat_mul_bound_strassen(double * bound, double * error, slong m, slong n, slong p, double A, double B, slong cutoff, slong nlimbs)
void _nfixed_complex_mat_mul_bound(double * bound, double * error, slong m, slong n, slong p, double A, double B, double C, double D, slong nlimbs)

For an $m \times n \times p$ product at precision `nlimbs` with entries
of $A$ and $B$ bounded by `A` and `B`, set `bound` to a bound for all
intermediate values of the corresponding algorithm, which must be less
than 1 for the computation to be valid, and `error` to a bound for the
output error in ulps. `_nfixed_mat_mul_bound` follows the algorithm
choice of `_nfixed_mat_mul`. The bounds are nondecreasing in `A` and
`B`. For complex matrices, the entries are bounded by $A + Bi$ and
$C + Di$.

::: details
A truncated product of two $n$-word numbers has an error of at most
$2n - 1$ ulps, so the classical algorithm has error at most
$n(2\,\mathit{nlimbs} - 1)$ ulps for an inner dimension $n$ and
intermediate values bounded by $nAB$. Waksman's and Strassen's
algorithms form sums and differences of entries before multiplying,
which increases both the bound on intermediate values (for example by
$6\lfloor n/2 \rfloor (A + B)^2$ for Waksman) and the error; the
fixed-point multiplication scales its inputs down until the bound is
below 1.
:::
:::

## Internal functions {#internal}

::: function
int _nfloat_underflow(nfloat_ptr res, int sgnbit, gr_ctx_t ctx)
int _nfloat_overflow(nfloat_ptr res, int sgnbit, gr_ctx_t ctx)

Set `res` to the result of an underflow (zero, returning `GR_SUCCESS`
only with `NFLOAT_ALLOW_UNDERFLOW`) or an overflow (a signed infinity
with `NFLOAT_ALLOW_INF`, otherwise zero and `GR_UNABLE`).
:::

::: macro
NFLOAT_HANDLE_UNDERFLOW(res, ctx)
NFLOAT_HANDLE_OVERFLOW(res, ctx)
NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx)

Return from the enclosing function with the result of
`_nfloat_underflow` or `_nfloat_overflow` if the exponent of `res` is
out of range.
:::

::: function
int _nfloat_cmp(nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)
int _nfloat_cmpabs(nfloat_srcptr x, nfloat_srcptr y, gr_ctx_t ctx)

Return $-1$, 0 or 1 according to the comparison of $x$ and $y$, or of
$|x|$ and $|y|$; neither may be NaN.
:::

::: function
int _nfloat_add_1(nfloat_ptr res, ulong x0, slong xexp, int xsgnbit, ulong y0, slong delta, gr_ctx_t ctx)
int _nfloat_sub_1(nfloat_ptr res, ulong x0, slong xexp, int xsgnbit, ulong y0, slong delta, gr_ctx_t ctx)
int _nfloat_add_2(nfloat_ptr res, nn_srcptr xd, slong xexp, int xsgnbit, nn_srcptr yd, slong delta, gr_ctx_t ctx)
int _nfloat_sub_2(nfloat_ptr res, nn_srcptr xd, slong xexp, int xsgnbit, nn_srcptr yd, slong delta, gr_ctx_t ctx)
int _nfloat_add_3(nfloat_ptr res, nn_srcptr x, slong xexp, int xsgnbit, nn_srcptr y, slong delta, gr_ctx_t ctx)
int _nfloat_sub_3(nfloat_ptr res, nn_srcptr x, slong xexp, int xsgnbit, nn_srcptr y, slong delta, gr_ctx_t ctx)
int _nfloat_add_4(nfloat_ptr res, nn_srcptr x, slong xexp, int xsgnbit, nn_srcptr y, slong delta, gr_ctx_t ctx)
int _nfloat_sub_4(nfloat_ptr res, nn_srcptr x, slong xexp, int xsgnbit, nn_srcptr y, slong delta, gr_ctx_t ctx)
int _nfloat_add_n(nfloat_ptr res, nn_srcptr xd, slong xexp, int xsgnbit, nn_srcptr yd, slong delta, slong nlimbs, gr_ctx_t ctx)
int _nfloat_sub_n(nfloat_ptr res, nn_srcptr xd, slong xexp, int xsgnbit, nn_srcptr yd, slong delta, slong nlimbs, gr_ctx_t ctx)

The fixed-size kernels of addition and subtraction of magnitudes: set
`res` to $(-1)^{\mathit{xsgnbit}} (|x| \pm |y|)$, where $x$ has
mantissa `x` and exponent `xexp`, and $y$ has mantissa `y` and exponent
`xexp - delta` with $\mathit{delta} \ge 0$ (so $|x| \ge |y|$ in
exponent, and for subtraction $|x| \ge |y|$). The mantissas have 1, 2,
3, 4 or `nlimbs` words.
:::

::: function
int _nfloat_complex_sqr_naive(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr a, nfloat_srcptr b, gr_ctx_t ctx)
int _nfloat_complex_sqr_standard(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr a, nfloat_srcptr b, gr_ctx_t ctx)
int _nfloat_complex_sqr_karatsuba(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr a, nfloat_srcptr b, gr_ctx_t ctx)
int _nfloat_complex_sqr(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr a, nfloat_srcptr b, gr_ctx_t ctx)
int _nfloat_complex_mul_naive(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr a, nfloat_srcptr b, nfloat_srcptr c, nfloat_srcptr d, gr_ctx_t ctx)
int _nfloat_complex_mul_standard(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr a, nfloat_srcptr b, nfloat_srcptr c, nfloat_srcptr d, gr_ctx_t ctx)
int _nfloat_complex_mul_karatsuba(nfloat_ptr res1, nfloat_ptr res2, nfloat_srcptr a, nfloat_srcptr b, nfloat_srcptr c, nfloat_srcptr d, gr_ctx_t ctx)

Set `res1` and `res2` to the real and imaginary parts of $(a + bi)^2$
or $(a + bi)(c + di)$ with real `nfloat` parts: `naive` with real
products and sums, `standard` with fused dot products of two terms,
and `karatsuba` with three real products instead of four.
`_nfloat_complex_sqr` uses the naive version when infinities or NaN are
enabled, real squarings when a part is zero, and otherwise the standard
version below 20 words and the Karatsuba version from 20 words.
:::
