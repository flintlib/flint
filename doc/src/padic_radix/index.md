---
title: padic_radix.h
header: padic_radix.h
description: $p$-adic numbers for word-size $p$ in a radix $p^e$ representation, with tracked error terms, absolute and relative precision, and the generic ring interface.
category: p-adic numbers
order: 5
---

This module implements $\mathbf{Q}_p$ for a word-size prime $p$. An
element is a ball $p^v u + O(p^N)$, or an exact element without error
term, where the unit $u$ is a `radix_integer_t` whose limbs are digits in
the radix $B = p^e$ ([](radix)); the limb radix is chosen as large as
possible, for example $7^{22}$ for $p = 7$ on a 64-bit machine. The
valuation $v$ and the precision $N$ are counted in powers of $p$, as in
`padic`. Unlike `padic`, every function propagates the error terms of its
inputs, so a result never claims more digits than are determined, and an
exact computation stays exact as long as the result fits the precision
of the context. Arithmetic on units is a truncation of limb arrays rather
than a division, which makes the module several times faster than
`padic` at high precision.

The module implements the generic ring interface (`gr.h`): the context
is a `gr_ctx_t`, functions return the status flags `GR_SUCCESS`,
`GR_DOMAIN` and `GR_UNABLE`, predicates return a `truth_t`, and
polynomials, matrices and power series over $\mathbf{Q}_p$ are available
as `gr_poly`, `gr_mat` and `gr_series`. The guide
[](padic_radix/getting-started) shows the precision model at work.

::: summary algorithms/p-adic/precision :::

## Precision {#precision}

A context has a default *relative precision* $r$ and *absolute
precision* $a$, either of which may be $+\infty$
(`PADIC_RADIX_PREC_INF`). Every result $p^v u + O(p^N)$ is truncated to
the precision $\min(N, a, v + r)$, and an exact result is truncated only
if it does not fit within $\min(a, v + r)$; a truncation is recorded in
the error term. With $a = +\infty$ the arithmetic is floating-point-like
with $r$ digits; with $r = +\infty$ it is fixed-point modulo $p^a$, as in
`padic`; with both infinite all arithmetic is exact in $\Z[1/p]$, and
operations whose result would be an infinite expansion (such as $1/3$)
return `GR_UNABLE`. The error terms are propagated by the rules of
[](algorithms/p-adic/precision#propagation): the precision of a sum is
the smaller absolute precision, that of a product, quotient or square
root the smaller relative precision, adjusted by the valuations.

A nonzero unit is never divisible by $p$; a zero unit has $v = 0$. For an
inexact element the unit is the nonnegative residue modulo
$p^{N - v}$. Exact elements may have negative units only if the context
has the flag `PADIC_RADIX_SIGNED`; otherwise a negative exact value is
stored as an inexact residue if the context precision is finite.

::: macro
PADIC_RADIX_EXACT

The value of $N$ that marks an exact element; it equals `WORD_MAX`.
:::

::: macro
PADIC_RADIX_PREC_INF

The precision value $+\infty$ for the context parameters, equal to
`WORD_MAX`.
:::

::: macro
PADIC_RADIX_ERR_MAX

The largest finite precision stored, `WORD_MAX / 4`. Larger precisions
are clamped to it, which is always sound; valuations beyond
$\pm$`PADIC_RADIX_ERR_MAX` and precisions below its negative cannot be
represented, and the operations that would produce them return
`GR_UNABLE`. The bound keeps sums of valuations and precisions from
overflowing a `slong`.
:::

## Context {#context}

::: function
int gr_ctx_init_padic_radix(gr_ctx_t ctx, ulong p, slong prec_rel, slong prec_abs, int flags)

Initialises `ctx` for $\mathbf{Q}_p$ with default relative precision
`prec_rel` and absolute precision `prec_abs`, either of which may be
`PADIC_RADIX_PREC_INF`, and the given combination of the flags below.
Negative relative precisions are raised to 0, and finite precisions are
clamped to $\pm$`PADIC_RADIX_ERR_MAX`. That $p$ is prime is only checked
in builds with assertions. Returns `GR_SUCCESS`. The context is a field
in the sense of `gr_ctx_is_field`, and is exact if and only if both
precisions are infinite.
:::

::: macro
PADIC_RADIX_SIGNED
PADIC_RADIX_NO_ERROR
PADIC_RADIX_DECIMAL
PADIC_RADIX_TEST_LIMITS

Context flags. `PADIC_RADIX_SIGNED` allows exact elements with negative
units. `PADIC_RADIX_NO_ERROR` is reserved for a future floating-point
mode without error terms and currently has no effect.
`PADIC_RADIX_DECIMAL` prints units as decimal integers rather than in
the default form of `radix_get_str_sum`. `PADIC_RADIX_TEST_LIMITS` makes
`padic_radix_randtest` produce valuations and precisions near
`PADIC_RADIX_ERR_MAX`.
:::

::: function
int gr_ctx_init_padic_radix_randtest(gr_ctx_t ctx, flint_rand_t state, slong maxprec)

Initialises `ctx` with a random prime and a random configuration for
testing: absolute precision only, relative precision only, both (with
precisions in $[1, \mathtt{maxprec}]$) or neither, and random flags
`PADIC_RADIX_SIGNED` and `PADIC_RADIX_TEST_LIMITS`.
:::

::: type
padic_radix_ctx_struct

The data of a context: `radix` (the `radix_struct` for digit radix $p$
and limb radix $p^e$), `p`, `prec_abs`, `prec_rel` and `flags`.
:::

::: macro
GR_PADIC_RADIX_CTX(ctx)

Pointer to the `padic_radix_ctx_struct` of a `gr_ctx_t`.
:::

::: macro
PADIC_RADIX_CTX_RADIX(ctx)

Pointer to the `radix_struct` used for the limb arithmetic of the
context.
:::

::: macro
PADIC_RADIX_CTX_PREC_ABS(ctx)
PADIC_RADIX_CTX_PREC_REL(ctx)
PADIC_RADIX_CTX_FLAGS(ctx)
PADIC_RADIX_CTX_SIGNED(ctx)
PADIC_RADIX_CTX_DECIMAL(ctx)

The precisions and flags of the context; the last two test single
flags.
:::

::: function
void padic_radix_ctx_clear(gr_ctx_t ctx)
int padic_radix_ctx_write(gr_stream_t out, gr_ctx_t ctx)

Implementations of `gr_ctx_clear` and `gr_ctx_write`; the description
written has the form `Radix 7-adic numbers (rel prec 10, abs prec inf)`.
:::

## Elements {#elements}

::: type
padic_radix_struct
padic_radix_t

An element $p^v u + O(p^N)$: the fields are `u` (a
`radix_integer_struct`), `v` and `N`, with $N =$ `PADIC_RADIX_EXACT` for
exact elements. A `padic_radix_t` is an array of one
`padic_radix_struct`.
:::

::: macro
PADIC_RADIX_UNIT(x)
PADIC_RADIX_VAL(x)
PADIC_RADIX_N(x)

The unit (as a pointer), the valuation and the precision of `x`.
:::

The following functions implement the generic ring methods; they are
also accessible through the corresponding `gr` functions. All of them
allow aliasing.

::: function
void padic_radix_init(padic_radix_t res, gr_ctx_t ctx)
void padic_radix_clear(padic_radix_t res, gr_ctx_t ctx)
void padic_radix_swap(padic_radix_t x, padic_radix_t y, gr_ctx_t ctx)
void padic_radix_set_shallow(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)

Initialise to the exact zero, clear, swap in constant time, and copy
the structure without copying the limbs.
:::

::: function
slong padic_radix_get_error(const padic_radix_t x, gr_ctx_t ctx)
truth_t padic_radix_is_exact(const padic_radix_t x, gr_ctx_t ctx)

Return the precision $N$ of `x` (`PADIC_RADIX_EXACT` for an exact
element), respectively whether `x` is exact.
:::

::: function
int padic_radix_randtest(padic_radix_t res, flint_rand_t state, gr_ctx_t ctx)

Sets `res` to a random element with small valuation, exact or with a
random error term (never with an error term in an exact context).
:::

::: function
int padic_radix_write(gr_stream_t out, const padic_radix_t x, gr_ctx_t ctx)

Writes `x` in the form `(u) * p^v + O(p^N)`, omitting the power of $p$
when $v = 0$ and the error term for exact elements, for example
`188316833 + O(7^10)`. The unit is written by `radix_get_str_sum`, or
as a decimal integer if the context has the flag `PADIC_RADIX_DECIMAL`.
Strings of this form, including `O(p^n)`, are accepted by `gr_set_str`.
:::

::: function
int padic_radix_zero(padic_radix_t res, gr_ctx_t ctx)
int padic_radix_one(padic_radix_t res, gr_ctx_t ctx)

Set `res` to the exact zero, respectively to one, which becomes
$0 + O(p^{t})$ if $t = \min(a, r)$ is not positive for the context
precisions $a$ and $r$.
:::

::: function
int padic_radix_set(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)

Sets `res` to `x`, truncated to the context precision.
:::

::: function
int padic_radix_set_ui(padic_radix_t res, ulong x, gr_ctx_t ctx)
int padic_radix_set_si(padic_radix_t res, slong x, gr_ctx_t ctx)
int padic_radix_set_fmpz(padic_radix_t res, const fmpz_t x, gr_ctx_t ctx)

Set `res` to the integer $x$, truncated to the context precision. A
negative $x$ in a context without `PADIC_RADIX_SIGNED` is stored as an
inexact residue, or gives `GR_UNABLE` if the context precision is
infinite.
:::

::: function
int padic_radix_exact_set_ui(padic_radix_t res, ulong x, gr_ctx_t ctx)
int padic_radix_exact_set_si(padic_radix_t res, slong x, gr_ctx_t ctx)
int padic_radix_exact_set_fmpz(padic_radix_t res, const fmpz_t x, gr_ctx_t ctx)

Set `res` to the exact integer $x$, ignoring the context precision; the
unit may be negative regardless of the flags.
:::

::: function
int padic_radix_get_fmpz(fmpz_t res, const padic_radix_t x, gr_ctx_t ctx)

Sets `res` to the value of `x` if `x` is an exact integer. Returns
`GR_DOMAIN` if $v < 0$ and `GR_UNABLE` if `x` is inexact or its
valuation exceeds the size limit of the context.
:::

::: function
int padic_radix_neg(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
int padic_radix_add(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx)
int padic_radix_sub(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx)
int padic_radix_mul(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx)

Set `res` to $-x$, $x + y$, $x - y$ and $xy$. The error term of a sum is
$O(p^{\min(N_x, N_y)})$ and that of a product
$O(p^{\min(N_x + v_y, N_y + v_x)})$; the result is then truncated to the
context precision. An exact zero factor gives an exact zero.
:::

::: function
int padic_radix_inv(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
int padic_radix_div(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx)

Set `res` to $1/x$, respectively $x/y$. The relative precision of the
result is the smaller relative precision of the operands, capped by the
context. Exact operands give an exact result when the quotient
terminates (for example $1/49$ or $6/3$). Return `GR_DOMAIN` when
dividing by an exact zero, and `GR_UNABLE` when the divisor is not known
to be nonzero or when the context is exact and the quotient has an
infinite expansion.
:::

::: function
truth_t padic_radix_is_invertible(const padic_radix_t x, gr_ctx_t ctx)

Returns whether `x` is nonzero: `T_TRUE` if it is exact and nonzero or
its leading digit is known, `T_FALSE` for an exact zero and `T_UNKNOWN`
otherwise.
:::

::: function
int padic_radix_sqrt(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
int padic_radix_rsqrt(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)

Set `res` to a square root $\sqrt x$, respectively $1/\sqrt x$, with
relative precision equal to that of $x$ (one digit less for $p = 2$),
capped by the context. The root is a fixed canonical branch, so that a
result at lower precision is the truncation of the result at higher
precision, and `rsqrt` gives the reciprocal of `sqrt`. An exact perfect
square gives its exact root. Return `GR_DOMAIN` if $x$ has odd
valuation or its unit is not a square (and for `rsqrt` if $x$ is an
exact zero), and `GR_UNABLE` if $x$ is not known to be nonzero, if
$p = 2$ and the unit is known to fewer than 3 digits, or if the context
is exact and the root is irrational. `sqrt` of an exact zero is zero.
:::

::: function
truth_t padic_radix_is_square(const padic_radix_t x, gr_ctx_t ctx)

Returns whether $x$ is a square: decided by the valuation and the unit
modulo $p$ (modulo 8 for $p = 2$), and `T_UNKNOWN` for an inexact zero
or, for $p = 2$, a unit known to fewer than 3 digits.
:::

::: function
truth_t padic_radix_is_zero(const padic_radix_t x, gr_ctx_t ctx)
truth_t padic_radix_is_one(const padic_radix_t x, gr_ctx_t ctx)
truth_t padic_radix_is_neg_one(const padic_radix_t x, gr_ctx_t ctx)
truth_t padic_radix_equal(const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx)

Ball predicates: `T_TRUE` only for exact elements with the property,
`T_FALSE` if no value in the ball has it (for equality, if the balls are
disjoint), and `T_UNKNOWN` otherwise.
:::

## Dot products {#dot}

::: function
int padic_radix_dot_strided(padic_radix_t res, const padic_radix_t initial, int subtract, const padic_radix_struct * vec1, slong stride1, const padic_radix_struct * vec2, slong stride2, slong len, gr_ctx_t ctx)
int padic_radix_dot_strided_naive(padic_radix_t res, const padic_radix_t initial, int subtract, const padic_radix_struct * vec1, slong stride1, const padic_radix_struct * vec2, slong stride2, slong len, gr_ctx_t ctx)
int padic_radix_dot_strided_delayed(padic_radix_t res, const padic_radix_t initial, int subtract, const padic_radix_struct * vec1, slong stride1, const padic_radix_struct * vec2, slong stride2, slong len, gr_ctx_t ctx)

Set `res` to $s \pm \sum_{i < len} a_i b_i$ with
$a_i$ = `vec1[i * stride1]`, $b_i$ = `vec2[i * stride2]`, $s$ =
`initial` (zero if `initial` is `NULL`) and the sign $-$ if `subtract`
is nonzero. The *naive* version calls `padic_radix_mul` and
`padic_radix_add` in a loop. The *delayed* version finds the minimum
valuation and the precision of the result first, then accumulates the
limb products of all terms in three-word accumulators, grouped by the
residue of the digit offset modulo $e$ and by sign, and propagates
carries, shifts and divisions by $p^e$ once at the end. The default
uses the naive version for $len \le 7$.
:::

::: function
int padic_radix_dot(padic_radix_t res, const padic_radix_t initial, int subtract, const padic_radix_struct * vec1, const padic_radix_struct * vec2, slong len, gr_ctx_t ctx)
int padic_radix_dot_rev(padic_radix_t res, const padic_radix_t initial, int subtract, const padic_radix_struct * vec1, const padic_radix_struct * vec2, slong len, gr_ctx_t ctx)

The contiguous dot product and the dot product with `vec2` reversed,
as wrappers of `padic_radix_dot_strided`. These implement the generic
dot product methods, which `gr_mat_mul` and the `gr_poly` products use.
:::

## Exponential and logarithm {#exp-log}

The algorithms are those of `padic.h` in radix arithmetic; see
[](algorithms/p-adic/exp-log).

::: function
slong _padic_radix_exp_bound(slong v, slong N, ulong p)

Returns the number $n$ of terms after which the exponential series for
an argument of valuation $v$ can be truncated modulo $p^N$,
$n = \max(1, \lceil ((p-1)N - 1)/((p-1)v - 1) \rceil)$. Assumes that the
series converges ($v \ge 1$, or $v \ge 2$ for $p = 2$).
:::

::: function
void _padic_radix_exp_rectangular(radix_integer_t rop, const radix_integer_t u, slong v, slong N, const radix_t radix)
void _padic_radix_exp_balanced(radix_integer_t rop, const radix_integer_t u, slong v, slong N, const radix_t radix)
void _padic_radix_exp(radix_integer_t rop, const radix_integer_t u, slong v, slong N, const radix_t radix)

Set `rop` to the nonnegative residue of $\exp(u p^v)$ modulo $p^N$, for
$v \ge 1$ ($v \ge 2$ if $p = 2$), not checked. The *rectangular*
version uses rectangular splitting, the *balanced* version the $p$-adic
bit-burst method with binary splitting. `_padic_radix_exp` uses the
rectangular version for $N$ below a cutoff between 80 and 200 that
decreases with the size of $p$.
:::

::: function
int padic_radix_exp_rectangular(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
int padic_radix_exp_balanced(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
int padic_radix_exp(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)

Set `res` to $\exp(x)$, with precision $\min(N_x, a, r)$ for the context
precisions $a$ and $r$ (the result is a unit). Return `GR_DOMAIN` if the
series diverges, and `GR_UNABLE` if convergence cannot be decided (an
inexact zero of too low precision) or the context precision is infinite.
:::

::: function
slong _padic_radix_log_bound(slong v, slong N, ulong p)

Returns the least $b \ge 2$ such that $iv - \operatorname{ord}_p(i) \ge
N$ for all $i \ge b$, so that the series for $\log(1 - y)$ with
$\operatorname{ord}_p(y) = v$ can be truncated after $b - 1$ terms.
:::

::: function
void _padic_radix_log_rectangular(radix_integer_t rop, const radix_integer_t y, slong N, const radix_t radix)
void _padic_radix_log_balanced(radix_integer_t rop, const radix_integer_t y, slong N, const radix_t radix)
void _padic_radix_log(radix_integer_t rop, const radix_integer_t y, slong N, const radix_t radix)

Set `rop` to the residue of $\log(1 - y)$ modulo $p^N$ for an integer
$y$ of valuation at least 1 (not checked). `_padic_radix_log` uses the
rectangular version for $N$ below a cutoff between 80 and 350 that
decreases with the size of $p$.
:::

::: function
int padic_radix_log_rectangular(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
int padic_radix_log_balanced(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
int padic_radix_log(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)

Set `res` to $\log(x)$, with absolute precision $\min(N_x, a, r)$.
Return `GR_DOMAIN` unless $x \equiv 1 \pmod p$ (for every $p$, including
$p = 2$), and `GR_UNABLE` if this cannot be decided or if $x$ is exact
and the context precision infinite.
:::

## Internal functions {#internal}

::: function
int _padic_radix_finalize(padic_radix_t res, gr_ctx_t ctx)

Canonicalises a freshly computed element (moving factors of $p$ from
the unit to the valuation), enforces the representability bounds of
`PADIC_RADIX_ERR_MAX` and truncates to the context precision. Returns
`GR_UNABLE` if the valuation or precision cannot be represented.
:::

::: function
int _padic_radix_add_sub_reference(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, int sub, gr_ctx_t ctx)
int _padic_radix_mul_reference(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx)

Simple reference implementations of addition ($x + y$, or $x - y$ if
`sub` is nonzero) and multiplication, used to test the optimised
functions.
:::
