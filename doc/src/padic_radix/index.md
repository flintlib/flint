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

Valuations and precisions may differ by up to about `WORD_MAX / 2`. In a
context with a finite precision, the intermediate results of every
operation are bounded in terms of the working precision (that of the
context, or less when the inputs are less precise), whatever the
exponents of the inputs: the fast algorithms for dot products, matrix and
polynomial products and characteristic polynomials give up in favour of
classical ones where their temporaries would be much larger than the
precision of the outputs requires (see the individual functions). (With
only an absolute precision $a$, an element of valuation $v$ has a working
precision of $a - v$ digits, which is large for $v$ far below $a$.) In a
context without finite precision, the sizes are those of the exact
results, which can be arbitrarily large (for example
$1 + p^{\mathrm{WORD\_MAX}/8}$); they are bounded by the `size_limit` of
the context, in digits (unlimited by default): exact results that would
exceed it give `GR_UNABLE`, and inexact results keep at most `size_limit`
digits of relative precision.

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
PADIC_RADIX_FAST_BOUNDS
PADIC_RADIX_FAST_LINALG

Context flags. `PADIC_RADIX_SIGNED` allows exact elements with negative
units. `PADIC_RADIX_NO_ERROR` is reserved for a future floating-point
mode without error terms and currently has no effect.
`PADIC_RADIX_DECIMAL` prints units as decimal integers rather than in
the default form of `radix_get_str_sum`. `PADIC_RADIX_TEST_LIMITS` makes
`padic_radix_randtest` produce valuations and precisions near
`PADIC_RADIX_ERR_MAX`. `PADIC_RADIX_FAST_BOUNDS` lets polynomial
products use quasi-linear lower bounds for the precisions of the
coefficients when the exact computation would cost more (irregular
precision patterns in long products, see [](#block)); the results may
then have less precision than with the default exact bounds, never more.
`PADIC_RADIX_FAST_LINALG` makes the linear algebra choose the fastest
algorithms rather than those with the optimal precision, which is useful
at high working precision, where losing a few digits costs less than the
optimal algorithms: the characteristic polynomial uses the Hessenberg
reduction ($O(n^3)$ operations, a few digits lost on structured matrices)
instead of the exact lifting method or Berkowitz's algorithm (6 to 15
times faster for structured matrices at 1000 to 10000 digits and
dimensions 30 to 80, losing up to 4 digits), and linear systems use the
LU decomposition. Determinants are unaffected: their
optimal algorithm is as fast as an LU decomposition. Results are valid
either way.
:::

::: function
int gr_ctx_init_padic_radix_randtest(gr_ctx_t ctx, flint_rand_t state, slong maxprec)

Initialises `ctx` with a random prime and a random configuration for
testing: absolute precision only, relative precision only, both (with
precisions in $[1, \mathtt{maxprec}]$) or neither, and random flags
`PADIC_RADIX_SIGNED`, `PADIC_RADIX_TEST_LIMITS` and
`PADIC_RADIX_FAST_LINALG`.
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
context precision. An exact zero factor gives an exact zero. An inexact
zero $O(p^{N})$ enters the error term of a product with valuation $N$,
so that for example $O(p^{N_x}) \cdot O(p^{N_y}) = O(p^{N_x + N_y})$.
:::

::: function
int padic_radix_mul_si(padic_radix_t res, const padic_radix_t x, slong c, gr_ctx_t ctx)
int padic_radix_mul_ui(padic_radix_t res, const padic_radix_t x, ulong c, gr_ctx_t ctx)
int padic_radix_mul_fmpz(padic_radix_t res, const padic_radix_t x, const fmpz_t c, gr_ctx_t ctx)

Set `res` to $xc$, with the result of converting $c$ with
`padic_radix_set_si` (etc.) followed by `padic_radix_mul`.
:::

::: function
int padic_radix_addmul(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx)
int padic_radix_submul(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx)
int padic_radix_addmul_si(padic_radix_t res, const padic_radix_t x, slong c, gr_ctx_t ctx)
int padic_radix_submul_si(padic_radix_t res, const padic_radix_t x, slong c, gr_ctx_t ctx)
int padic_radix_addmul_ui(padic_radix_t res, const padic_radix_t x, ulong c, gr_ctx_t ctx)
int padic_radix_submul_ui(padic_radix_t res, const padic_radix_t x, ulong c, gr_ctx_t ctx)
int padic_radix_addmul_fmpz(padic_radix_t res, const padic_radix_t x, const fmpz_t c, gr_ctx_t ctx)
int padic_radix_submul_fmpz(padic_radix_t res, const padic_radix_t x, const fmpz_t c, gr_ctx_t ctx)
int _padic_radix_addmul_sub(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, int sub, gr_ctx_t ctx)

Set `res` to $\mathrm{res} \pm xy$ (respectively $\pm xc$), with exactly the
result of a multiplication into a temporary followed by an addition or
subtraction. Product digits that cannot reach the result (above its
precision) are not computed, and sums of at most two limbs are formed
inline without temporaries.
:::

::: function
int padic_radix_pow_ui(padic_radix_t res, const padic_radix_t x, ulong n, gr_ctx_t ctx)
int padic_radix_pow_si(padic_radix_t res, const padic_radix_t x, slong n, gr_ctx_t ctx)
int padic_radix_pow_fmpz(padic_radix_t res, const padic_radix_t x, const fmpz_t n, gr_ctx_t ctx)

Set `res` to $x^n$, computing the power of the unit with truncated
multiplications and rounding once. For $x = u p^v + O(p^N)$ and $n \ge 1$
the error term is $O(p^{(n - 1)v + N + v_p(n)})$, which is valid since
every term $\binom{n}{k} a^{n - k} h^k$, $k \ge 1$, of $(a + h)^n - a^n$
has valuation at least $(n - 1)v + N + v_p(n)$; for an inexact zero
$O(p^N)$ it is $O(p^{nN})$. Negative exponents invert first. An exact
input keeps an exact result when it fits the context precision.
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
int _padic_radix_vec_mul_scalar(padic_radix_struct * res, const padic_radix_struct * vec, slong len, const padic_radix_t c, gr_ctx_t ctx)
int _padic_radix_vec_mul_scalar_si(padic_radix_struct * res, const padic_radix_struct * vec, slong len, slong c, gr_ctx_t ctx)
int _padic_radix_vec_mul_scalar_ui(padic_radix_struct * res, const padic_radix_struct * vec, slong len, ulong c, gr_ctx_t ctx)
int _padic_radix_vec_mul_scalar_fmpz(padic_radix_struct * res, const padic_radix_struct * vec, slong len, const fmpz_t c, gr_ctx_t ctx)
int _padic_radix_vec_addmul_scalar(padic_radix_struct * res, const padic_radix_struct * vec, slong len, const padic_radix_t c, gr_ctx_t ctx)
int _padic_radix_vec_submul_scalar(padic_radix_struct * res, const padic_radix_struct * vec, slong len, const padic_radix_t c, gr_ctx_t ctx)
int _padic_radix_vec_addmul_scalar_si(padic_radix_struct * res, const padic_radix_struct * vec, slong len, slong c, gr_ctx_t ctx)
int _padic_radix_vec_submul_scalar_si(padic_radix_struct * res, const padic_radix_struct * vec, slong len, slong c, gr_ctx_t ctx)
int _padic_radix_vec_addmul_scalar_fmpz(padic_radix_struct * res, const padic_radix_struct * vec, slong len, const fmpz_t c, gr_ctx_t ctx)

Elementwise $\mathrm{res}_i = \mathrm{vec}_i c$ and
$\mathrm{res}_i = \mathrm{res}_i \pm \mathrm{vec}_i c$, with the results of
the corresponding scalar functions; the scalar is converted once.
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
residue of the digit offset modulo $e$ and by sign; at the end the
residue classes are combined in one pass (each accumulator multiplied by
$p^r$ for its residue $r$) and the carries are propagated once. The
default uses the naive version for at most 4 terms at a precision of one
limb, 3 terms at two limbs and 2 terms above (relative precision, or
absolute if there is none), and, at precisions of up to two limbs, also
for up to 7 terms whose valuations $v(a_i) + v(b_i)$ are not all equal.

When all terms are exact, the delayed version keeps the result exact if
it fits in the context precision (as `padic_radix_mul` and
`padic_radix_add` do), provided that the a priori extent of the terms,
$\max_i (v(a_i) + e |a_i| + v(b_i) + e |b_i|)$ for units of $|a_i|$ and
$|b_i|$ limbs, is at most one limb beyond the context precision;
otherwise the sum is rounded to the context precision. In a context
without finite precision, when the window of the delayed version (the
extent of the exact terms, or the precision they determine) exceeds the
`size_limit` of the context, it calls the naive version, which builds
only the partial sums (and still gets exact cancellations of terms with
huge valuations right).
:::

::: function
int padic_radix_dot(padic_radix_t res, const padic_radix_t initial, int subtract, const padic_radix_struct * vec1, const padic_radix_struct * vec2, slong len, gr_ctx_t ctx)
int padic_radix_dot_rev(padic_radix_t res, const padic_radix_t initial, int subtract, const padic_radix_struct * vec1, const padic_radix_struct * vec2, slong len, gr_ctx_t ctx)

The contiguous dot product and the dot product with `vec2` reversed,
as wrappers of `padic_radix_dot_strided`. These implement the generic
dot product methods, which `gr_mat_mul` and the `gr_poly` products use.
:::

## Matrix and polynomial products {#block}

The block products use a packed *limb-aligned* form internally: an element
$u p^v$ is stored relative to the limb frame $w = \lfloor v/e \rfloor$ as the
residue $u p^{v - ew} \bmod B^L$, so that all terms of an output land at
whole-limb offsets of a single accumulator, without digit shifts. Each
input is converted once, the precision of every output is computed up
front from arrays of exponents (min-plus products or convolutions), and
at precisions of up to two limbs the accumulation runs in registers.
Every output is computed with exactly the semantics of
`padic_radix_dot_strided_delayed` (the same precision, and the same
value modulo it, including exact sums of exact terms). The packed versions return `GR_UNABLE` without
modifying the output when they do not apply: in contexts without a
finite precision (where exact results need full products) and for
elements with exponents beyond $\pm\mathrm{WORD\_MAX}/16$.

::: function
int padic_radix_mat_mul_packed(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
int padic_radix_mat_mul_fmpz(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
int padic_radix_mat_mul_mpn_mod(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
int padic_radix_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)

Set $C = AB$. The packed version accumulates the limb products of every
output directly. The *fmpz* version scales row $i$ of $A$ and column $j$
of $B$ to base frames $r_i = \min_t w(A_{it})$ and $c_j = \min_t w(B_{tj})$,
truncates the residues to the largest output window, multiplies the
resulting integer matrices with `nmod_mat_mul` (modulo $B$, when every
window is a single limb) or `fmpz_mat_mul` (multimodular and blocked
algorithms), and reduces every output modulo $B^s$ for its window of $s$
limbs before the conversion to radix $B$. The *mpn_mod* version is the
same with the products done modulo $B^S$ for the largest window of $S$
limbs by `mpn_mod_mat_mul` (classical with delayed reduction, Waksman or
multimodular), and returns `GR_UNABLE` if $B^S$ does not have between 2 and
`MPN_MOD_MAX_LIMBS` binary limbs; `padic_radix_mat_mul` uses it in that
range, the nmod version for one-limb windows and the fmpz version above.
Outputs whose terms are all exact are summed exactly, in windows covering
the a priori extents of their terms plus a carry and a sign limb, when
`padic_radix_dot_strided_delayed` would keep them exact. Its precisions are computed
with blockwise minima of the exponents, in $O(mn \cdot \text{blocks})$
operations for outputs whose row of $A$ or column of $B$ has exponents
that are constant on blocks of the inner index (e.g. uniform valuations
and precisions), and with the exact $O(k)$ minimum for the others; the
results are identical to those of the packed version. The base frames
$r_i + c_j$ can be far below the frame of an output whose terms all have
higher valuations; the fmpz and mpn_mod versions return `GR_UNABLE` when
the largest window in the base frames exceeds four times (plus 8 limbs)
the largest window of an output in its own frame, since the integers
would then be much larger than the precision requires (arbitrarily large
for valuations spread beyond the precision).

`padic_radix_mat_mul` implements the generic method `gr_mat_mul`: the
classical algorithm for inner dimension at most 3, the modular versions
above tuned crossovers (estimated from an upper bound for the output
windows: from dimension 8 when all windows are one limb, about 32 at two
limbs, 16 at 3 to 8 limbs, 28 up to 15 limbs and 48 beyond, raised when the
limb frames of the entries are spread out and when the exact precision
bookkeeping is needed), and the packed version otherwise.
:::

::: function
int _padic_radix_poly_mulmid_KS(padic_radix_struct * res, const padic_radix_struct * poly1, slong len1, const padic_radix_struct * poly2, slong len2, slong nlo, slong nhi, gr_ctx_t ctx)
int _padic_radix_poly_mulmid_mod(padic_radix_struct * res, const padic_radix_struct * poly1, slong len1, const padic_radix_struct * poly2, slong len2, slong nlo, slong nhi, gr_ctx_t ctx)

Middle product by Kronecker substitution in radix $B$: the coefficients,
scaled to the base frames of the two polynomials and truncated to the
largest output window of $S$ limbs, are placed in slots of $2S + O(1)$
limbs of a single radix integer, and one `radix_mulmid` (which uses the
small-prime FFT directly in radix $B$ for large operands, without
conversions to binary) gives every output in its own slot. The *mod*
version instead multiplies the residues modulo $B^S$, avoiding the
padding of the slots: with `_nmod_poly_mulmid` when $S = 1$ and with
`_mpn_mod_poly_mulmid` (classical, Karatsuba or the small-prime FFT, as
tuned for `mpn_mod`) when $B^S$ has at most `MPN_MOD_MAX_LIMBS` binary
limbs, converting to binary and back; otherwise it returns `GR_UNABLE`.
The output precisions come from min-plus convolutions of the exponent
sequences (below), and coefficients whose terms are all exact are summed
exactly (in windows covering the extents of the terms) when
`padic_radix_dot_strided_delayed` would keep them exact. The results are
those of `_padic_radix_poly_mulmid_packed`, except with the context flag
`PADIC_RADIX_FAST_BOUNDS`, which allows smaller (but never larger)
precisions. Both return `GR_UNABLE` in the same cases as the packed
version. Outputs whose window in the base frames exceeds four times (plus
8 limbs) the largest window of an output in its own frame (all of whose
terms have valuations far above the base frames, e.g. next to a few
coefficients of huge valuation) are computed separately by
`padic_radix_dot_strided_delayed` when there are few of them (costing at
most about the transform and a sixteenth of the classical algorithm);
otherwise, that is when the valuations vary beyond the precision along
the outputs (alternating, trending or random valuations spread over more
digits than the precision), both return `GR_UNABLE`. High-valuation
coefficients that are beyond the precision of every output they
contribute to cost nothing.
:::

The per-coefficient precisions of a polynomial product are min-plus
convolutions $c_k = \min_i (a_i + b_{k-i})$ of exponent sequences, with
$+\infty$ for absent terms (exact zeros, or exact coefficients in the
error terms). They are always computed exactly by default:

* directly, over the pairs of finite entries only, when there are few of
  these (sparse sequences, e.g. with exact zero padding, cost only the
  number of pairs);
* by the lower convex hulls (the Minkowski sum of their slopes, in linear
  time) when both sequences are convex without gaps (in particular
  constant or linear: uniform valuations and precisions);
* by bivariate Kronecker substitution when, after subtracting a common
  linear trend, the finite values span a small range $R$: the minimum is
  the lowest power of $y$ present in the coefficient of $x^k$ in
  $(\sum_i x^i y^{a_i})(\sum_j x^j y^{b_j})$, which one integer
  multiplication of $O(nR \log n)$ bits provides;
* otherwise by bands of values: with the values split into 32 bands, the
  convolution is the minimum over the pairs of bands of the convolutions
  restricted to those bands, each done directly or by Kronecker
  substitution over the actual range of values of the band (small for
  clustered values: block-dependent precisions, alternating or periodic
  patterns). The pairs are processed by increasing sum of the band
  levels, so that outputs become final as soon as their current value is
  below the level of the remaining pairs, and pairs whose output range only
  contains final outputs (or outputs without any finite term) are skipped.

The cost is linear or quasi-linear in the regular cases and at most about
that of the quadratic algorithm (often much less, thanks to the early
termination) for irregular ones. With `PADIC_RADIX_FAST_BOUNDS`, when the
exact computation would cost more than a Kronecker substitution of range
64, the result may instead be a lower bound: the maximum of the convex
hull bound and the Kronecker substitution on values rounded down to
multiples of $\delta = \lceil R / 64 \rceil$ (within $2\delta$), with the
outputs near the ends computed exactly and those without finite terms
recognised exactly.

::: function
int _padic_radix_poly_mulmid_packed(padic_radix_struct * res, const padic_radix_struct * poly1, slong len1, const padic_radix_struct * poly2, slong len2, slong nlo, slong nhi, gr_ctx_t ctx)
int _padic_radix_poly_mulmid(padic_radix_struct * res, const padic_radix_struct * poly1, slong len1, const padic_radix_struct * poly2, slong len2, slong nlo, slong nhi, gr_ctx_t ctx)
int _padic_radix_poly_mullow(padic_radix_struct * res, const padic_radix_struct * poly1, slong len1, const padic_radix_struct * poly2, slong len2, slong n, gr_ctx_t ctx)

Set `res` to the coefficients $n_{lo}, \ldots, n_{hi} - 1$ (respectively
$0, \ldots, n - 1$) of the product of the polynomials. Coefficients that
are a single product are computed with `padic_radix_mul`. The unpacked
versions implement the generic methods. For both lengths at least 4, from
an estimate $S$ of the output windows in limbs (the size of the units,
plus one if the valuations are not multiples of $e$, plus the spread of
the limb frames): the mod version when $S = 1$; for $S = 2$ from a
length of about 64, for $S \le 10$ from about 24, and the radix KS
version for $S > 10$ from about 32 (these crossovers raised when the
frames are spread out); the packed version otherwise. Before choosing a
Kronecker substitution route when the frames are spread out, the windows
of 16 sampled outputs are estimated from the context precision, and the
packed version is used if the base frames inflate at least three of them
as above; if a route still gives up, the packed version is used without
trying the other. Shorter products use
the classical algorithm.
:::

## Linear algebra {#linear-algebra}

::: function
int padic_radix_mat_find_nonzero_pivot(slong * pivot_row, gr_mat_t mat, slong start_row, slong end_row, slong column, gr_ctx_t ctx)

Implements the generic method `gr_mat_find_nonzero_pivot`, used by the
generic LU decomposition, solving, inversion and reduced row echelon form
algorithms: chooses a row whose entry in the given column has minimal
valuation, so that the elimination only involves multipliers of
nonnegative valuation. Returns `GR_UNABLE` if there is no nonzero entry
but some inexact zero, and `GR_DOMAIN` if all entries are exact zeros.
:::

::: function
int padic_radix_mat_det_lu_full(padic_radix_t res, const gr_mat_t A, gr_ctx_t ctx)
int padic_radix_mat_det_levels(padic_radix_t res, const gr_mat_t A, gr_ctx_t ctx)
int padic_radix_mat_det(padic_radix_t res, const gr_mat_t A, gr_ctx_t ctx)

Set `res` to the determinant of the square matrix `A`. The *lu_full*
version uses Gaussian elimination with full pivoting on valuations (an
entry of minimal valuation in the whole remaining submatrix). Then the
multipliers and the pivot rows have valuations at least that of the pivot,
the Schur complements keep the absolute precision of the input, and the
pivots have the valuations $a_1 \le \cdots \le a_n$ of the elementary
divisors. For a matrix with entries known to $O(p^N)$ (and $a_n < N$) the
determinant is obtained to precision $N + a_1 + \cdots + a_{n-1}$, which is
optimal: it is that of the differential
$d \det(A) = \operatorname{tr}(\operatorname{adj}(A)\, dA)$. In comparison,
division-free algorithms (Berkowitz, cofactor expansion) only give
precision about $N$, losing up to $v(\det A) - a_n$ digits, and Gaussian
elimination with the first nonzero pivot loses precision growing with the
dimension for small $p$, even for generic matrices (about 19 digits at
$n = 60$ for $p = 7$). If the remaining submatrix consists of
zeros, the determinant of that block is bounded using the sums of the row
and column minima of the precisions.

The *levels* version is a blocked form of full pivoting. At the minimal
valuation $a$ of the current Schur complement $S$, a rank-revealing LU
decomposition of $(S / p^a) \bmod p$ over $\mathbb{F}_p$ selects $r$ rows
and columns such that $S_{11} = p^a U$ with $U$ invertible modulo $p$: these
$r$ pivots all have valuation $a$, as in full pivoting, and are eliminated
at once, with $\det S = \pm \det(S_{11}) \det(S_{22} - S_{21} S_{11}^{-1} S_{12})$,
using the LU decomposition of $S_{11}$, triangular solves and one matrix
product. The new Schur complement has entries of valuation at least
$a + 1$. The precision is again optimal, and the cost is that of a blocked
LU decomposition (a single level for a matrix invertible modulo $p$).

`padic_radix_mat_det` implements the generic method `gr_mat_det`: the
cofactor formula for $n \le 2$, the *lu_full* version up to a dimension
depending on the precision (40 at one limb, 24 at two limbs, 16 beyond),
and the *levels* version above.
:::

::: function
int padic_radix_mat_nonsingular_solve(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
int padic_radix_mat_inv(gr_mat_t X, const gr_mat_t A, gr_ctx_t ctx)

Set $X = A^{-1} B$ (respectively $X = A^{-1}$), returning `GR_DOMAIN` if
$A$ is singular and `GR_UNABLE` if this cannot be decided. The algorithm
is Gauss-Jordan elimination with full pivoting by valuation levels on
$[A \mid B]$: at each level the pivot block (as in
`padic_radix_mat_det_levels`) is applied to its rows, and its columns are
eliminated from all the other rows, including the pivot rows of the
previous levels, with one matrix product. This order of operations (that
of Gauss-Jordan elimination with full pivoting) gives every entry of the
solution with the optimal precision, that of the differential
$dX = A^{-1}(dB - dA\, X)$; LU decompositions with column pivots and back
substitution lose some digits when $A$ is not invertible modulo $p$.
When $A / p^a$ is invertible modulo $p$ for the minimal valuation $a$ of
its entries, all pivots have valuation $a$ and the generic LU solve (with
the pivots of `padic_radix_mat_find_nonzero_pivot`) is used directly, as it
is always with the context flag `PADIC_RADIX_FAST_LINALG`. These
functions implement the generic methods `gr_mat_nonsingular_solve` and
`gr_mat_inv`.
:::

::: function
int _padic_radix_mat_charpoly_berkowitz(padic_radix_struct * cp, const gr_mat_t mat, gr_ctx_t ctx)
int padic_radix_mat_charpoly_hessenberg(padic_radix_struct * cp, const gr_mat_t mat, gr_ctx_t ctx)
int padic_radix_mat_charpoly_lift(padic_radix_struct * cp, const gr_mat_t mat, gr_ctx_t ctx)
int padic_radix_mat_charpoly(padic_radix_struct * cp, const gr_mat_t mat, gr_ctx_t ctx)

Set `cp` to the characteristic polynomial $\chi(x) = \det(xI - A)$ of the
square matrix `mat`.

The *berkowitz* version is the division-free Berkowitz algorithm, with
the matrix-vector products done on the matrix packed once (falling back to
`_gr_mat_charpoly_berkowitz` for dimension at most 5, in contexts without
a relative precision, and for matrices with exact nonzero entries). For
integral entries known to $O(p^N)$ it gives every coefficient to
precision about $N$, which is optimal when the matrix is cyclic modulo $p$
but can lose many digits otherwise (dozens for the low coefficients of
moderately structured matrices of dimension 16).

The *hessenberg* version reduces the matrix to Hessenberg form by
similarity transformations with pivots of minimal valuation (so that the
transformations are in $\mathrm{GL}_n(\mathbb{Z}_p)$), in $O(n^3)$
operations; it loses a few digits.

The *lift* version gives the optimal precision. With the matrix scaled by
a power of $p$ to make it integral, $M$ the exact lift of the entries and
$\operatorname{adj}(xI - M) = \sum_k B_k x^k$, a perturbation $E$ with
$v(E_{ij}) \ge N_{ij}$ changes $c_k$ by $-\operatorname{tr}(B_k E)$ plus
terms of order at least two in $E$, whose coefficients (minors of
$xI - M$) are integral. Hence $c_k(A) = c_k(M) + O(p^{N_k})$ with
$N_k = \min(2 N_{\min}, \min_{ij} (N_{ij} + v((B_k)_{ji})))$, the precision of
the differential of $\chi$. The coefficients $c_k(M)$ are computed modulo
$p^{2 N_{\min}}$ by the *hessenberg* version at a raised working precision
(retried with more digits if needed, with the *berkowitz* version as a
last resort), and the $B_k$ by the recurrence $B_{n-1} = I$,
$B_{k-1} = M B_k + c_k I$, modulo one word first and modulo $p^{N_{\min}}$
if valuations at that bound matter, with $n$ matrix products. Exact
matrices give the exact characteristic polynomial (multimodular), rounded
to the context precision.

The lifted entries are reduced modulo the largest working precision, so
that entries of high valuation cost nothing. The working precision
$2(N_{\min} + s)$, where $p^s$ is the scaling, is large when an entry has
a much smaller valuation (or error term) than the others; the *lift*
version returns `GR_UNABLE` when $N_{\min} + s > 64(R + n) + 1024$, where
$R$ is the largest relative precision of an inexact entry. For exact
matrices it returns `GR_UNABLE` when the exact characteristic polynomial
would exceed the `size_limit` of the context, when the entries of $M$ have more than
$4 r + 256$ digits for a finite context precision $r$, or when they are
dominated by powers of $p$ (valuations above four times the size of the
units plus 256 digits), which Berkowitz's algorithm in
$\mathbf{Q}_p$ handles without expanding them.

`padic_radix_mat_charpoly` implements the generic method. For matrices
with integral entries of uniform precision $N$ that are cyclic modulo $p$
(tested with a random Krylov basis over $\mathbb{F}_p$), the optimal
precision is $N$ for every coefficient: the *hessenberg* version is tried
first and accepted if it attains it (as for matrices that are generic
modulo $p$), and otherwise the *berkowitz* version, which always does, is
used. Other matrices use the *lift*
version, together with Danilevsky's algorithm (keeping the more precise
result for each coefficient) where the bound $2 N_{\min}$ limits it,
and the *berkowitz* version when the *lift* version gives up. With the
context flag `PADIC_RADIX_FAST_LINALG`, the *hessenberg* version is used
for all matrices.
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
`GR_UNABLE` if the valuation or precision cannot be represented, in
which case `res` is set to $O(p^{-\mathrm{ERR\_MAX}})$ so that further
arithmetic on it stays within the representable range.
:::

::: function
int _padic_radix_add_sub_reference(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, int sub, gr_ctx_t ctx)
int _padic_radix_mul_reference(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx)

Simple reference implementations of addition ($x + y$, or $x - y$ if
`sub` is nonzero) and multiplication, used to test the optimised
functions.
:::
