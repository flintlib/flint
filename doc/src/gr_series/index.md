---
title: gr_series.h
header: gr_series.h
description: Generic rings of formal power series over any base ring - truncated power series modulo x^n, and power series with an error term O(x^n) tracked per element.
category: Generic rings
order: 8
---

This module constructs two generic rings of formal power series over a
base ring $R$ given by a `gr_ctx_t`, whose elements can be used with the
whole `gr` interface (`gr_add`, `gr_mul`, `gr_exp`, `gr_gamma`, and so
on), in matrices, as coefficients of polynomials, and in any other
generic algorithm:

- **Truncated power series** `gr_series_mod`: the ring $R[[x]] / (x^n)$
  for a fixed $n$. Elements are polynomials (`gr_poly_t`) of length at
  most $n$ and are exact: the ring is a quotient ring in which $x^n = 0$.
- **Power series** `gr_series`: the ring $R[[x]]$. An element is a
  polynomial together with an error term $O(x^m)$, which stands for an
  unknown tail $c_m x^m + c_{m+1} x^{m+1} + \cdots$; the exponent $m$ is
  stored per element, and polynomials of low degree can be exact (no
  error term). The context has a default precision: tails from that
  order on are replaced by an error term.

These are different algebraic structures. $\mathbb{Q}[[x]]$ is an integral
domain, but $\mathbb{Q}[[x]] / (x^n)$ is not for $n \ge 2$. In $R[[x]] / (x^n)$
every element is known exactly, so $0 = 0$ is decided; in $R[[x]]$ the
element $0 + O(x^n)$ might hide a nonzero tail, so $0 + O(x^n) = 0$ and
even $0 + O(x^n) = 0 + O(x^n)$ are `T_UNKNOWN`. Following the enclosure
semantics of inexact rings ([](algorithms/generic/overview)), every
operation on `gr_series` elements returns a result whose error term
covers all the values the inputs can represent. When either ring would
do, `gr_series_mod` has less overhead.

The arithmetic is done by the power series functions of `gr_poly`
(`gr_poly_mullow`, `gr_poly_inv_series`, `gr_poly_exp_series` and so on)
over the base ring, so their algorithms and the ring-specific
multiplication of the base ring apply. The guide
[](gr_series/getting-started) shows both rings in a short program.

## Contexts {#contexts}

::: function
void gr_series_mod_ctx_init(gr_ctx_t ctx, gr_ctx_t base_ring, slong n)

Initialises `ctx` as the ring $R[[x]] / (x^n)$ of truncated power series
over `base_ring`. Elements have type `gr_poly_t`. All inputs are assumed
to have length at most $n$, and all outputs are truncated to length
$n$. The truncation order is fixed; elements of a context with order $m
\le n$ can be used as inputs in a context with order $n$, but not the
other way round. For $n = 0$ this is the zero ring.
:::

::: function
void gr_series_ctx_init(gr_ctx_t ctx, gr_ctx_t base_ring, slong prec)

Initialises `ctx` as the ring $R[[x]]$ of power series over `base_ring`.
Elements have type `gr_series_t`. Each element carries its own error
term $O(x^m)$; the default precision `prec` is the order from which
tails are replaced by an error term in the results of operations. The
default precision can be changed (`GR_SERIES_PREC`), and elements can
be shared between contexts with the same base ring and different
precisions.
:::

The contexts are cleared with `gr_ctx_clear`. The generator is named
`x` unless changed with `gr_ctx_set_gen_name`. The base ring must
outlive the context.

## Truncated power series {#mod}

::: macro
GR_SERIES_MOD_ELEM_CTX(ctx)
GR_SERIES_MOD_N(ctx)

The context of the base ring, and the truncation order $n$, of a
`gr_series_mod` context.
:::

::: type
series_mod_ctx_t

The data of a `gr_series_mod` context, stored in the context object:
the base ring `base_ring`, the order `n` and the name `var` of the
generator.
:::

The following functions implement the methods of the `gr_series_mod`
context. They are normally reached through the generic functions
(`gr_series_mod_add` is called by `gr_add`, and so on), but can be
called directly with the series context.

::: function
void gr_series_mod_ctx_clear(gr_ctx_t ctx)
int gr_series_mod_ctx_write(gr_stream_t out, gr_ctx_t ctx)
int gr_series_mod_ctx_set_gen_name(gr_ctx_t ctx, const char * s)
int gr_series_mod_ctx_set_gen_names(gr_ctx_t ctx, const char ** s)
int gr_series_mod_gens_recursive(gr_vec_t vec, gr_ctx_t ctx)

Context methods: clearing, printing a description, naming the
generator, and listing the generators of the tower of rings (the
generator $x$ after those of the base ring).
:::

::: function
truth_t gr_series_mod_ctx_is_ring(gr_ctx_t ctx)
truth_t gr_series_mod_ctx_is_commutative_ring(gr_ctx_t ctx)
truth_t gr_series_mod_ctx_is_integral_domain(gr_ctx_t ctx)
truth_t gr_series_mod_ctx_is_field(gr_ctx_t ctx)
truth_t gr_series_mod_ctx_is_rational_vector_space(gr_ctx_t ctx)
truth_t gr_series_mod_ctx_is_real_vector_space(gr_ctx_t ctx)
truth_t gr_series_mod_ctx_is_complex_vector_space(gr_ctx_t ctx)
truth_t gr_series_mod_ctx_is_finite_characteristic(gr_ctx_t ctx)

Structural predicates. The ring inherits the properties of the base
ring, except that it is an integral domain or a field only for $n = 1$
(when it is isomorphic to $R$), and that for $n = 0$ it is the zero
ring.
:::

::: function
void gr_series_mod_init(gr_poly_t res, gr_ctx_t ctx)
void gr_series_mod_clear(gr_poly_t res, gr_ctx_t ctx)
void gr_series_mod_swap(gr_poly_t x, gr_poly_t y, gr_ctx_t ctx)
void gr_series_mod_set_shallow(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_randtest(gr_poly_t res, flint_rand_t state, gr_ctx_t ctx)
int gr_series_mod_write(gr_stream_t out, const gr_poly_t x, gr_ctx_t ctx)

Element methods: memory management, a random element for testing, and
output in the form `f (mod x^n)`.
:::

::: function
int gr_series_mod_zero(gr_poly_t res, gr_ctx_t ctx)
int gr_series_mod_one(gr_poly_t res, gr_ctx_t ctx)
int gr_series_mod_gen(gr_poly_t res, gr_ctx_t ctx)
int gr_series_mod_set(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_set_other(gr_poly_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
int gr_series_mod_set_si(gr_poly_t res, slong c, gr_ctx_t ctx)
int gr_series_mod_set_ui(gr_poly_t res, ulong c, gr_ctx_t ctx)
int gr_series_mod_set_fmpz(gr_poly_t res, const fmpz_t c, gr_ctx_t ctx)
int gr_series_mod_set_fmpq(gr_poly_t res, const fmpq_t c, gr_ctx_t ctx)

Assignment: $0$, $1$, the generator $x$ (which is $0$ for $n \le 1$), a
copy, and conversions into the ring (with truncation to length $n$).
:::

Error terms are created with `gr_big_o_base_fmpz` and the series
context: for a constant base it sets the result to the constant
$O(\mathit{base}^{\mathit{exp}})$ of the base ring (for example a ball
$[0 \pm 2^{-e}]$), and for a nonconstant base it returns `GR_UNABLE`,
since elements of $R[[x]]/(x^n)$ have no error terms in $x$.

::: function
truth_t gr_series_mod_is_zero(const gr_poly_t x, gr_ctx_t ctx)
truth_t gr_series_mod_is_one(const gr_poly_t x, gr_ctx_t ctx)
truth_t gr_series_mod_equal(const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx)

Predicates, decided coefficientwise in the base ring.
:::

::: function
int gr_series_mod_neg(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_add(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx)
int gr_series_mod_sub(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx)
int gr_series_mod_mul(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx)
int gr_series_mod_inv(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_div(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx)

Arithmetic modulo $x^n$: products by `gr_poly_mullow`, inverses and
quotients by `gr_poly_inv_series` and `gr_poly_div_series`, which
require an invertible constant term of the divisor (`GR_DOMAIN` if it is
zero or not invertible).
:::

::: function
int gr_series_mod_sqrt(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_rsqrt(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_exp(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_log(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_sin(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_cos(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_sin_pi(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_cos_pi(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_sin_cos(gr_poly_t res1, gr_poly_t res2, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_sin_cos_pi(gr_poly_t res1, gr_poly_t res2, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_tan(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_tanh(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_cot(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_coth(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_tan_pi(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_cot_pi(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_asin(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_acos(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_atan(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_asinh(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_acosh(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_mod_atanh(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)

Elementary functions modulo $x^n$, computed with the corresponding
`gr_poly` series functions ([](gr_poly#special)). The value at the
constant term must exist in the base ring, and the coefficients involve
divisions by small integers, so these functions fail in small positive
characteristic.
:::

## Power series {#series}

::: type
gr_series_struct
gr_series_t

A power series with an error term: the structure contains a
`gr_poly_struct` `poly` (the known coefficients) and a `slong` `error`,
the exponent $m$ of the error term $O(x^m)$, or `GR_SERIES_ERR_EXACT` if
the series is exact. The polynomial has length at most $m$.
:::

::: type
gr_series_vec_struct
gr_series_vec_t

A vector of `gr_series_struct`, with the layout of a `gr_vec_t` whose
elements are power series, used for the parameters of
`gr_series_hypgeom_pfq`.
:::

::: type
series_ctx_t

The data of a `gr_series` context: the base ring `base_ring`, the
default precision `prec` and the name `var` of the generator.
:::

::: macro
GR_SERIES_POLY(x)
GR_SERIES_ERROR(x)

The polynomial part and the error exponent of a `gr_series_t`.
:::

::: macro
GR_SERIES_ERR_EXACT
GR_SERIES_ERR_MAX

`GR_SERIES_ERR_EXACT` (`WORD_MAX`) is the error exponent of an exact
series. `GR_SERIES_ERR_MAX` (`WORD_MAX / 4`) is the largest exponent
allowed in an error term $O(x^m)$.
:::

::: macro
GR_SERIES_CTX(ctx)
GR_SERIES_ELEM_CTX(ctx)
GR_SERIES_PREC(ctx)
GR_SERIES_MOD_CTX(ctx)

The context data of a `gr_series` context, its base ring and its
default precision (which may be assigned to), and the context data of a
`gr_series_mod` context.
:::

### Error terms {#error}

::: function
slong _gr_series_get_error(const gr_series_t x, gr_ctx_t ctx)

Returns the exponent $m$ of the error term $O(x^m)$ of `x`, or
`GR_SERIES_ERR_EXACT` if `x` is exact.
:::

::: function
truth_t _gr_series_is_exact(const gr_series_t x, gr_ctx_t ctx)

Returns whether `x` has no error term. This does not check whether the
coefficients themselves are exact.
:::

::: function
void _gr_series_set_error(gr_series_t x, slong err, gr_ctx_t ctx)

Adds an error term $O(x^m)$ to `x` in place, discarding the
coefficients of order $m$ and higher. The exponent is clamped to $0 \le
m \le$ `GR_SERIES_ERR_MAX`.
:::

::: function
void _gr_series_make_exact(gr_series_t x, gr_ctx_t ctx)

Removes the error term of `x`, keeping its coefficients.
:::

### Methods {#methods}

The following functions implement the methods of the `gr_series`
context and are normally reached through the generic functions. In
all of them, the error term of the result is the smallest of the error
terms of the inputs and the default precision, whenever the exact
result would be longer than the default precision; for example the sum
$1 + x + O(x^3)$ and $x^2 + O(x^5)$ is $1 + x + x^2 + O(x^3)$.

::: function
void gr_series_ctx_clear(gr_ctx_t ctx)
int gr_series_ctx_write(gr_stream_t out, gr_ctx_t ctx)
int gr_series_ctx_set_gen_name(gr_ctx_t ctx, const char * s)
int gr_series_ctx_set_gen_names(gr_ctx_t ctx, const char ** s)
int gr_series_gens_recursive(gr_vec_t vec, gr_ctx_t ctx)

Context methods: clearing, printing a description, naming the
generator, and listing the generators of the tower of rings.
:::

::: function
truth_t gr_series_ctx_is_ring(gr_ctx_t ctx)
truth_t gr_series_ctx_is_commutative_ring(gr_ctx_t ctx)
truth_t gr_series_ctx_is_approx_commutative_ring(gr_ctx_t ctx)
truth_t gr_series_ctx_is_integral_domain(gr_ctx_t ctx)
truth_t gr_series_ctx_is_rational_vector_space(gr_ctx_t ctx)
truth_t gr_series_ctx_is_real_vector_space(gr_ctx_t ctx)
truth_t gr_series_ctx_is_complex_vector_space(gr_ctx_t ctx)
truth_t gr_series_ctx_is_finite_characteristic(gr_ctx_t ctx)

Structural predicates, inherited from the base ring: $R[[x]]$ is a
ring, a commutative ring or an integral domain exactly when $R$ is.
:::

::: function
void gr_series_init(gr_series_t res, gr_ctx_t ctx)
void gr_series_clear(gr_series_t res, gr_ctx_t ctx)
void gr_series_swap(gr_series_t x, gr_series_t y, gr_ctx_t ctx)
void gr_series_set_shallow(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_randtest(gr_series_t res, flint_rand_t state, gr_ctx_t ctx)
int gr_series_write(gr_stream_t out, const gr_series_t x, gr_ctx_t ctx)

Element methods: memory management (an initialised series is the exact
zero), a random element for testing, and output in the form `f +
O(x^m)`.
:::

::: function
int gr_series_zero(gr_series_t res, gr_ctx_t ctx)
int gr_series_one(gr_series_t res, gr_ctx_t ctx)
int gr_series_gen(gr_series_t res, gr_ctx_t ctx)
int gr_series_set(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_set_gr_poly(gr_series_t res, const gr_poly_t x, gr_ctx_t ctx)
int gr_series_set_scalar(gr_series_t res, gr_srcptr x, gr_ctx_t ctx)
int gr_series_set_si(gr_series_t res, slong c, gr_ctx_t ctx)
int gr_series_set_ui(gr_series_t res, ulong c, gr_ctx_t ctx)
int gr_series_set_fmpz(gr_series_t res, const fmpz_t c, gr_ctx_t ctx)
int gr_series_set_fmpq(gr_series_t res, const fmpq_t c, gr_ctx_t ctx)
int gr_series_set_other(gr_series_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)

Assignment of exact values (a polynomial longer than the default
precision is truncated and given an error term), of the generator, and
of copies. `gr_series_set_other` converts elements of the base ring,
polynomials and power series over other rings in a variable of the same
name (keeping their error terms), and anything else that converts into
the base ring, as a constant.
:::

Error terms are created with `gr_big_o_base_fmpz` and the series
context, which sets the result to $O(\mathit{base}^{\mathit{exp}})$, for
example $O(x^{10})$ for the base $x$ and the exponent 10.

::: function
int gr_series_set_interval_mid_rad(gr_series_t res, const gr_series_t m, const gr_series_t r, gr_ctx_t ctx)

Sets `res` to the series whose coefficients are the enclosures $[m_i \pm
r_i]$ (`gr_set_interval_mid_rad` in the base ring) and whose error
term is the smaller of those of $m$ and $r$.
:::

::: function
truth_t gr_series_is_zero(const gr_series_t x, gr_ctx_t ctx)
truth_t gr_series_is_one(const gr_series_t x, gr_ctx_t ctx)
truth_t gr_series_equal(const gr_series_t x, const gr_series_t y, gr_ctx_t ctx)
truth_t gr_series_coeff_is_zero(const gr_series_t x, slong i, gr_ctx_t ctx)

Predicates. A known coefficient that differs decides `T_FALSE`, but a
series with an error term can never be shown to be zero, one or equal
to another series: the answer is then `T_UNKNOWN`. A coefficient of
order beyond the error term is `T_UNKNOWN`.
:::

::: function
int gr_series_neg(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_add(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx)
int gr_series_sub(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx)
int gr_series_mul(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx)

Ring operations. The product of an exact zero with anything is the
exact zero; otherwise the error term of a product is the smaller of the
error terms of the factors (and of the default precision, if the exact
product is longer).
:::

::: function
int gr_series_inv(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_div(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx)
int gr_series_divexact(gr_series_t res, const gr_series_t x, const gr_series_t y, gr_ctx_t ctx)

Inverse and quotients in $R[[x]]$. The inverse needs an invertible
constant term: it returns `GR_DOMAIN` for the exact zero and `GR_UNABLE`
for a series whose constant term is hidden in the error term. Division
first removes the common power of $x$ from numerator and denominator (a
quotient such as $x / x = 1$ exists in $R[[x]]$), returning `GR_DOMAIN`
if the denominator has larger valuation than the numerator and
`GR_UNABLE` if the valuations cannot be decided.
:::

::: function
int gr_series_sqrt(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_rsqrt(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_exp(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_log(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_sin(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_cos(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_sin_pi(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_cos_pi(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_sin_cos(gr_series_t res1, gr_series_t res2, const gr_series_t x, gr_ctx_t ctx)
int gr_series_sin_cos_pi(gr_series_t res1, gr_series_t res2, const gr_series_t x, gr_ctx_t ctx)
int gr_series_tan(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_tanh(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_cot(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_coth(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_tan_pi(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_cot_pi(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_asin(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_acos(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_atan(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_asinh(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_acosh(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_atanh(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)

Elementary functions of power series, computed by the `gr_poly` series
functions to the precision allowed by the error term of the input and
the default precision. If the constant term is unknown (error term
$O(x^0)$), `GR_UNABLE` is returned.
:::

::: function
int gr_series_gamma(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_rgamma(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_lgamma(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_digamma(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_erf(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_erfc(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_erfi(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_exp_integral_ei(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_cos_integral(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_cosh_integral(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_sin_integral(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_sinh_integral(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_fresnel(gr_series_t res1, gr_series_t res2, const gr_series_t x, int normalized, gr_ctx_t ctx)
int gr_series_fresnel_s(gr_series_t res, const gr_series_t x, int normalized, gr_ctx_t ctx)
int gr_series_fresnel_c(gr_series_t res, const gr_series_t x, int normalized, gr_ctx_t ctx)
int gr_series_bessel_j(gr_series_t res, const gr_srcptr nu, const gr_series_t z, gr_ctx_t ctx)
int gr_series_airy(gr_series_t res1, gr_series_t res2, gr_series_t res3, gr_series_t res4, const gr_series_t x, gr_ctx_t ctx)
int gr_series_airy_ai(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_airy_ai_prime(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_airy_bi(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_airy_bi_prime(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_log_integral(gr_series_t res, const gr_series_t x, int offset, gr_ctx_t ctx)
int gr_series_gamma_upper(gr_series_t res, const gr_series_t s, const gr_series_t x, int regularized, gr_ctx_t ctx)
int gr_series_gamma_lower(gr_series_t res, const gr_series_t s, const gr_series_t x, int regularized, gr_ctx_t ctx)
int gr_series_beta_lower(gr_series_t res, const gr_series_t a, const gr_series_t b, const gr_series_t x, int regularized, gr_ctx_t ctx)
int gr_series_polylog(gr_series_t res, const gr_series_t s, const gr_series_t z, gr_ctx_t ctx)
int gr_series_hurwitz_zeta(gr_series_t res, const gr_series_t s, const gr_series_t z, gr_ctx_t ctx)
int gr_series_dirichlet_l(gr_series_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const gr_series_t x, gr_ctx_t ctx)
int gr_series_dirichlet_hardy_theta(gr_series_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const gr_series_t x, gr_ctx_t ctx)
int gr_series_dirichlet_hardy_z(gr_series_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const gr_series_t x, gr_ctx_t ctx)
int gr_series_jacobi_theta(gr_series_t res1, gr_series_t res2, gr_series_t res3, gr_series_t res4, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx)
int gr_series_jacobi_theta_1(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx)
int gr_series_jacobi_theta_2(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx)
int gr_series_jacobi_theta_3(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx)
int gr_series_jacobi_theta_4(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx)
int gr_series_agm1(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_elliptic_k(gr_series_t res, const gr_series_t x, gr_ctx_t ctx)
int gr_series_weierstrass_p(gr_series_t res, const gr_series_t x, const gr_series_t tau, gr_ctx_t ctx)
int gr_series_hypgeom_pfq(gr_series_t res, const gr_series_vec_t a, const gr_series_vec_t b, const gr_series_t x, int regularized, gr_ctx_t ctx)

Special functions of power series. These are currently implemented
only over the base rings `arb` and `acb`, through the series functions
of `arb_hypgeom`, `acb_hypgeom`, `acb_dirichlet`, `acb_elliptic` and
`acb_modular` (for example `arb_hypgeom_gamma_series`), and return
`GR_UNABLE` for other base rings. The flags `normalized`, `offset` and
`regularized` have the meanings of the corresponding `arb_hypgeom`
functions. For the functions of two or more series arguments, all
arguments except one must currently be constant series.
:::
