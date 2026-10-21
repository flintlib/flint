---
title: gr_ore_poly.h
header: gr_ore_poly.h
description: Dense univariate Ore polynomials (linear differential, difference, q-difference, Mahler and Frobenius operators) over generic rings, as a generic noncommutative ring.
category: Generic rings
order: 10
---

This module implements Ore polynomials $L = \sum_i a_i D^i$ with
coefficients in a ring $R$ given by a generic context
[@Ore1933; @BroPet1996]. The multiplication is twisted by an endomorphism
$\sigma$ of $R$ and a $\sigma$-derivation $\delta$ (an additive map with
$\delta(ab) = \sigma(a)\delta(b) + \delta(a) b$) through the commutation
rule

$$
D a = \sigma(a) D + \delta(a), \qquad a \in R.
$$

With $R = K[x]$ and $\sigma = \mathrm{id}$, $\delta = d/dx$, Ore
polynomials are linear differential operators; with $\sigma(x) = x +
1$, $\delta = 0$, linear recurrence (shift) operators; and so on for
the predefined algebras listed below. The context of an Ore polynomial
ring is a `gr_ctx_t`, so Ore polynomials can be used with the generic
functions (`gr_mul` composes operators) and in other generic
structures.

The module is under development: it provides memory management,
additive arithmetic, multiplication and Euclidean division.

## Example {#example}

The program below multiplies differential operators over $\mathbb{Q}[x]$,
where $D x = x D + 1$, and recurrence operators, where $S x = (x + 1)
S$. Coefficients are written on the left of the powers of the
generator.

```c example name=gr_ore_poly_operators
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_poly.h>
#include <flint/gr_ore_poly.h>

int main(void)
{
    gr_ctx_t QQ, QQx, Dx, Sx;
    gr_ore_poly_t L, M, P;

    gr_ctx_init_fmpq(QQ);
    gr_ctx_init_gr_poly(QQx, QQ);           /* base ring Q[x] */

    /* differential operators: D x = x D + 1 */
    gr_ore_poly_ctx_init(Dx, QQx, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(L, Dx);
    gr_ore_poly_init(M, Dx);
    gr_ore_poly_init(P, Dx);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(L, "D", Dx));
    GR_MUST_SUCCEED(gr_ore_poly_set_str(M, "x", Dx));
    GR_MUST_SUCCEED(gr_ore_poly_mul(P, L, M, Dx));
    flint_printf("D * x         = "); gr_ore_poly_print(P, Dx); flint_printf("\n");
    GR_MUST_SUCCEED(gr_ore_poly_set_str(L, "x*D^2 + D - x", Dx));
    GR_MUST_SUCCEED(gr_ore_poly_set_str(M, "D + x^2", Dx));
    GR_MUST_SUCCEED(gr_ore_poly_mul(P, L, M, Dx));
    flint_printf("L * M         = "); gr_ore_poly_print(P, Dx); flint_printf("\n");
    gr_ore_poly_clear(L, Dx);
    gr_ore_poly_clear(M, Dx);
    gr_ore_poly_clear(P, Dx);
    gr_ore_poly_ctx_clear(Dx);

    /* recurrence operators: S x = (x + 1) S */
    gr_ore_poly_ctx_init(Sx, QQx, 0, ORE_ALGEBRA_FORWARD_SHIFT);
    GR_MUST_SUCCEED(_gr_ore_poly_ctx_set_gen_name(Sx, "S"));
    gr_ore_poly_init(L, Sx);
    gr_ore_poly_init(P, Sx);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(L, "S - x - 1", Sx));
    GR_MUST_SUCCEED(gr_ore_poly_mul(P, L, L, Sx));
    flint_printf("(S - x - 1)^2 = "); gr_ore_poly_print(P, Sx); flint_printf("\n");
    gr_ore_poly_clear(L, Sx);
    gr_ore_poly_clear(P, Sx);
    gr_ore_poly_ctx_clear(Sx);

    gr_ctx_clear(QQx);
    gr_ctx_clear(QQ);
    return 0;
}
```
```output
D * x         = 1 + x*D
L * M         = (4*x - x^3) + (-x + 5*x^2)*D + (1 + x^3)*D^2 + x*D^3
(S - x - 1)^2 = (1 + 2*x + x^2) + (-3 - 2*x)*S + S^2
```

## Ore algebras {#algebras}

::: type
ore_algebra_t

The predefined Ore algebra structures, an enumeration with the values
below. Except for the commutative and Frobenius algebras, they act on a
distinguished generator $x$ of the base ring, which must currently be a
univariate polynomial ring (created with `gr_ctx_init_gr_poly`);
otherwise $\sigma$ and $\delta$ return `GR_UNABLE`.

- `ORE_ALGEBRA_CUSTOM`: $\sigma$ and $\delta$ given by a user function
  (`gr_ore_poly_ctx_init_custom`).
- `ORE_ALGEBRA_COMMUTATIVE`: $\sigma = \mathrm{id}$, $\delta = 0$;
  ordinary polynomials in $D$.
- `ORE_ALGEBRA_DERIVATIVE`: $\sigma = \mathrm{id}$, $\delta = d/dx$;
  linear differential operators.
- `ORE_ALGEBRA_EULER_DERIVATIVE`: $\sigma = \mathrm{id}$, $\delta = x
  \, d/dx$ (the Euler operator $\theta$).
- `ORE_ALGEBRA_FORWARD_SHIFT`: $\sigma(a)(x) = a(x + 1)$, $\delta = 0$;
  linear recurrence operators.
- `ORE_ALGEBRA_FORWARD_DIFFERENCE`: $\sigma(a)(x) = a(x + 1)$,
  $\delta(a) = a(x + 1) - a(x)$ (the operator $\Delta$).
- `ORE_ALGEBRA_BACKWARD_SHIFT`: $\sigma(a)(x) = a(x - 1)$, $\delta =
  0$.
- `ORE_ALGEBRA_BACKWARD_DIFFERENCE`: $\sigma(a)(x) = a(x - 1)$,
  $\delta(a) = a(x) - a(x - 1)$ (the operator $\nabla$).
- `ORE_ALGEBRA_Q_SHIFT`: $\sigma(a)(x) = a(qx)$, $\delta = 0$; linear
  $q$-difference operators (`gr_ore_poly_ctx_init_q_shift`).
- `ORE_ALGEBRA_MAHLER`: $\sigma(a)(x) = a(x^b)$, $\delta = 0$; linear
  Mahler operators with base $b$ (`gr_ore_poly_ctx_init_mahler`).
- `ORE_ALGEBRA_FROBENIUS`: $\sigma$ the Frobenius automorphism $a
  \mapsto a^p$ of a finite field ($\delta = 0$), computed with
  `gr_fq_frobenius`.

The value `ORE_POLY_NUM_ALGEBRAS` is the number of algebra types.
:::

::: function
ore_algebra_t ore_algebra_randtest(flint_rand_t state)

Returns a random Ore algebra type, favouring the derivative and the
forward shift, for testing.
:::

## Types and conventions {#conventions}

::: type
gr_ore_poly_struct
gr_ore_poly_t

An Ore polynomial: an array `coeffs` of `alloc` coefficients in the base
ring, of which the first `length` are used; the coefficient of $D^i$ is
at index $i$. The layout is that of `gr_poly_t`, so the functions of
`gr_poly` that do not multiply (copying, truncation, addition, scalar
operations, and the container functions) can be applied to the
coefficient arrays with the base ring context.
:::

::: type
gr_ore_poly_ctx_struct
gr_ore_poly_ctx_t

The context of an Ore polynomial ring, a `gr_ctx_t`. Its data (type
`_gr_ore_poly_ctx_struct`, compatible with the context data of
`gr_poly` rings) holds the base ring, a degree limit, the name of the
generator (`D` by default), the algebra type, the function computing
$\sigma$ and $\delta$, and algebra-specific data (type
`gr_ore_poly_ore_data_t`: the index of the distinguished generator
`base_var`, $\sigma(x)$ and $q$ for $q$-shifts, or the Mahler base).
:::

::: macro
GR_ORE_POLY_CTX(ctx)
GR_ORE_POLY_ELEM_CTX(ctx)
GR_ORE_POLY_ORE_DATA(ctx)

The context data, the base ring, and the algebra-specific data of an
Ore polynomial context.
:::

The context argument of the functions is the context of the Ore
polynomial ring (not the base ring as in `gr_poly`). As in `gr_poly`, Ore
polynomials are weakly normalised (a leading coefficient is removed
only if it is provably zero, so the length bounds the order plus one
over inexact rings), functions returning `int` return status codes, and
the functions with a leading underscore work on coefficient arrays
without memory management.

## Contexts {#contexts}

::: function
void gr_ore_poly_ctx_init(gr_ore_poly_ctx_t ctx, gr_ctx_t base_ring, slong base_var, const ore_algebra_t which_algebra)
void gr_ctx_init_gr_ore_poly(gr_ctx_t ctx, gr_ctx_t base_ring, slong base_var, const ore_algebra_t which_algebra)

Initialise `ctx` as the ring of Ore polynomials over `base_ring` with
the predefined algebra `which_algebra`, acting on the generator with
index `base_var` of the base ring. These are for the algebras without
parameters; the $q$-shift and Mahler algebras have their own
constructors. The two functions are the same.
:::

::: function
int gr_ore_poly_ctx_init_q_shift(gr_ore_poly_ctx_t ctx, gr_ctx_t base_ring, slong base_var, gr_srcptr q)
int gr_ore_poly_ctx_init_mahler(gr_ore_poly_ctx_t ctx, gr_ctx_t base_ring, slong base_var, slong mahler_base)

Initialise `ctx` for the $q$-shift algebra $\sigma(a)(x) = a(qx)$ with
$q$ an element of the base ring (`GR_UNABLE` unless the base ring is a
univariate polynomial ring), or the Mahler algebra $\sigma(a)(x) =
a(x^b)$ with $b$ = `mahler_base` (`GR_DOMAIN` for $b = 0$).
:::

::: function
void gr_ore_poly_ctx_init_custom(gr_ore_poly_ctx_t ctx, gr_ctx_t base_ring, const gr_ore_poly_sigma_delta_t sigma_delta, void * ore_data)
void * gr_ore_poly_ctx_data_ptr(gr_ore_poly_ctx_t ctx)

`gr_ore_poly_ctx_init_custom` initialises `ctx` with a user-defined
algebra: `sigma_delta` computes $\sigma$ and $\delta$ (see
`gr_ore_poly_sigma_delta`), and `ore_data` is a pointer to data for it,
which `sigma_delta` can retrieve with `gr_ore_poly_ctx_data_ptr(ctx)`.
The data is not freed by `gr_ore_poly_ctx_clear`.
:::

::: function
void gr_ore_poly_ctx_init_randtest(gr_ore_poly_ctx_t ctx, flint_rand_t state, gr_ctx_t base_ring)
void gr_ore_poly_ctx_init_randtest2(gr_ctx_t base_ring, gr_ore_poly_ctx_t ctx, flint_rand_t state)

Initialise `ctx` with a random algebra over `base_ring`, or (the second
function) over a random base ring that it initialises in `base_ring`,
for testing.
:::

::: function
void gr_ore_poly_ctx_clear(gr_ore_poly_ctx_t ctx)

Clears the context and the algebra data it owns.
:::

::: function
int _gr_ore_poly_ctx_set_gen_name(gr_ctx_t ctx, const char * s)
int _gr_ore_poly_ctx_set_gen_names(gr_ctx_t ctx, const char ** s)
int gr_ore_poly_gens_recursive(gr_vec_t vec, gr_ore_poly_ctx_t ctx)

Set the name of the generator $D$ to `s` (respectively `s[0]`), and set
`vec` to the generators of the base ring (as constant Ore polynomials)
followed by $D$.
:::

::: function
int gr_ore_poly_ctx_write(gr_stream_t out, gr_ore_poly_ctx_t ctx)
truth_t gr_ore_poly_ctx_is_ring(gr_ore_poly_ctx_t ctx)
truth_t gr_ore_poly_ctx_is_zero_ring(gr_ore_poly_ctx_t ctx)
truth_t gr_ore_poly_ctx_is_commutative_ring(gr_ore_poly_ctx_t ctx)
truth_t gr_ore_poly_ctx_is_integral_domain(gr_ore_poly_ctx_t ctx)
truth_t gr_ore_poly_ctx_is_threadsafe(gr_ore_poly_ctx_t ctx)

Context methods. Whether the ring is a ring, the zero ring or
thread-safe is inherited from the base ring; commutativity and the
absence of zero divisors are currently reported as `T_UNKNOWN`.
:::

## Memory management {#memory}

::: function
void gr_ore_poly_init(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
void gr_ore_poly_init2(gr_ore_poly_t poly, slong len, gr_ore_poly_ctx_t ctx)
void gr_ore_poly_clear(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)

Initialise `poly` as zero (with space for `len` coefficients), and
clear it.
:::

::: function
gr_ptr gr_ore_poly_coeff_ptr(gr_ore_poly_t poly, slong i, gr_ore_poly_ctx_t ctx)
gr_srcptr gr_ore_poly_coeff_srcptr(const gr_ore_poly_t poly, slong i, gr_ore_poly_ctx_t ctx)
slong gr_ore_poly_length(const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)

Return a pointer to the coefficient of $D^i$, and the length.
:::

::: function
void gr_ore_poly_swap(gr_ore_poly_t poly1, gr_ore_poly_t poly2, gr_ore_poly_ctx_t ctx)
void gr_ore_poly_set_shallow(gr_ore_poly_t res, const gr_ore_poly_t x, const gr_ore_poly_ctx_t ctx)

Swap two Ore polynomials in constant time, and make `res` a shallow copy
of `x` sharing its coefficients.
:::

::: function
void gr_ore_poly_fit_length(gr_ore_poly_t poly, slong len, gr_ore_poly_ctx_t ctx)
void _gr_ore_poly_set_length(gr_ore_poly_t poly, slong len, gr_ore_poly_ctx_t ctx)
void _gr_ore_poly_normalise(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)

Ensure space for `len` coefficients, set the length without
normalising, and remove provably zero leading coefficients, as the
`gr_poly` functions of the same names.
:::

## Basic manipulation {#basic}

::: function
int gr_ore_poly_set(gr_ore_poly_t res, const gr_ore_poly_t src, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_truncate(gr_ore_poly_t poly, const gr_ore_poly_t src, slong newlen, gr_ore_poly_ctx_t ctx)

Set `res` to a copy of `src`, or to `src` with the terms of order
`newlen` and higher removed.
:::

::: function
int gr_ore_poly_zero(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_one(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_neg_one(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_gen(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)

Set `poly` to $0$, $1$, $-1$ or the generator $D$.
:::

::: function
int gr_ore_poly_set_si(gr_ore_poly_t poly, slong x, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_set_ui(gr_ore_poly_t poly, ulong x, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_set_fmpz(gr_ore_poly_t poly, const fmpz_t x, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_set_fmpq(gr_ore_poly_t poly, const fmpq_t x, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_set_other(gr_ore_poly_t poly, gr_srcptr x, gr_ctx_t x_ctx, gr_ore_poly_ctx_t ctx)

Set `poly` to a constant converted into the base ring. The function
`gr_ore_poly_set_other` accepts an Ore polynomial with the same context
(copied), an element of the base ring, an element of the coefficient
ring of the base ring when the base ring is a univariate polynomial
ring, or an integer or rational number (`fmpz` and `fmpq` contexts),
each giving a constant; it returns `GR_UNABLE` for any other ring.
:::

::: function
truth_t _gr_ore_poly_equal(gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ore_poly_ctx_t ctx)
truth_t gr_ore_poly_equal(const gr_ore_poly_t poly1, const gr_ore_poly_t poly2, gr_ore_poly_ctx_t ctx)
truth_t gr_ore_poly_is_zero(const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
truth_t gr_ore_poly_is_one(const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
truth_t gr_ore_poly_is_gen(const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)

Comparisons, decided coefficientwise in the base ring. The underscore
version requires $\mathit{len1} \ge \mathit{len2}$.
:::

## Input and output {#io}

::: function
int _gr_ore_poly_write(gr_stream_t out, gr_srcptr poly, slong len, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_write(gr_stream_t out, const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
int _gr_ore_poly_get_str(char ** res, gr_srcptr f, slong len, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_get_str(char ** res, const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_print(const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)

Write the Ore polynomial as an expression in the generator name of the
context, with the coefficients on the left of the powers of $D$, to a
stream, to a new string (to be freed with `flint_free`), or to standard
output.
:::

::: function
int _gr_ore_poly_set_str(gr_ptr res, const char * s, slong len, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_set_str(gr_ore_poly_t res, const char * s, gr_ore_poly_ctx_t ctx)

Parse an Ore polynomial from an expression in the generator name of the
context and the generator names of the base ring; products are
evaluated in the Ore algebra, so `D*x` and `x*D` differ. The underscore
version writes exactly `len` coefficients (zero-padding) and returns
`GR_UNABLE` if the result is longer; intermediate results may be
longer.
:::

::: function
int gr_ore_poly_randtest(gr_ore_poly_t poly, flint_rand_t state, slong len, gr_ore_poly_ctx_t ctx)
int _gr_ore_poly_randtest_default(gr_ore_poly_t res, flint_rand_t state, gr_ore_poly_ctx_t ctx)

Set `poly` to a random Ore polynomial of length at most `len`
(respectively a random length below 5), for testing.
:::

## Action on the base ring {#action}

::: type
gr_ore_poly_sigma_delta_t

The type of a function `int f(gr_ptr sigma, gr_ptr delta, gr_srcptr a,
gr_ore_poly_ctx_struct * ctx)` that sets `sigma` to $\sigma(a)$ and
`delta` to $\delta(a)$ for an element $a$ of the base ring; either
output may be `NULL`, in which case it is not computed.
:::

::: function
int gr_ore_poly_sigma_delta(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_sigma(gr_ptr res, gr_srcptr a, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_delta(gr_ptr res, gr_srcptr a, gr_ore_poly_ctx_t ctx)

Compute $\sigma(a)$ and $\delta(a)$, or one of them, by calling the
function stored in the context.
:::

::: function
int sigma_delta_unable(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_struct * ctx)
int sigma_delta_commutative(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_struct * ctx)
int sigma_delta_derivative(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_struct * ctx)
int sigma_delta_euler_derivative(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_struct * ctx)
int sigma_delta_forward_shift(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_struct * ctx)
int sigma_delta_backward_shift(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_struct * ctx)
int sigma_delta_forward_difference(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_struct * ctx)
int sigma_delta_backward_difference(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_struct * ctx)
int sigma_delta_compose(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_struct * ctx)
int sigma_delta_mahler(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_struct * ctx)
int sigma_delta_frobenius(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_struct * ctx)

The functions implementing $\sigma$ and $\delta$ for the predefined
algebras, with the signature `gr_ore_poly_sigma_delta_t`.
`sigma_delta_unable` always returns `GR_UNABLE`;
`sigma_delta_commutative` is the identity and zero; the derivative and
Euler derivative versions differentiate a polynomial (with
`gr_poly_derivative`, then a multiplication by $x$); the shift and
difference versions compute Taylor shifts $a(x \pm 1)$ (with
`gr_poly_taylor_shift`) and their differences with $a$;
`sigma_delta_compose` (used for $q$-shifts) composes $a$ with the stored
$\sigma(x) = qx$; `sigma_delta_mahler` substitutes $x^b$ (an inflation);
and `sigma_delta_frobenius` applies `gr_fq_frobenius`. All except the
commutative and Frobenius versions require a base ring of univariate
polynomials and return `GR_UNABLE` otherwise.
:::

::: variable
_gr_ore_poly_default_sigma_delta

The array of the `sigma_delta` functions of the predefined algebras,
indexed by `ore_algebra_t` (`sigma_delta_unable` for
`ORE_ALGEBRA_CUSTOM`).
:::

## Arithmetic {#arithmetic}

::: function
int gr_ore_poly_neg(gr_ore_poly_t res, const gr_ore_poly_t src, gr_ore_poly_ctx_t ctx)
int _gr_ore_poly_add(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_add(gr_ore_poly_t res, const gr_ore_poly_t poly1, const gr_ore_poly_t poly2, gr_ore_poly_ctx_t ctx)
int _gr_ore_poly_sub(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_sub(gr_ore_poly_t res, const gr_ore_poly_t poly1, const gr_ore_poly_t poly2, gr_ore_poly_ctx_t ctx)

Negation, sum and difference, coefficientwise.
:::

::: function
int gr_ore_poly_add_ui(gr_ore_poly_t res, const gr_ore_poly_t poly, ulong c, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_add_si(gr_ore_poly_t res, const gr_ore_poly_t poly, slong c, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_add_fmpz(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpz_t c, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_add_fmpq(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpq_t c, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_add_other(gr_ore_poly_t res, const gr_ore_poly_t poly, gr_srcptr x, gr_ctx_t x_ctx, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_sub_ui(gr_ore_poly_t res, const gr_ore_poly_t poly, ulong c, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_sub_si(gr_ore_poly_t res, const gr_ore_poly_t poly, slong c, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_sub_fmpz(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpz_t c, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_sub_fmpq(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpq_t c, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_sub_other(gr_ore_poly_t res, const gr_ore_poly_t poly, gr_srcptr x, gr_ctx_t x_ctx, gr_ore_poly_ctx_t ctx)

Add or subtract a constant $c$. The `_other` versions accept the same
kinds of elements as `gr_ore_poly_set_other`, including an Ore
polynomial with the same context.
:::

::: function
int gr_ore_poly_mul_ui(gr_ore_poly_t res, const gr_ore_poly_t poly, ulong c, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_mul_si(gr_ore_poly_t res, const gr_ore_poly_t poly, slong c, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_mul_fmpz(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpz_t c, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_mul_fmpq(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpq_t c, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_other_mul(gr_ore_poly_t res, gr_srcptr x, gr_ctx_t x_ctx, const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)

Multiply by an integer or rational constant (which commutes with $D$),
or multiply on the left by $x$, that is compute $x L$: an Ore
polynomial with the same context is multiplied with `gr_ore_poly_mul`,
and the other kinds of elements accepted by `gr_ore_poly_set_other` are
constants $c$ with $c L = \sum_i (c a_i) D^i$.
:::

::: function
int _gr_ore_poly_lmul_gen(gr_ptr res, gr_srcptr poly, slong len, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_lmul_gen(gr_ore_poly_t res, const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)

Set `res` to the left product $D L = \sum_i (\sigma(a_i) D^{i+1} +
\delta(a_i) D^i)$. The underscore version requires $\mathit{len} \ge 1$
and writes $\mathit{len} + 1$ coefficients.
:::

::: function
int _gr_ore_poly_mul(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_mul(gr_ore_poly_t res, const gr_ore_poly_t poly1, const gr_ore_poly_t poly2, gr_ore_poly_ctx_t ctx)

Set `res` to the product $L_1 L_2$ in the Ore algebra, computed as
$\sum_i a_i (D^i L_2)$ with the powers $D^i L_2$ obtained by repeated
left multiplication by $D$: $O(\mathit{len1} \cdot \mathit{len2})$
applications of $\sigma$ and $\delta$ and multiplications in the base
ring. The underscore version requires positive lengths, writes
$\mathit{len1} + \mathit{len2} - 1$ coefficients and does not allow
aliasing.
:::

::: function
int _gr_ore_poly_divrem(gr_ptr Q, gr_ptr R, gr_srcptr U, slong lenU, gr_srcptr V, slong lenV, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_divrem(gr_ore_poly_t Q, gr_ore_poly_t R, const gr_ore_poly_t U, gr_ore_poly_t V, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_div(gr_ore_poly_t Q, const gr_ore_poly_t U, gr_ore_poly_t V, gr_ore_poly_ctx_t ctx)
int gr_ore_poly_rem(gr_ore_poly_t R, const gr_ore_poly_t U, gr_ore_poly_t V, gr_ore_poly_ctx_t ctx)

Compute the right Euclidean division $U = QV + R$ with
$\operatorname{ord} R < \operatorname{ord} V$, which exists and is
unique when the leading coefficient of $V$ and its images under the
powers of $\sigma$ are invertible (over a field with $\sigma$ injective,
for instance). Each step divides the leading coefficient of the current
remainder by $\sigma^k(\operatorname{lc}(V))$ and subtracts $c D^k V$,
for $k$ from $\operatorname{ord} U - \operatorname{ord} V$ down to 0;
the eliminated coefficients of $R$ are then set to exact zeros (in
inexact rings they would otherwise only be approximately zero).
The non-underscore functions return `GR_DOMAIN` if $V = 0$; the
underscore version requires $\mathit{lenU} \ge \mathit{lenV} \ge 1$ and
writes $\mathit{lenU} - \mathit{lenV} + 1$ coefficients of $Q$ and
$\mathit{lenU}$ of $R$. (In FLINT 3.6 and
earlier the last step, $k = 0$, was skipped, so that for example
$D^2 + D + 1$ divided by $D$ gave $Q = D$ and $R = 1 + D$.)
:::
