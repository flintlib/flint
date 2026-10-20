---
title: Using ca
summary: How to write programs with exact real and complex numbers: numbers as elements of automatically constructed fields, special values, contexts and memory, arithmetic, printing and conversion, three-valued predicates, polynomials and matrices, and the ca contexts of the generic rings interface.
order: 1
---

This guide introduces the programming interface of [](ca), FLINT's
exact real and complex numbers (formerly the library Calcium
[@Joh2021calcium]). What can and cannot be decided about these numbers
is explained in
[](algorithms/exact/zero-testing#in-practice), and the choice between
`ca_t` and the other exact representations in
[](algorithms/exact/overview).

## Numbers as elements of fields {#fields}

A `ca_t` represents a number as an element of a field

$$
K = \mathbb{Q}(a_1, \ldots, a_n)
$$

obtained by adjoining *extension numbers* $a_k$ to the rational numbers.
An extension number is either an algebraic number, given by its minimal
polynomial and an isolating ball ([](qqbar)), or a symbolic expression
such as $\pi$, $e^{z}$, $\log z$, $\sqrt{z}$ or $\Gamma(z)$ whose
arguments are again `ca_t` numbers, possibly in other fields. The fields
are constructed automatically: adding $\sqrt{2}$ and $\pi$ creates the
field $\mathbb{Q}(\sqrt 2, \pi)$, and the sum is the element $a + b$
with $a = \pi$, $b = \sqrt 2$.

To compute in $K$, FLINT uses the isomorphism

$$
\mathbb{Q}(a_1, \ldots, a_n) \cong \operatorname{Frac}(\mathbb{Q}[X_1, \ldots, X_n] / I)
$$

with the *formal field* on the right, where $I$ is the ideal of all
polynomial relations between the extension numbers. The relations may
be algebraic, as $b^2 - 2 = 0$ for $b = \sqrt 2$, or involve
transcendental numbers, as $e^{\pi} \cdot e^{-\pi} = 1$. An element is
stored as a fraction of multivariate polynomials, reduced modulo the
relations that are known; with a single algebraic generator, $K$ is a
number field $\mathbb{Q}[X]/\langle f \rangle$ and arithmetic is
fast number field arithmetic ([](nf_elem)). The connection between the
formal field and the actual complex numbers is maintained by evaluating
the generators in ball arithmetic when needed. The construction is
described in [](algorithms/exact/field-towers).

The following program shows the representation: a number is printed as
a numerical approximation followed by its formal representation and the
definitions of the generators, and relations between the generators make
results simplify automatically. The last two lines show special values.

```c example name=ca_using_fields
#include <flint/ca.h>

int main(void)
{
    ca_ctx_t ctx;
    ca_t x, y, z;

    ca_ctx_init(ctx);
    ca_init(x, ctx);
    ca_init(y, ctx);
    ca_init(z, ctx);

    /* x = sqrt(2) + pi, an element of Q(sqrt(2), pi) */
    ca_sqrt_ui(x, 2, ctx);
    ca_pi(y, ctx);
    ca_add(x, x, y, ctx);
    ca_print(x, ctx); flint_printf("\n");

    /* x^2 - 2 sqrt(2) pi - pi^2 = 2 */
    ca_mul(z, x, x, ctx);
    ca_sqrt_ui(x, 2, ctx);
    ca_mul(x, x, y, ctx);
    ca_mul_ui(x, x, 2, ctx);
    ca_sub(z, z, x, ctx);
    ca_mul(x, y, y, ctx);
    ca_sub(z, z, x, ctx);
    ca_print(z, ctx); flint_printf("\n");

    /* exp(pi) * exp(-pi) = 1 */
    ca_exp(x, y, ctx);
    ca_neg(y, y, ctx);
    ca_exp(y, y, ctx);
    ca_mul(z, x, y, ctx);
    ca_print(z, ctx); flint_printf("\n");

    /* 1 / 0 and log(0) */
    ca_zero(x, ctx);
    ca_inv(y, x, ctx);
    ca_print(y, ctx); flint_printf("\n");
    ca_log(y, x, ctx);
    ca_print(y, ctx); flint_printf("\n");

    ca_clear(x, ctx);
    ca_clear(y, ctx);
    ca_clear(z, ctx);
    ca_ctx_clear(ctx);
    return 0;
}
```

```output
4.55581 {a+b where a = 3.14159 [Pi], b = 1.41421 [b^2-2=0]}
2
1
UnsignedInfinity
-Infinity
```

## Special values {#special-values}

Every operation on `ca_t` numbers returns a value, even when the
mathematical result is not a complex number. Besides numbers, a `ca_t`
can hold the signed infinities $c \cdot \infty$ (for a complex sign $c$,
for example $+\infty$ and $-\infty$), the unsigned infinity $\tilde
\infty$, the value *Undefined* (for example $0/0$ or $\exp(\tilde\infty)$),
and the value *Unknown*, which means that the result is some value,
possibly special, that could not be determined. *Unknown* typically
arises when an operation depends on a zero test that fails: if $x$
cannot be decided to be zero or nonzero, $1/x$ is *Unknown*. The rules
are described in [](ca#special-values).

## Contexts and memory {#contexts}

All `ca_t` numbers belong to a context `ca_ctx_t`. The context stores
the extension numbers and fields that have been created, so that every
occurrence of $\pi$ or $\sqrt 2$ refers to the same object, and caches
their numerical enclosures and the relations between them. Every
function takes the context as its last argument, and numbers must be
initialised and cleared with the context they belong to:

```c
ca_ctx_t ctx;
ca_t x;

ca_ctx_init(ctx);
ca_init(x, ctx);
/* ... */
ca_clear(x, ctx);
ca_ctx_clear(ctx);
```

All numbers must be cleared before their context. Numbers of different
contexts cannot be mixed; `ca_transfer` copies a number to another
context. Since even reading a number may modify the context (for
example by refining a cached enclosure), a context must not be shared
between threads. The context also holds options that bound the work done
by zero tests and simplifications (see [](ca/options)).

## Arithmetic, printing and conversion {#arithmetic}

Numbers are created with assignment functions (`ca_set_si`,
`ca_set_fmpq`, `ca_set_qqbar`, `ca_pi`, `ca_i`, ...) and combined with
arithmetic and elementary functions (`ca_add`, `ca_div`, `ca_sqrt`,
`ca_exp`, `ca_log`, `ca_sin`, `ca_atan`, ...); aliasing of inputs and
output is allowed. `ca_print` shows a numerical approximation, the
formal representation and, in braces, the generators with their
definitions; brackets give either the minimal polynomial of an
algebraic generator or the symbolic definition of a transcendental one.
For output in other forms there are `ca_get_decimal_str` (a decimal
approximation with a given number of digits, computed to that accuracy)
and `ca_get_fexpr` (a symbolic expression, see [](fexpr)). `ca_get_acb`
gives a ball enclosure at any precision, and `ca_get_qqbar`,
`ca_get_fmpq` and `ca_get_fmpz` convert to exact types when the number
can be proved to belong to them.

Predicates return a `truth_t` value, `T_TRUE`, `T_FALSE` or `T_UNKNOWN`
(see [](calcium#truth) and [](ca/check-is-zero)); functions such as
`ca_floor` and `ca_sgn` use them internally and return *Unknown* if they
fail.

```c example name=ca_using_basics
#include <flint/calcium.h>
#include <flint/ca.h>
#include <flint/qqbar.h>
#include <flint/fexpr.h>

int main(void)
{
    ca_ctx_t ctx;
    ca_t phi, x, y;
    qqbar_t a;
    fexpr_t e;
    char * s;

    ca_ctx_init(ctx);
    ca_init(phi, ctx);
    ca_init(x, ctx);
    ca_init(y, ctx);
    qqbar_init(a);
    fexpr_init(e);

    /* phi = (1 + sqrt(5)) / 2 */
    ca_sqrt_ui(phi, 5, ctx);
    ca_add_ui(phi, phi, 1, ctx);
    ca_div_ui(phi, phi, 2, ctx);
    flint_printf("phi = "); ca_print(phi, ctx); flint_printf("\n");

    /* phi^2 - phi - 1 is simplified to 0 */
    ca_mul(x, phi, phi, ctx);
    ca_sub(x, x, phi, ctx);
    ca_sub_ui(x, x, 1, ctx);
    flint_printf("phi^2 - phi - 1 = "); ca_print(x, ctx); flint_printf("\n");

    /* conversions */
    s = ca_get_decimal_str(phi, 40, 0, ctx);
    flint_printf("40 digits: %s\n", s);
    flint_free(s);
    if (ca_get_qqbar(a, phi, ctx))
    {
        flint_printf("qqbar: "); qqbar_print(a); flint_printf("\n");
    }
    ca_get_fexpr(e, phi, 0, ctx);
    flint_printf("fexpr: "); fexpr_print(e); flint_printf("\n");

    /* predicates return truth_t */
    ca_pi(x, ctx);
    ca_sqr(y, phi, ctx);
    ca_add_ui(y, y, 1, ctx);       /* phi^2 + 1 = phi + 2 = 3.618... */
    flint_printf("pi < phi^2 + 1: "); truth_print(ca_check_lt(x, y, ctx)); flint_printf("\n");
    ca_exp(x, x, ctx);
    flint_printf("exp(pi) is real: "); truth_print(ca_check_is_real(x, ctx)); flint_printf("\n");
    ca_floor(y, x, ctx);
    flint_printf("floor(exp(pi)) = "); ca_print(y, ctx); flint_printf("\n");

    fexpr_clear(e);
    qqbar_clear(a);
    ca_clear(phi, ctx);
    ca_clear(x, ctx);
    ca_clear(y, ctx);
    ca_ctx_clear(ctx);
    return 0;
}
```

```output
phi = 1.61803 {(a+1)/2 where a = 2.23607 [a^2-5=0]}
phi^2 - phi - 1 = 0
40 digits: 1.618033988749894848204586834365638117720
qqbar: deg 2 [-1, -1, 1] [1.6180339887498948482045868343656381177 +/- 3.69e-38]
fexpr: Where(Div(Add(a_1, 1), 2), Def(a_1, Sqrt(5)))
pi < phi^2 + 1: T_TRUE
exp(pi) is real: T_TRUE
floor(exp(pi)) = 23
```

`ca.h` includes `fexpr.h`, but `qqbar.h` declares its conversions to
and from [](fexpr), such as `qqbar_get_fexpr_formula`, only when
`fexpr.h` is included before it.

## Polynomials and matrices {#poly-mat}

[](ca_vec), [](ca_poly) and [](ca_mat) provide vectors, univariate
polynomials and matrices with `ca_t` entries, with the same context
argument. Functions whose result depends on zero tests report failure:
`ca_poly_gcd`, `ca_poly_roots` and `ca_mat_rank` return 0, and
`ca_mat_inv` (which returns `T_FALSE` for a singular matrix),
`ca_mat_check_is_zero` and similar predicates return `T_UNKNOWN`. Roots of polynomials are found
when the polynomial factors into pieces of degree at most 3 or with
rational coefficients (see [](ca_poly/roots)); eigenvalues, Jordan
forms and matrix exponentials build on them (see
[](algorithms/exact/linear-algebra)).

## The generic rings interface {#gr}

The generic rings interface [](gr) provides fields of exact numbers
with `ca_t` elements: `gr_ctx_init_real_ca` ($\mathbb{R}$),
`gr_ctx_init_complex_ca` ($\mathbb{C}$), `gr_ctx_init_real_algebraic_ca`
and `gr_ctx_init_complex_algebraic_ca` (the algebraic numbers) and
`gr_ctx_init_complex_extended_ca` ($\mathbb{C}$ with the special values)
(see [](gr_domains#real-complex)). Each such context owns a `ca_ctx_t`.
Generic predicates return `truth_t` values as above, and functions return
`GR_DOMAIN` when the result is not in the domain (for example $\sqrt{-1}$
in the real numbers) and `GR_UNABLE` when a zero test fails. All generic
polynomial, series and matrix algorithms then work with exact numbers.

```c example name=ca_using_gr
#include <flint/calcium.h>
#include <flint/gr.h>

int main(void)
{
    gr_ctx_t R;
    gr_ptr x, y;
    int status;

    gr_ctx_init_real_ca(R);
    GR_TMP_INIT2(x, y, R);

    /* sqrt(2) sqrt(3) - sqrt(6) */
    status = gr_set_ui(x, 2, R);
    status |= gr_sqrt(x, x, R);
    status |= gr_set_ui(y, 3, R);
    status |= gr_sqrt(y, y, R);
    status |= gr_mul(x, x, y, R);
    status |= gr_set_ui(y, 6, R);
    status |= gr_sqrt(y, y, R);
    status |= gr_sub(x, x, y, R);
    flint_printf("status %d, value ", status);
    status |= gr_println(x, R);
    flint_printf("is zero: "); truth_print(gr_is_zero(x, R)); flint_printf("\n");

    /* sqrt(-1) is not a real number */
    status = gr_set_si(x, -1, R);
    status |= gr_sqrt(y, x, R);
    flint_printf("sqrt(-1): status %d (GR_DOMAIN = %d)\n", status, GR_DOMAIN);

    GR_TMP_CLEAR2(x, y, R);
    gr_ctx_clear(R);
    return 0;
}
```

```output
status 0, value 0
is zero: T_TRUE
sqrt(-1): status 1 (GR_DOMAIN = 1)
```
