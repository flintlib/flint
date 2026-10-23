---
title: fq_poly_factor.h
header: fq_poly_factor.h
template-headers: fq_poly_factor_templates
description: Factorisation, irreducibility and squarefreeness testing and root finding for polynomials over finite fields $\mathbb{F}_q$ of arbitrary characteristic.
category: Finite fields
order: 9
---

This module factors polynomials over a finite field $\mathbb{F}_q$ in the
representation of [](fq) into irreducible factors, tests irreducibility
and squarefreeness, and finds roots. All algorithms are those of the
generic module [](gr_poly), which this module wraps: squarefree
factorisation, distinct-degree factorisation, equal-degree splitting, and
the complete algorithms of Cantor and Zassenhaus [@CanZas1981], Berlekamp
[@Ber1970] and Kaltofen and Shoup [@KalSho1998]. The interface is
generated from `fq_poly_factor_templates.h` and is the same for
`fq_nmod_poly_factor` and `fq_zech_poly_factor`; `fq_default_poly_factor`
provides it for any finite field.

::: summary algorithms/polynomials/factoring-finite-fields :::

## Factorisation structure {#type}

::: type
fq_poly_factor_struct
fq_poly_factor_t

A factorisation: an array `poly` of `alloc` initialised polynomials and an
array `exp` of exponents, of which the first `num` pairs are used. The
factored polynomial is the product of `poly[i]^exp[i]` (times a leading
coefficient that is returned separately). `fq_poly_factor_t` is an array
of one structure.
:::

::: function
void fq_poly_factor_init(fq_poly_factor_t fac, const fq_ctx_t ctx)

Initialises `fac` as an empty factorisation, with room for 5 factors.
:::

::: function
void fq_poly_factor_clear(fq_poly_factor_t fac, const fq_ctx_t ctx)

Frees the memory used by `fac`.
:::

::: function
void fq_poly_factor_realloc(fq_poly_factor_t fac, slong alloc, const fq_ctx_t ctx)
void fq_poly_factor_fit_length(fq_poly_factor_t fac, slong len, const fq_ctx_t ctx)

`fq_poly_factor_realloc` changes the allocation to `alloc` factors
(clearing factors beyond it; `alloc` = 0 reinitialises `fac`).
`fq_poly_factor_fit_length` ensures room for `len` factors, at least
doubling the allocation when it grows.
:::

::: function
void fq_poly_factor_set(fq_poly_factor_t res, const fq_poly_factor_t fac, const fq_ctx_t ctx)

Sets `res` to a copy of `fac`.
:::

::: function
void fq_poly_factor_insert(fq_poly_factor_t fac, const fq_poly_t poly, slong exp, const fq_ctx_t ctx)

Multiplies the factorisation by `poly^exp`: if `poly` is already one of
the factors, `exp` is added to its exponent, otherwise `poly` is appended
with exponent `exp`. Constant polynomials are ignored.
:::

::: function
void fq_poly_factor_concat(fq_poly_factor_t res, const fq_poly_factor_t fac, const fq_ctx_t ctx)

Inserts all factors of `fac` into `res` with `fq_poly_factor_insert`.
:::

::: function
void fq_poly_factor_pow(fq_poly_factor_t fac, slong exp, const fq_ctx_t ctx)

Multiplies all exponents by `exp`, raising the factored polynomial to the
power `exp`.
:::

::: function
void fq_poly_factor_get_poly(fq_poly_t z, const fq_poly_factor_t fac, slong i, const fq_ctx_t ctx)

Sets `z` to the factor of index $i$, $0 \le i < \mathit{num}$.
:::

::: function
void fq_poly_factor_print(const fq_poly_factor_t fac, const fq_ctx_t ctx)
void fq_poly_factor_print_pretty(const fq_poly_factor_t fac, const char * var, const fq_ctx_t ctx)

Print the factors to `stdout`, one per line, each followed by ` ^ ` and
its exponent, in the format of `fq_poly_print`, respectively of
`fq_poly_print_pretty` with the variable name `var`.
:::

## Squarefree and irreducible polynomials {#tests}

::: function
int _fq_poly_is_squarefree(const fq_struct * f, slong len, const fq_ctx_t ctx)
int fq_poly_is_squarefree(const fq_poly_t f, const fq_ctx_t ctx)

Return whether $f$ is squarefree, that is, not divisible by the square of
a polynomial of positive degree. Polynomials of length at most 2 other
than 0 are squarefree, and 0 is not. This is a wrapper of
`gr_poly_is_squarefree`, which computes $\gcd(f, f')$ and handles
$p$-th powers.
:::

::: function
void fq_poly_factor_squarefree(fq_poly_factor_t res, const fq_poly_t f, const fq_ctx_t ctx)

Sets `res` to a squarefree factorisation of $f$: pairwise coprime monic
squarefree polynomials $g_i$ with distinct exponents $e_i$ such that
$f = c \prod g_i^{e_i}$, where $c$ is the leading coefficient of $f$.
Previous contents of `res` are discarded; for $f$ of length at most 1 the
factorisation is empty.
:::

::: function
int fq_poly_is_irreducible(const fq_poly_t f, const fq_ctx_t ctx)
int fq_poly_is_irreducible_ddf(const fq_poly_t f, const fq_ctx_t ctx)
int fq_poly_is_irreducible_ben_or(const fq_poly_t f, const fq_ctx_t ctx)

Return whether $f$ is irreducible. Unlike `gr_poly_is_irreducible`, these
functions return 1 for all polynomials of length at most 2, including
constants and zero.

::: details
`fq_poly_is_irreducible` calls `gr_poly_is_irreducible`, which uses the
distinct-degree test below degree 600 and Rabin's test above it.
`fq_poly_is_irreducible_ddf` uses the coarse phase of the baby-step
giant-step distinct-degree factorisation, which stops at the first factor
found, and `fq_poly_is_irreducible_ben_or` Ben-Or's test, computing
$\gcd(f, x^{q^i} - x)$ for $i \le \deg f / 2$ [@BenOr1981]. See
[](algorithms/polynomials/factoring-finite-fields).
:::
:::

## Factorisation {#factorisation}

::: function
void fq_poly_factor(fq_poly_factor_t result, fq_t leading_coeff, const fq_poly_t input, const fq_ctx_t ctx)

Factors `input` as $c \prod g_i^{e_i}$ with distinct monic irreducible
$g_i$: sets `result` to the factors and exponents (discarding previous
contents) and `leading_coeff` to $c$. For the zero polynomial $c = 0$ and
the factorisation is empty; for a constant, $c$ is the constant.

::: details
This is a wrapper of `_gr_poly_factor_finite_field` with the default
algorithm. The polynomial is first deflated (written as a polynomial in
$x^k$ for the largest possible $k$), then a squarefree factorisation is
computed, and each squarefree factor is factored with the
Cantor–Zassenhaus algorithm if its degree $n$ satisfies
$n < 125 / (\log_2 q)^{1.3}$, and with
the Kaltofen–Shoup algorithm otherwise. The factors of the deflated
polynomial are then inflated and factored again.
:::
:::

::: function
void fq_poly_factor_with_cantor_zassenhaus(fq_poly_factor_t result, fq_t leading_coeff, const fq_poly_t input, const fq_ctx_t ctx)
void fq_poly_factor_with_berlekamp(fq_poly_factor_t result, fq_t leading_coeff, const fq_poly_t input, const fq_ctx_t ctx)
void fq_poly_factor_with_kaltofen_shoup(fq_poly_factor_t result, fq_t leading_coeff, const fq_poly_t input, const fq_ctx_t ctx)

As `fq_poly_factor`, but factoring the squarefree parts with the
Cantor–Zassenhaus algorithm (distinct-degree factorisation followed by
equal-degree splitting), Berlekamp's algorithm (a basis of the kernel of
the Frobenius map minus the identity, followed by random splitting), or
the Kaltofen–Shoup algorithm (distinct-degree factorisation with
baby-step giant-step Frobenius computations by modular composition).
:::

::: function
void fq_poly_factor_cantor_zassenhaus(fq_poly_factor_t res, const fq_poly_t f, const fq_ctx_t ctx)
void fq_poly_factor_berlekamp(fq_poly_factor_t factors, const fq_poly_t f, const fq_ctx_t ctx)
void fq_poly_factor_kaltofen_shoup(fq_poly_factor_t res, const fq_poly_t poly, const fq_ctx_t ctx)

As the `_with_` functions, without returning the leading coefficient:
the output is the factorisation of the monic polynomial $f / c$. Previous
contents of the output are discarded.
:::

::: function
void fq_poly_factor_distinct_deg(fq_poly_factor_t res, const fq_poly_t poly, slong * const * degs, const fq_ctx_t ctx)

Distinct-degree factorisation of the squarefree polynomial `poly`:
appends to `res`, with exponent 1, polynomials $g_1, g_2, \ldots$ such
that $\operatorname{monic}(\mathit{poly}) = \prod g_i$ and each $g_i$ is
the product of all irreducible factors of degree $d_i$, and writes the
degrees $d_i$ to `(*degs)[0]`, `(*degs)[1]`, .... The caller must
provide room for $\deg(\mathit{poly})$ entries in `*degs`. The result is
undefined if `poly` is not squarefree. This is a wrapper of
`gr_poly_factor_distinct_deg`.
:::

::: function
void fq_poly_factor_equal_deg(fq_poly_factor_t factors, const fq_poly_t pol, slong d, const fq_ctx_t ctx)

Equal-degree factorisation: appends to `factors`, with exponent 1, the
irreducible factors of `pol`, which must be a monic squarefree product of
irreducible polynomials of degree $d$. The factors are found by random
splitting (Cantor–Zassenhaus), using an internal random state.
:::

::: function
void fq_poly_factor_split_single(fq_poly_t linfactor, const fq_poly_t input, const fq_ctx_t ctx)

Sets `linfactor` to one monic linear factor $x - a$ of `input`, which must
split into linear factors (not necessarily distinct). If `input` has
degree 1 it is copied. Throws `FLINT_ERROR` if no root is found. This is
the step used by `fq_embed_gens` to find a root of a defining
polynomial in an extension field.
:::

::: function
void fq_poly_iterated_frobenius_preinv(fq_poly_t * rop, slong n, const fq_poly_t v, const fq_poly_t vinv, const fq_ctx_t ctx)

Sets `rop[i]` to $x^{q^i} \bmod v$ for $0 \le i < n$, where `vinv` is the
inverse of the reversal of $v$ (see `fq_poly_reverse` and
`fq_poly_inv_series`). The array `rop` must contain $n$ initialised
polynomials. If `FQ_POLY_ITERATED_FROBENIUS_CUTOFF` returns 1, the
powers after the first are computed by modular composition with a
precomputed Brent–Kung matrix ($x^{q^{i}} = x^{q^{i-1}}(x^q)$), otherwise
by exponentiation.
:::

::: function
int FQ_POLY_ITERATED_FROBENIUS_CUTOFF(const fq_ctx_t ctx, slong length)

Returns 1 if `fq_poly_iterated_frobenius_preinv` should use modular
composition for a modulus of length `length`, which is the case if
$\log_2 q < 3(\lfloor \sqrt{\mathit{length}} \rfloor + 1)$.
(`FQ_NMOD_POLY_ITERATED_FROBENIUS_CUTOFF` and
`FQ_ZECH_POLY_ITERATED_FROBENIUS_CUTOFF` use $2 \log_2 q$ instead.)
:::

## Roots {#roots}

::: function
void fq_poly_roots(fq_poly_factor_t r, const fq_poly_t f, int with_multiplicity, const fq_ctx_t ctx)

Sets `r` to the roots of the nonzero polynomial $f$ in $\mathbb{F}_q$,
stored as the monic linear factors $x - a$. If `with_multiplicity` is
nonzero, the exponents are the multiplicities of the roots; otherwise
they are 1. Previous contents of `r` are discarded. This is a wrapper of
`gr_poly_roots_finite_field`, which splits $\gcd(f, x^q - x)$ by
equal-degree factorisation; see
[](algorithms/polynomials/root-finding).
:::

## Example {#example}

Factoring $x^8 - 1$ over $\mathbb{F}_3$, where it has the two roots
$\pm 1$ and three irreducible quadratic factors, and over
$\mathbb{F}_{3^2}$, whose multiplicative group of order 8 consists of
the roots of $x^8 - 1$, so that the polynomial splits into linear
factors:

```c example name=fq_poly_factor_cyclotomic
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fq.h>
#include <flint/fq_poly.h>
#include <flint/fq_poly_factor.h>

void factor_x8_minus_1(slong d)
{
    fmpz_t p;
    fq_ctx_t ctx;
    fq_poly_t f;
    fq_poly_factor_t fac;
    fq_t c;
    slong i, deg[20] = { 0 };

    fmpz_init_set_ui(p, 3);
    fq_ctx_init(ctx, p, d, "a");
    fq_poly_init(f, ctx);
    fq_poly_factor_init(fac, ctx);
    fq_init(c, ctx);

    fq_one(c, ctx);
    fq_poly_set_coeff(f, 8, c, ctx);
    fq_neg(c, c, ctx);
    fq_poly_set_coeff(f, 0, c, ctx);              /* x^8 - 1 */

    fq_poly_factor(fac, c, f, ctx);
    for (i = 0; i < fac->num; i++)
        deg[fq_poly_degree(fac->poly + i, ctx)]++;
    flint_printf("GF(3^%wd): %wd factors, of degree", d, fac->num);
    for (i = 1; i < 20; i++)
        if (deg[i] != 0)
            flint_printf(" %wd (%wd times)", i, deg[i]);
    flint_printf("\n");

    fq_poly_roots(fac, f, 1, ctx);
    flint_printf("GF(3^%wd): %wd roots\n", d, fac->num);

    fq_clear(c, ctx);
    fq_poly_factor_clear(fac, ctx);
    fq_poly_clear(f, ctx);
    fq_ctx_clear(ctx);
    fmpz_clear(p);
}

int main(void)
{
    factor_x8_minus_1(1);
    factor_x8_minus_1(2);
    return 0;
}
```
```output
GF(3^1): 5 factors, of degree 1 (2 times) 2 (3 times)
GF(3^1): 2 roots
GF(3^2): 8 factors, of degree 1 (8 times)
GF(3^2): 8 roots
```
