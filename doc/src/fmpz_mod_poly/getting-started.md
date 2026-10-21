---
title: Getting started with fmpz_mod_poly
summary: How to set up a modulus context and polynomials over Z/nZ for a multiword modulus, how to use precomputed inverses for arithmetic modulo a fixed polynomial, and how the _f functions discover factors of a composite modulus.
order: 1
---

An `fmpz_mod_poly_t` is a polynomial over $\Z/n\Z$ whose coefficients are
`fmpz` integers in $[0, n)$. Unlike `nmod_poly_t`, it does not store its
modulus: the modulus lives in an `fmpz_mod_ctx_t` (see `fmpz_mod`), which
is created once and passed to every function, and which must outlive all
polynomials used with it. The functions follow the conventions of
`nmod_poly` (see [](nmod_poly/getting-started)): outputs may alias
inputs, underscore functions work on coefficient arrays, and the
functions that divide require invertible leading coefficients.

## A first program {#first}

The following program works modulo the Mersenne prime $p = 2^{127} - 1$.
It builds $f = (x - 3)(x - 5)(x^2 + 1)$ and computes $\gcd(f, x^p - x)$,
the product of the linear factors of $f$
([](algorithms/polynomials/root-finding)). Since $p \equiv 3 \pmod 4$, the
factor $x^2 + 1$ is irreducible and does not appear. The power $x^p \bmod
f$ is computed with `fmpz_mod_poly_powmod_x_fmpz_preinv`, which needs the
inverse of the reversal of $f$ as a power series; with it, every
reduction modulo $f$ costs two multiplications.

```c example name=fmpz_mod_poly_first
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_poly.h>

int main(void)
{
    fmpz_t p;
    fmpz_mod_ctx_t ctx;
    fmpz_mod_poly_t f, finv, t, h;

    fmpz_init(p);
    fmpz_ui_pow_ui(p, 2, 127);
    fmpz_sub_ui(p, p, 1);                       /* the prime 2^127 - 1 */
    fmpz_mod_ctx_init(ctx, p);

    fmpz_mod_poly_init(f, ctx);
    fmpz_mod_poly_init(finv, ctx);
    fmpz_mod_poly_init(t, ctx);
    fmpz_mod_poly_init(h, ctx);

    /* f = (x - 3)(x - 5)(x^2 + 1) */
    fmpz_mod_poly_set_coeff_ui(f, 1, 1, ctx);
    fmpz_mod_poly_set_coeff_si(f, 0, -3, ctx);
    fmpz_mod_poly_set_coeff_ui(t, 1, 1, ctx);
    fmpz_mod_poly_set_coeff_si(t, 0, -5, ctx);
    fmpz_mod_poly_mul(f, f, t, ctx);
    fmpz_mod_poly_zero(t, ctx);
    fmpz_mod_poly_set_coeff_ui(t, 2, 1, ctx);
    fmpz_mod_poly_set_coeff_ui(t, 0, 1, ctx);
    fmpz_mod_poly_mul(f, f, t, ctx);
    flint_printf("f = "); fmpz_mod_poly_print_pretty(f, "x", ctx); flint_printf("\n");

    /* finv = inverse of the reversal of f modulo x^len(f) */
    fmpz_mod_poly_reverse(finv, f, f->length, ctx);
    fmpz_mod_poly_inv_series(finv, finv, f->length, ctx);

    /* h = x^p mod f, then gcd(f, h - x) = product of the linear factors */
    fmpz_mod_poly_powmod_x_fmpz_preinv(h, p, f, finv, ctx);
    fmpz_mod_poly_zero(t, ctx);
    fmpz_mod_poly_set_coeff_ui(t, 1, 1, ctx);
    fmpz_mod_poly_sub(h, h, t, ctx);
    fmpz_mod_poly_gcd(h, f, h, ctx);
    flint_printf("gcd(f, x^p - x) = "); fmpz_mod_poly_print_pretty(h, "x", ctx); flint_printf("\n");

    fmpz_mod_poly_clear(f, ctx);
    fmpz_mod_poly_clear(finv, ctx);
    fmpz_mod_poly_clear(t, ctx);
    fmpz_mod_poly_clear(h, ctx);
    fmpz_mod_ctx_clear(ctx);
    fmpz_clear(p);
    return 0;
}
```
```output
f = x^4+170141183460469231731687303715884105719*x^3+16*x^2+170141183460469231731687303715884105719*x+15
gcd(f, x^p - x) = x^2+170141183460469231731687303715884105719*x+15
```

Coefficients are printed as their representatives in $[0, p)$: the
coefficient of $x$ in the gcd is $p - 8$, so the gcd is $x^2 - 8x + 15 =
(x - 3)(x - 5)$.

## Arithmetic modulo a fixed polynomial {#preinv}

Many functions come in a plain version and a `_preinv` version that
takes the inverse $\mathit{finv}$ of the reversal $\operatorname{rev}(f)
= x^{\deg f} f(1/x)$ modulo $x^{\deg f + 1}$: `fmpz_mod_poly_mulmod_preinv`,
the `powmod` functions, the division functions `_newton_n_preinv` and the
Brent–Kung modular composition. Computing $\mathit{finv}$ costs about as
much as two multiplications ([](algorithms/polynomials/newton)), so it
pays off as soon as more than a few reductions modulo the same $f$ are
needed.

## Composite moduli {#composite}

Most functions work for any modulus, but those that invert elements
(division, GCDs, inverses) fail when they meet a non-invertible leading
coefficient: the plain versions then raise an exception. The versions
with suffix `_f` instead return a nontrivial factor of $n$ in their first
argument, which is set to 1 when the computation succeeds. This is useful
when the factorisation of $n$ is unknown, as in elliptic curve methods
or in algorithms that work "as if $n$ were prime" and either succeed or
split $n$.

Here $n = 1000003 \cdot 998244353$, and the Euclidean algorithm for
$\gcd((x + 1)(x + 2), (x + 1)(1000003x + 3))$ needs to invert the leading
coefficient 1000003 of the second polynomial:

```c example name=fmpz_mod_poly_gcd_f
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_poly.h>

int main(void)
{
    fmpz_t p, q, n, f;
    fmpz_mod_ctx_t ctx;
    fmpz_mod_poly_t a, b, g;

    fmpz_init(p); fmpz_init(q); fmpz_init(n); fmpz_init(f);
    fmpz_set_ui(p, 1000003);
    fmpz_set_ui(q, 998244353);
    fmpz_mul(n, p, q);                          /* n = p q, assumed unknown */
    fmpz_mod_ctx_init(ctx, n);

    fmpz_mod_poly_init(a, ctx);
    fmpz_mod_poly_init(b, ctx);
    fmpz_mod_poly_init(g, ctx);

    /* a = (x + 1)(x + 2), b = (x + 1)(p x + 3) */
    fmpz_mod_poly_set_coeff_ui(a, 2, 1, ctx);
    fmpz_mod_poly_set_coeff_ui(a, 1, 3, ctx);
    fmpz_mod_poly_set_coeff_ui(a, 0, 2, ctx);
    fmpz_mod_poly_set_coeff_fmpz(b, 2, p, ctx);
    fmpz_add_ui(f, p, 3);
    fmpz_mod_poly_set_coeff_fmpz(b, 1, f, ctx);
    fmpz_mod_poly_set_coeff_ui(b, 0, 3, ctx);

    fmpz_mod_poly_gcd_f(f, g, a, b, ctx);
    if (fmpz_is_one(f))
    {
        flint_printf("gcd = "); fmpz_mod_poly_print_pretty(g, "x", ctx); flint_printf("\n");
    }
    else
    {
        flint_printf("factor of n found: "); fmpz_print(f); flint_printf("\n");
    }

    fmpz_mod_poly_clear(a, ctx);
    fmpz_mod_poly_clear(b, ctx);
    fmpz_mod_poly_clear(g, ctx);
    fmpz_mod_ctx_clear(ctx);
    fmpz_clear(p); fmpz_clear(q); fmpz_clear(n); fmpz_clear(f);
    return 0;
}
```
```output
factor of n found: 1000003
```

## Choosing a polynomial type {#choice}

For moduli $n < 2^{\bits}$, use `nmod_poly`, whose coefficients are words
and whose multiplication uses dedicated FFTs, which makes it considerably
faster. For larger moduli `fmpz_mod_poly` is the standard choice. The
generic `gr_poly` type over an `fmpz_mod` or `mpn_mod` context provides
the same algorithms; the `mpn_mod` context, which stores coefficients as limb arrays of fixed
size, is designed to be faster than `fmpz` coefficients for moduli of a
few words.
