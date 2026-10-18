---
title: Getting started with fmpq_poly
summary: How polynomials over Q are represented as an integer numerator with a common denominator, what canonical form means, how to build and inspect polynomials, and how the underscore functions and weak canonicalisation work.
order: 1
---

An `fmpq_poly_t` represents a polynomial with rational coefficients as an
integer polynomial divided by a positive integer, for example
$$
\tfrac12 x^2 + \tfrac13 x + \tfrac16 = \frac{3x^2 + 2x + 1}{6}.
$$
The numerator is stored as an array of `fmpz` coefficients, like an
`fmpz_poly_t`, and the denominator as one `fmpz`. Arithmetic is therefore
integer polynomial arithmetic plus a little bookkeeping with integer
gcds, which is much faster than working with a vector of `fmpq`
coefficients, each with its own denominator.

## Canonical form {#canonical}

The representation is made unique by requiring that the denominator is
positive and coprime to the content of the numerator (the gcd of its
coefficients), and that the numerator is normalised (no leading zero
coefficients); the zero polynomial is $0/1$. All functions without a
leading underscore expect and produce this *canonical form*. Functions
that set individual coefficients, such as `fmpq_poly_set_coeff_fmpq`,
restore it automatically. When the numerator and denominator are
modified directly, through `fmpq_poly_numref` and `fmpq_poly_denref`,
`fmpq_poly_canonicalise` must be called afterwards.

```c example name=fmpq_poly_basics
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq_poly.h>

int main(void)
{
    fmpq_poly_t f, g, h;
    fmpq_t c;
    fmpz_poly_t num;
    fmpz_t den;

    fmpq_poly_init(f);
    fmpq_poly_init(g);
    fmpq_poly_init(h);
    fmpq_init(c);
    fmpz_poly_init(num);
    fmpz_init(den);

    /* f = x^2/2 + x/3 + 1/6, stored as (3x^2 + 2x + 1)/6 */
    fmpq_set_si(c, 1, 2);
    fmpq_poly_set_coeff_fmpq(f, 2, c);
    fmpq_set_si(c, 1, 3);
    fmpq_poly_set_coeff_fmpq(f, 1, c);
    fmpq_set_si(c, 1, 6);
    fmpq_poly_set_coeff_fmpq(f, 0, c);
    fmpq_poly_get_numerator(num, f);
    fmpq_poly_get_denominator(den, f);
    flint_printf("f = "); fmpq_poly_print_pretty(f, "x");
    flint_printf("  numerator "); fmpz_poly_print_pretty(num, "x");
    flint_printf(", denominator "); fmpz_print(den); flint_printf("\n");
    flint_printf("plain format: "); fmpq_poly_print(f); flint_printf("\n");

    /* the gcd over Q is monic */
    fmpq_poly_set_str(g, "3  -2 0 6");               /* 6x^2 - 2 */
    fmpq_poly_set_str(h, "2  1/3 1");                /* x + 1/3 */
    fmpq_poly_mul(g, g, h);
    fmpq_poly_scalar_mul_si(h, h, 3);
    fmpq_poly_mul(h, h, f);
    fmpq_poly_gcd(f, g, h);
    flint_printf("gcd = "); fmpq_poly_print_pretty(f, "x"); flint_printf("\n");

    /* power series: exp(x/2) to order x^6 */
    fmpq_poly_set_str(g, "2  0 1/2");
    fmpq_poly_exp_series(h, g, 6);
    flint_printf("exp(x/2) = "); fmpq_poly_print_pretty(h, "x"); flint_printf(" + O(x^6)\n");

    fmpq_poly_clear(f);
    fmpq_poly_clear(g);
    fmpq_poly_clear(h);
    fmpq_clear(c);
    fmpz_poly_clear(num);
    fmpz_clear(den);
    return 0;
}
```
```output
f = 1/2*x^2+1/3*x+1/6  numerator 3*x^2+2*x+1, denominator 6
plain format: 3  1/6 1/3 1/2
gcd = x+1/3
exp(x/2) = 1/3840*x^5+1/384*x^4+1/48*x^3+1/8*x^2+1/2*x+1 + O(x^6)
```

The plain format, used by `fmpq_poly_print` and `fmpq_poly_set_str`, is
the length followed by two spaces and the coefficients in increasing
degree. Over the field $\mathbb{Q}$, `fmpq_poly_gcd` returns the monic gcd.

## Underscore functions {#underscore}

The functions with a leading underscore work on a numerator array `poly`,
a denominator `den` and a length `len`, and leave memory management to the
caller. They usually require nonzero lengths, inputs in lowest terms and
without zero padding, and no aliasing between inputs and outputs, and
they produce outputs in lowest terms unless stated otherwise. Some of
them, and the non-underscore functions with suffix `_can`, take a flag
`can`: with `can = 0` the result is not reduced to lowest terms. Skipping
the reduction saves the integer gcds in a long sequence of additions,
after which a single `fmpq_poly_canonicalise` restores the canonical
form; this is called *weak canonicalisation*.

## Algorithms {#algorithms}

Since $f = a/d$ with $a \in \Z[x]$, most operations reduce to `fmpz_poly`:
a product $(a/d)(b/e)$ is $ab/(de)$ after cancelling the gcds of the
contents of $a$ and $b$ with $e$ and $d$; division with remainder is
pseudo-division over $\Z$; the gcd is the gcd of the primitive parts of
the numerators, made monic; and the resultant is computed from that of
the numerators. Power series functions such as `fmpq_poly_exp_series`
use fraction-free basecase algorithms and Newton iteration
([](algorithms/polynomials/newton)), keeping a single denominator for the
whole series.
