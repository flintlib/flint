---
title: Factoring polynomials over the integers
functions: fmpz_poly_factor fmpz_poly_factor_zassenhaus fmpz_poly_factor_squarefree
---

`fmpz_poly_factor` writes a polynomial $f \in \Z[x]$ as $f = c
\prod_i f_i^{e_i}$, where $c$ is the content of $f$ with the sign of its
leading coefficient and the $f_i$ are the distinct irreducible
primitive factors of positive degree, with positive leading
coefficients. The factors are stored in `fac->p[i]` with
multiplicities `fac->exp[i]` for $0 \le i <$ `fac->num`, and $c$ in
`fac->c`. The integer content is not factored; use `fmpz_factor` on
`fac->c` for that. `fmpz_poly_factor_squarefree` computes only the
squarefree factorisation, which is much cheaper.

::: summary algorithms/polynomials/factoring-integers :::

## How it works {#how}

After removing powers of $x$, deflating, and computing a squarefree
factorisation, each squarefree part $g$ is factored modulo a small
prime $p$; of the first three primes for which $g$ stays squarefree,
the one with fewest modular factors is used. The modular factors are
lifted to factors modulo $p^a$ by [Hensel lifting](algorithms/polynomials/hensel-lifting).
A modular factorisation into $r$ factors corresponds to at most $r$
true factors, each the product of a subset of the modular factors.

- If $r \le 8$, the Zassenhaus algorithm tries the subsets in order of
  increasing size, after lifting far enough that a true factor is
  recovered exactly from its image modulo $p^a$ (the Mignotte bound),
  and skips subsets whose degree is impossible modulo one of the three
  primes.
- If $r > 8$, van Hoeij's algorithm finds the subsets as the short
  vectors of a lattice built from the logarithmic derivatives of the
  modular factors, so that the time is polynomial in $r$ instead of
  exponential. The lifting starts at a precision much below the
  Mignotte bound and doubles until the lattice reveals the factors.

Most polynomials are irreducible, and most irreducible polynomials of
degree $n$ have about $\log n$ factors modulo a random prime, so the
typical cost is that of the modular factorisation and a short lift.
The hard cases are polynomials with many modular factors modulo every
prime, such as the Swinnerton-Dyer polynomials (factors of degree at
most 2 modulo every prime) or products of many polynomials with
abelian Galois groups; for those, van Hoeij's algorithm is essential.

## Examples {#examples}

$x^{12} - 1$ is the product of the cyclotomic polynomials $\Phi_d$
for $d \mid 12$. Deflation reduces the problem to factoring $x - 1$
and certifying or factoring the inflations. A repeated factor and a
content show up in the factorisation of $2 (x^2 - 1)^3 (x^2 + 1)$:

```c example name=fmpz_poly_factor_basic
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>

static void print_factors(const fmpz_poly_factor_t fac)
{
    slong i;

    fmpz_print(&fac->c);
    for (i = 0; i < fac->num; i++)
    {
        flint_printf(" * (");
        fmpz_poly_print_pretty(fac->p + i, "x");
        flint_printf(")^%wd", fac->exp[i]);
    }
    flint_printf("\n");
}

int main(void)
{
    fmpz_poly_t f, g;
    fmpz_poly_factor_t fac;

    fmpz_poly_init(f);
    fmpz_poly_init(g);

    fmpz_poly_set_coeff_si(f, 12, 1);
    fmpz_poly_set_coeff_si(f, 0, -1);         /* x^12 - 1 */
    fmpz_poly_factor_init(fac);
    fmpz_poly_factor(fac, f);
    flint_printf("x^12 - 1: %wd factors\n", fac->num);
    fmpz_poly_factor_clear(fac);

    fmpz_poly_set_str(f, "3  -1 0 1");        /* x^2 - 1 */
    fmpz_poly_pow(f, f, 3);
    fmpz_poly_set_str(g, "3  2 0 2");         /* 2x^2 + 2 */
    fmpz_poly_mul(f, f, g);

    fmpz_poly_factor_init(fac);
    fmpz_poly_factor(fac, f);
    print_factors(fac);
    fmpz_poly_factor_clear(fac);

    fmpz_poly_factor_init(fac);
    fmpz_poly_factor_squarefree(fac, f);
    print_factors(fac);
    fmpz_poly_factor_clear(fac);

    fmpz_poly_clear(f);
    fmpz_poly_clear(g);
    return 0;
}
```
```output
x^12 - 1: 6 factors
2 * (x^2+1)^1 * (x-1)^3 * (x+1)^3
2 * (x^2+1)^1 * (x^2-1)^3
```

The Swinnerton-Dyer polynomial $S_4$ of degree 16 is irreducible but
splits into eight quadratic factors modulo every prime, and the
product $S_4 \cdot S_4(x + 1)$ of degree 32 has at least 16 modular factors, so
that the Zassenhaus algorithm would have to try tens of thousands of subsets;
`fmpz_poly_factor` uses van Hoeij's algorithm and finds the two
factors at once.

```c example name=fmpz_poly_factor_swinnerton_dyer
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>

int main(void)
{
    fmpz_poly_t f, g;
    fmpz_poly_factor_t fac;
    fmpz_t one;

    fmpz_poly_init(f);
    fmpz_poly_init(g);
    fmpz_init_set_ui(one, 1);

    fmpz_poly_swinnerton_dyer(f, 4);
    fmpz_poly_factor_init(fac);
    fmpz_poly_factor(fac, f);
    flint_printf("S_4: degree %wd, %wd factor(s)\n", fmpz_poly_degree(f), fac->num);
    fmpz_poly_factor_clear(fac);

    fmpz_poly_taylor_shift(g, f, one);
    fmpz_poly_mul(f, f, g);
    fmpz_poly_factor_init(fac);
    fmpz_poly_factor(fac, f);
    flint_printf("S_4(x) S_4(x+1): degree %wd, %wd factors of degrees %wd and %wd\n",
        fmpz_poly_degree(f), fac->num,
        fmpz_poly_degree(fac->p + 0), fmpz_poly_degree(fac->p + 1));
    fmpz_poly_factor_clear(fac);

    fmpz_poly_clear(f);
    fmpz_poly_clear(g);
    fmpz_clear(one);
    return 0;
}
```
```output
S_4: degree 16, 1 factor(s)
S_4(x) S_4(x+1): degree 32, 2 factors of degrees 16 and 16
```
