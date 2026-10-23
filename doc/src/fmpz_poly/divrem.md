---
title: Division and pseudo-division
functions: fmpz_poly_divrem fmpz_poly_div fmpz_poly_rem fmpz_poly_divexact fmpz_poly_divides fmpz_poly_remove fmpz_poly_pseudo_divrem fmpz_poly_pseudo_div
---

$\Z[x]$ is not a Euclidean ring: $x$ divided by $2x$ has no quotient
and remainder with $\deg R < 1$ in $\Z[x]$. FLINT therefore offers
three kinds of division. *Exact division* (`fmpz_poly_divexact`,
`fmpz_poly_divides`) computes $A/B$ when $B$ divides $A$. *Division
with remainder* (`fmpz_poly_divrem`, `fmpz_poly_div`, `fmpz_poly_rem`)
computes integer quotient coefficients from the top down and agrees with
division in $\mathbb{Q}[x]$ when the leading coefficient of $B$ is $\pm 1$
or the division is exact. *Pseudo-division* (`fmpz_poly_pseudo_divrem`)
multiplies $A$ by a power of the leading coefficient of $B$ so that
the division in $\mathbb{Q}[x]$ becomes integral. For division in
$\mathbb{Q}[x]$ itself, use `fmpq_poly_divrem`.

## Division with remainder {#divrem}

Let $\ell$ be the leading coefficient of $B$, of degree $m$. The
quotient coefficients are determined from the top: while the partial
remainder $R$ has a coefficient $r_k$ of degree $k \ge m$ (starting
from the top), the quotient coefficient of $x^{k-m}$ is $q = \lfloor
r_k / \ell \rfloor$ if $|r_k| \ge |\ell|$ and $q = 0$ otherwise, and
$R$ is replaced by $R - q x^{k-m} B$. The result satisfies

$$
A = BQ + R, \qquad |r_k| < |\ell| \text{ for all } k \ge m.
$$

If $\ell = \pm 1$ every step is exact, and $Q$ and $R$ are the
quotient and remainder of $\mathbb{Q}[x]$, with $\deg R < m$. The same
holds if $B$ divides $A$, since then all steps of the division in
$\mathbb{Q}[x]$ are integral. Otherwise $R$ may have degree $\ge m$.
The divide-and-conquer algorithm, used from $\operatorname{len}(B) \ge
6$ for `fmpz_poly_divrem`, gives the same $Q$ and $R$ as the basecase
algorithm; its cost is $O(M(n) \log n)$ coefficient operations for a
product cost $M(n)$ [@Mul2000].

For repeated division by the same polynomial with leading coefficient
$\pm 1$, precompute an inverse with `fmpz_poly_preinvert` and use
`fmpz_poly_divrem_preinv`, which replaces the division by two truncated
products ([](algorithms/polynomials/newton)), or precompute the powers
of $x$ modulo $B$ with `fmpz_poly_powers_precompute` for reducing
polynomials of length less than $2 \operatorname{len}(B)$.

## Exact division {#exact}

`fmpz_poly_divexact` assumes that $B$ divides $A$ and does not check
it. For short quotients or short divisors it computes the quotient
from both ends at once: the low half of $Q$ from the constant coefficients upwards, as a power
series quotient, and the high half from the top down. Each coefficient
is an exact integer division, and each half only involves the part of
$B$ that it needs, which halves the work of the basecase algorithm.
Long divisions use the divide-and-conquer quotient.

`fmpz_poly_divides` decides divisibility and computes the quotient.
Most non-divisors are rejected cheaply: the constant coefficient of
$B$ must divide that of $A$, $B(1)$ must divide $A(1)$, and the
division with `exact` set stops at the first quotient coefficient that
is not an integer. `fmpz_poly_remove` divides out the largest power of
a polynomial, using the values at 1 to bound the exponent first.

## Pseudo-division {#pseudo}

Pseudo-division finds $Q, R \in \Z[x]$ and $d \ge 0$ with

$$
\ell^d A = BQ + R, \qquad \deg R < \deg B.
$$

With the classical choice $d = \deg A - \deg B + 1$ the division of
$\ell^d A$ by $B$ is always integral; this is what
`fmpz_poly_pseudo_divrem_cohen` and `fmpz_poly_pseudo_rem_cohen` do,
following Algorithm 3.1.2 of [@Coh1996], as the subresultant algorithm
requires. `fmpz_poly_pseudo_divrem` multiplies by $\ell$ only when a
quotient coefficient would otherwise not be an integer, so $d$ is
smaller when the division is close to exact, and it has a
divide-and-conquer implementation for long polynomials. In both cases
$Q/\ell^d$ and $R/\ell^d$ are the quotient and remainder of $A$ by $B$
in $\mathbb{Q}[x]$.

## Examples {#examples}

Division by $2x + 1$, which is not monic: the remainder of $x^3 + x +
1$ has degree 2, while pseudo-division gives the quotient and
remainder over $\mathbb{Q}$ scaled by $2^3$. Exact division and the
divisibility test work for any divisor.

```c example name=fmpz_poly_divrem_kinds
#include <flint/flint.h>
#include <flint/fmpz_poly.h>

int main(void)
{
    fmpz_poly_t a, b, q, r;
    ulong d;

    fmpz_poly_init(a);
    fmpz_poly_init(b);
    fmpz_poly_init(q);
    fmpz_poly_init(r);

    fmpz_poly_set_str(a, "4  1 1 0 1");       /* x^3 + x + 1 */
    fmpz_poly_set_str(b, "2  1 2");           /* 2x + 1 */

    fmpz_poly_divrem(q, r, a, b);
    flint_printf("divrem:        q = "); fmpz_poly_print_pretty(q, "x");
    flint_printf(", r = "); fmpz_poly_print_pretty(r, "x"); flint_printf("\n");

    fmpz_poly_pseudo_divrem(q, r, &d, a, b);
    flint_printf("pseudo_divrem: q = "); fmpz_poly_print_pretty(q, "x");
    flint_printf(", r = "); fmpz_poly_print_pretty(r, "x");
    flint_printf(", d = %wu\n", d);

    fmpz_poly_mul(a, a, b);                   /* now b divides a */
    fmpz_poly_divexact(q, a, b);
    flint_printf("divexact:      q = "); fmpz_poly_print_pretty(q, "x"); flint_printf("\n");

    fmpz_poly_set_coeff_si(a, 0, 2);
    flint_printf("divides:       %d\n", fmpz_poly_divides(q, a, b));

    fmpz_poly_clear(a);
    fmpz_poly_clear(b);
    fmpz_poly_clear(q);
    fmpz_poly_clear(r);
    return 0;
}
```
```output
divrem:        q = 0, r = x^3+x+1
pseudo_divrem: q = 4*x^2-2*x+5, r = 3, d = 3
divexact:      q = x^3+x+1
divides:       0
```

The first quotient is zero because the leading coefficient 1 of $A$ is
smaller in absolute value than 2. Over $\mathbb{Q}$, $x^3 + x + 1 =
(2x + 1)(x^2/2 - x/4 + 5/8) + 3/8$, which is the pseudo-division
divided by $2^3$.

With a monic divisor, division with remainder is the division of
$\mathbb{Q}[x]$; the example checks $A = BQ + R$ for random
polynomials, once with the default algorithm and once with a
precomputed inverse.

```c example name=fmpz_poly_divrem_monic
#include <flint/flint.h>
#include <flint/fmpz_poly.h>

int main(void)
{
    flint_rand_t state;
    fmpz_poly_t a, b, binv, q, r, t;

    flint_rand_init(state);
    fmpz_poly_init(a);
    fmpz_poly_init(b);
    fmpz_poly_init(binv);
    fmpz_poly_init(q);
    fmpz_poly_init(r);
    fmpz_poly_init(t);

    fmpz_poly_randtest(a, state, 300, 100);
    fmpz_poly_randtest(b, state, 120, 100);
    fmpz_poly_set_coeff_ui(b, 120, 1);        /* monic of degree 120 */

    fmpz_poly_divrem(q, r, a, b);
    fmpz_poly_mul(t, b, q);
    fmpz_poly_add(t, t, r);
    flint_printf("A == BQ + R: %d, deg R = %wd\n",
                 fmpz_poly_equal(t, a), fmpz_poly_degree(r));

    fmpz_poly_preinvert(binv, b);
    fmpz_poly_divrem_preinv(t, r, a, b, binv);
    flint_printf("preinv agrees: %d\n", fmpz_poly_equal(t, q));

    fmpz_poly_clear(a);
    fmpz_poly_clear(b);
    fmpz_poly_clear(binv);
    fmpz_poly_clear(q);
    fmpz_poly_clear(r);
    fmpz_poly_clear(t);
    flint_rand_clear(state);
    return 0;
}
```
```output
A == BQ + R: 1, deg R = 119
preinv agrees: 1
```
