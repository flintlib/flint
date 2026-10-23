---
title: Power series
functions: fmpz_poly_inv_series fmpz_poly_inv_series_basecase fmpz_poly_div_series fmpz_poly_div_series_basecase fmpz_poly_compose_series fmpz_poly_compose_series_horner fmpz_poly_revert_series fmpz_poly_sqrt_series
---

A polynomial of length at most $n$ also represents a power series in
$\Z[[x]]$ known modulo $x^n$. The functions below compute
inverses, quotients, compositions, reversions and square roots of
such series. Over $\Z$ these operations exist only under conditions
on the leading terms: $1/Q$ requires $Q(0) = \pm 1$, $A/B$ requires
the quotient to be integral (for example $B(0) = \pm 1$), reversion
requires $Q(0) = 0$ and $Q'(0) = \pm 1$, and a square root must have
integer coefficients. Truncated products (`fmpz_poly_mullow`) and
powers (`fmpz_poly_pow_trunc`) are in the multiplication and powering
sections. Transcendental functions of series with integer
coefficients, such as $\exp$ and $\log$, leave $\Z[[x]]$; they are in
`fmpq_poly`.

::: summary algorithms/polynomials/newton :::

## Algorithms {#algorithms}

*Inversion* uses the recurrence $b_0 = q_0$, $b_k = -q_0 \sum_{j=1}^k
q_j b_{k-j}$ (since $q_0^{-1} = q_0$), which costs $O(n^2)$ coefficient
operations but runs in word arithmetic as long as the coefficients
stay small, and Newton iteration $b \mapsto b + b(1 - Qb)$, which
doubles the precision with two truncated products and costs a small
multiple of one product of length $n$. `fmpz_poly_inv_series` uses
Newton iteration when both $Q$ and $n$ are at least 64. Binomials
$q_0 + q_k x^k$ are inverted directly.

*Division* $A/B$ uses the analogous recurrence for short inputs, a
multiplication by the inverse of $B$ for long inputs with $B(0) = \pm
1$, and otherwise divide-and-conquer division of the reversed
polynomials, which also works when $B(0) \ne \pm 1$ but all divisions
by $B(0)$ are exact. An inexact division raises an exception.

*Composition* $g(h) \bmod x^n$ with $h(0) = 0$ uses Horner's rule for
outer polynomials of length at most 10, the baby-step giant-step
algorithm of Brent and Kung [@BrentKung1978] for $n < 500$, and the
algorithm of Kinoshita and Li [@KL2024], which reduces composition to
bivariate products and costs $O(M(n) \log n)$, for larger $n$.
*Reversion* uses fast Lagrange inversion (baby-step giant-step) for
$n < 2500$, and Newton iteration,
whose compositions use the method of [@KL2024], above.

The *square root* of a series with even valuation $2v$ is computed
from that of the series divided by $x^{2v}$, and it is determined only
to precision $n - v$. The square root is computed as the square root
with remainder of the reversed polynomial, divide and conquer, and
exists only if the resulting coefficients are integers.

## Examples {#examples}

The inverse of $1 - x - x^2$ generates the Fibonacci numbers, the
reversion of $x - x^2$ generates the Catalan numbers, and $\sqrt{1 -
4x} = 1 - 2 \sum_{k \ge 1} C_{k-1} x^k$. A quotient with $B(0) = 3$ is
fine when it is integral.

```c example name=fmpz_poly_series_basic
#include <flint/flint.h>
#include <flint/fmpz_poly.h>

int main(void)
{
    fmpz_poly_t a, b, c;

    fmpz_poly_init(a);
    fmpz_poly_init(b);
    fmpz_poly_init(c);

    fmpz_poly_set_str(a, "3  1 -1 -1");       /* 1 - x - x^2 */
    fmpz_poly_inv_series(b, a, 10);
    flint_printf("1/(1 - x - x^2) = "); fmpz_poly_print(b); flint_printf("\n");

    fmpz_poly_set_str(a, "3  0 1 -1");        /* x - x^2 */
    fmpz_poly_revert_series(b, a, 10);
    flint_printf("reversion       = "); fmpz_poly_print(b); flint_printf("\n");

    fmpz_poly_set_str(a, "2  1 -4");          /* 1 - 4x */
    flint_printf("sqrt exists: %d, ", fmpz_poly_sqrt_series(b, a, 8));
    fmpz_poly_print(b); flint_printf("\n");

    fmpz_poly_set_str(a, "2  3 6");           /* 3 + 6x */
    fmpz_poly_set_str(c, "2  3 -3");          /* 3 - 3x */
    fmpz_poly_div_series(b, a, c, 6);
    flint_printf("(3 + 6x)/(3 - 3x) = "); fmpz_poly_print(b); flint_printf("\n");

    /* 1/(1 - y) at y = x + x^2 */
    fmpz_poly_set_str(a, "5  1 1 1 1 1");     /* 1/(1 - y) mod y^5 */
    fmpz_poly_set_str(c, "3  0 1 1");         /* x + x^2 */
    fmpz_poly_compose_series(b, a, c, 5);
    flint_printf("1/(1 - x - x^2) = "); fmpz_poly_print(b); flint_printf("\n");

    fmpz_poly_clear(a);
    fmpz_poly_clear(b);
    fmpz_poly_clear(c);
    return 0;
}
```
```output
1/(1 - x - x^2) = 10  1 1 2 3 5 8 13 21 34 55
reversion       = 10  0 1 1 2 5 14 42 132 429 1430
sqrt exists: 1, 8  1 -2 -2 -4 -10 -28 -84 -264
(3 + 6x)/(3 - 3x) = 6  1 3 3 3 3 3
1/(1 - x - x^2) = 5  1 1 2 3 5
```

Newton iteration and the recurrence give the same inverse; the
example inverts a series with 10-bit coefficients to length 300 both
ways, and checks $Q \cdot Q^{-1} \equiv 1$ at length 600, where
Newton iteration is used. The coefficients of the inverse grow
linearly with the index, here to between 3000 and 4000 bits.

```c example name=fmpz_poly_series_newton
#include <flint/flint.h>
#include <flint/fmpz_poly.h>

int main(void)
{
    flint_rand_t state;
    fmpz_poly_t q, b1, b2, t;
    slong n = 600;

    flint_rand_init(state);
    fmpz_poly_init(q);
    fmpz_poly_init(b1);
    fmpz_poly_init(b2);
    fmpz_poly_init(t);

    fmpz_poly_randtest(q, state, n, 10);
    fmpz_poly_set_coeff_si(q, 0, -1);

    fmpz_poly_inv_series_newton(b1, q, 300);
    fmpz_poly_inv_series_basecase(b2, q, 300);
    flint_printf("newton == basecase: %d\n", fmpz_poly_equal(b1, b2));

    fmpz_poly_inv_series(b1, q, n);
    fmpz_poly_mullow(t, q, b1, n);
    flint_printf("q * (1/q) == 1: %d\n", fmpz_poly_is_one(t));
    flint_printf("bits: %wd\n", FLINT_ABS(fmpz_poly_max_bits(b1)) / 1000 * 1000);

    fmpz_poly_clear(q);
    fmpz_poly_clear(b1);
    fmpz_poly_clear(b2);
    fmpz_poly_clear(t);
    flint_rand_clear(state);
    return 0;
}
```
```output
newton == basecase: 1
q * (1/q) == 1: 1
bits: 3000
```
