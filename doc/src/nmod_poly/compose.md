---
title: Composition and Taylor shifts
functions: nmod_poly_compose nmod_poly_compose_horner nmod_poly_taylor_shift nmod_poly_taylor_shift_horner nmod_poly_taylor_shift_convolution
---

The composition $f(g)$ of polynomials of degrees $d$ and $e$ has degree
$de$, so even writing it down costs $de$ operations. `nmod_poly_compose`
computes it in $O(M(de) \log d)$ operations, close to this bound. The
most common special case, $g = x + c$, is a Taylor shift, which changes
the basis from powers of $x$ to powers of $x - c$ and costs only one
multiplication of length $d$. For compositions modulo a third polynomial
see [](nmod_poly/compose_mod), and for power series composition
`nmod_poly_compose_series`.

::: summary algorithms/polynomials/composition :::

## Algorithms {#algorithms}

`nmod_poly_compose` handles constant $g$ by an evaluation, and $g = a x^k
+ c$ by a Taylor shift followed by scaling the coefficients by powers of
$a$ and spreading them out. Otherwise it uses Horner's rule, $f(g) =
(\cdots(f_d g + f_{d-1}) g + \cdots) g + f_0$, for $d \le 6$, and the
divide-and-conquer algorithm for larger $d$: with $f = f_{\mathrm{lo}} +
x^{k} f_{\mathrm{hi}}$ for $k$ a power of two about $d/2$,
$$
f(g) = f_{\mathrm{lo}}(g) + g^{k} f_{\mathrm{hi}}(g),
$$
where the powers $g^{2^j}$ are computed once. The last level multiplies
polynomials of length about $de/2$, which dominates the cost.

The Taylor shift $f(x + c)$ by Horner's rule takes $d^2/2$ multiply-adds
(only additions for $c = \pm 1$) and works for every modulus. The
convolution method [@ASU1975] [@GG1997] uses
$$
f(x + c) = \sum_j \frac{x^j}{j!} \sum_{k \ge j} \bigl(k!\, f_k\bigr)
\frac{c^{k-j}}{(k-j)!},
$$
in which the inner sums are the coefficients of one product of length $d
+ 1$, after reversing $(k! f_k)$. It divides by $d!$, so it needs a prime
modulus larger than $d$ (or one without prime factors up to $d$); for
other moduli it throws `FLINT_IMPINV`. `nmod_poly_taylor_shift` uses
Horner's rule for lengths below 100, for lengths larger than the modulus,
and for $c = \pm 1$ below length 1000, and the convolution method
otherwise. A composite modulus with a small prime factor and a length
between 100 and the modulus therefore needs
`nmod_poly_taylor_shift_horner`.

## Examples {#examples}

Composition agrees with evaluation, $f(g)(c) = f(g(c))$, and the Taylor
shift by $c$ moves the roots by $-c$:

```c example name=nmod_poly_compose_basic
#include <flint/flint.h>
#include <flint/nmod_poly.h>

int main(void)
{
    flint_rand_t state;
    nmod_poly_t f, g, h;
    ulong p = 1000003, c = 12345, roots[3] = { 1, 2, 3 };

    flint_rand_init(state);
    nmod_poly_init(f, p);
    nmod_poly_init(g, p);
    nmod_poly_init(h, p);

    nmod_poly_rand(f, state, 200);
    nmod_poly_rand(g, state, 30);
    nmod_poly_compose(h, f, g);
    flint_printf("deg f(g) = %wd\n", nmod_poly_degree(h));
    flint_printf("f(g)(c) == f(g(c)): %d\n",
        nmod_poly_evaluate_nmod(h, c) ==
        nmod_poly_evaluate_nmod(f, nmod_poly_evaluate_nmod(g, c)));

    /* (x - 1)(x - 2)(x - 3) shifted by 10 has roots -9, -8, -7 */
    nmod_poly_product_roots_nmod_vec(f, roots, 3);
    nmod_poly_taylor_shift(g, f, 10);
    flint_printf("f(x + 10) = "); nmod_poly_print_pretty(g, "x");
    flint_printf("\ng(-9) = %wu\n", nmod_poly_evaluate_nmod(g, p - 9));

    /* the two Taylor shift algorithms agree */
    nmod_poly_rand(f, state, 5000);
    nmod_poly_taylor_shift_horner(g, f, c);
    nmod_poly_taylor_shift_convolution(h, f, c);
    flint_printf("horner == convolution: %d\n", nmod_poly_equal(g, h));

    nmod_poly_clear(f);
    nmod_poly_clear(g);
    nmod_poly_clear(h);
    flint_rand_clear(state);
    return 0;
}
```
```output
deg f(g) = 5771
f(g)(c) == f(g(c)): 1
f(x + 10) = x^3+24*x^2+191*x+504
g(-9) = 0
horner == convolution: 1
```
