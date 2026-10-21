---
title: Composition and Taylor shifts
functions: fmpz_poly_compose fmpz_poly_compose_horner fmpz_poly_taylor_shift fmpz_poly_taylor_shift_horner fmpz_poly_taylor_shift_divconquer fmpz_poly_taylor_shift_multi_mod
---

`fmpz_poly_compose` computes $g(h(x))$, a polynomial of degree $\deg g
\cdot \deg h$. The *Taylor shift* $f(x + c)$ is the special case $h = x
+ c$, which occurs often enough (for example in root isolation) to have
its own functions.
For composition of power series, see `fmpz_poly_compose_series`.

::: summary algorithms/polynomials/composition :::

## Algorithms {#algorithms}

Horner's rule $g(h) = (\cdots(g_m h +
g_{m-1}) h + \cdots) h + g_0$ performs $m = \deg g$ products of
increasingly unbalanced size. The divide-and-conquer algorithm groups
the coefficients of $g$ in pairs, $g_{2i} + g_{2i+1} h$, then combines
pairs with $h^2$, quadruples with $h^4$ and so on, so that all products
are balanced; with fast multiplication its cost is $O(M(mn) \log m)$
for $n = \deg h$, against $O(m M(mn))$ for Horner's rule.
`fmpz_poly_compose` uses Horner's rule for $\operatorname{len}(g) \le
4$ and divide and conquer otherwise, and reduces inner polynomials of
the form $a x^k + b$ to a Taylor shift, a scaling and an inflation.

For the Taylor shift, Horner's scheme $f \mapsto f(x + c)$ costs
$O(n^2)$ additions and multiplications by $c$ (only additions for $c =
\pm 1$), on coefficients that grow to $O(n \log |c| + n + b)$ bits.
The divide-and-conquer algorithm writes $f = f_0 + x^m f_1$, shifts the
halves recursively and computes $(x + c)^m$ by the binomial theorem, so
that the work is in a few large products [@GG1997]. `fmpz_poly_taylor_shift`
uses Horner's scheme below length 64 and divide and conquer otherwise;
the divide-and-conquer function itself returns to Horner's scheme
below a length between 100 and 1000 that grows with the coefficient
size, and runs the two recursive halves in parallel when threads are
available. `fmpz_poly_taylor_shift_multi_mod` computes the shift
modulo many word-size primes in parallel and reconstructs the result,
which is competitive for long polynomials with small coefficients
and many threads.

## Examples {#examples}

Composition with $h = x^2 + 1$ and a Taylor shift by $-1$; the shift
of a polynomial by $c$ and then by $-c$ returns the original:

```c example name=fmpz_poly_compose_shift
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

int main(void)
{
    fmpz_poly_t g, h, r, s;
    fmpz_t c;

    fmpz_poly_init(g);
    fmpz_poly_init(h);
    fmpz_poly_init(r);
    fmpz_poly_init(s);
    fmpz_init(c);

    fmpz_poly_set_str(g, "3  1 -3 2");        /* 2x^2 - 3x + 1 */
    fmpz_poly_set_str(h, "3  1 0 1");         /* x^2 + 1 */
    fmpz_poly_compose(r, g, h);
    flint_printf("g(h(x)) = "); fmpz_poly_print_pretty(r, "x"); flint_printf("\n");

    fmpz_set_si(c, -1);
    fmpz_poly_taylor_shift(r, g, c);
    flint_printf("g(x - 1) = "); fmpz_poly_print_pretty(r, "x"); flint_printf("\n");

    /* the shift agrees with composition by x + c */
    fmpz_poly_chebyshev_t(g, 100);
    fmpz_set_si(c, 7);
    fmpz_poly_taylor_shift(r, g, c);
    fmpz_poly_set_str(h, "2  7 1");
    fmpz_poly_compose(s, g, h);
    flint_printf("shift == compose: %d\n", fmpz_poly_equal(r, s));
    fmpz_neg(c, c);
    fmpz_poly_taylor_shift(r, r, c);
    flint_printf("shift back: %d\n", fmpz_poly_equal(r, g));

    fmpz_poly_clear(g);
    fmpz_poly_clear(h);
    fmpz_poly_clear(r);
    fmpz_poly_clear(s);
    fmpz_clear(c);
    return 0;
}
```
```output
g(h(x)) = 2*x^4+x^2
g(x - 1) = 2*x^2-7*x+6
shift == compose: 1
shift back: 1
```
