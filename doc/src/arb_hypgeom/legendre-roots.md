---
title: Roots of Legendre polynomials and Gauss–Legendre nodes
functions: arb_hypgeom_legendre_p_ui_root _arb_hypgeom_legendre_p_ui_root_bitburst _arb_hypgeom_legendre_p_ui_root_taylor
---

`arb_hypgeom_legendre_p_ui_root` computes an enclosure of the $k$-th
root $x_k$ of the Legendre polynomial $P_n$ and, optionally, the
weight $w_k = 2/((1-x_k^2) P_n'(x_k)^2)$, so that
$$
\int_{-1}^{1} f(x)\,dx \approx \sum_{k=0}^{n-1} w_k f(x_k)
$$
is the $n$-point Gauss–Legendre quadrature rule, exact for polynomials
of degree less than $2n$. The roots are numbered in decreasing order,
$1 > x_0 > x_1 > \cdots > x_{n-1} > -1$, which is the order of
increasing $\theta$ for the roots of $P_n(\cos\theta)$. Only the roots
with $k < \lceil n/2 \rceil$ need to be computed, since $x_{n-1-k} =
-x_k$ and $w_{n-1-k} = w_k$; the function uses this symmetry itself.
An index $k \ge n$ raises an exception. The function is the one used
by `acb_calc_integrate` to generate its quadrature rules
([](algorithms/ball-arithmetic/integration)). The algorithms are
described in [](algorithms/special-functions/legendre#roots).

::: summary algorithms/special-functions/legendre :::

## Example {#example}

The 20-point rule integrates $e^x$ over $[-1, 1]$ to more than 30
digits; the error of the rule itself is about $10^{-60}$:

```c example name=arb_hypgeom_gauss_legendre
#include <flint/arb.h>
#include <flint/arb_hypgeom.h>

int main(void)
{
    arb_t x, w, s, t, wsum;
    ulong n = 20, k;
    slong prec = 128;

    arb_init(x);
    arb_init(w);
    arb_init(s);
    arb_init(t);
    arb_init(wsum);

    for (k = 0; k < n; k++)
    {
        arb_hypgeom_legendre_p_ui_root(x, w, n, k, prec);
        if (k == 0)
        {
            flint_printf("x_0 = "); arb_printn(x, 30, 0); flint_printf("\n");
            flint_printf("w_0 = "); arb_printn(w, 30, 0); flint_printf("\n");
        }
        /* s += w_k exp(x_k) */
        arb_exp(t, x, prec);
        arb_addmul(s, w, t, prec);
        arb_add(wsum, wsum, w, prec);
    }

    flint_printf("sum of weights = "); arb_printn(wsum, 30, 0); flint_printf("\n");
    flint_printf("quadrature     = "); arb_printn(s, 30, 0); flint_printf("\n");

    /* exact value e - 1/e */
    arb_const_e(t, prec);
    arb_inv(x, t, prec);
    arb_sub(t, t, x, prec);
    flint_printf("e - 1/e        = "); arb_printn(t, 30, 0); flint_printf("\n");

    arb_clear(x);
    arb_clear(w);
    arb_clear(s);
    arb_clear(t);
    arb_clear(wsum);
    flint_cleanup();
    return 0;
}
```
```output
x_0 = [0.993128599185094924786122388471 +/- 3.21e-31]
w_0 = [0.0176140071391521183118619623519 +/- 4.72e-32]
sum of weights = [2.00000000000000000000000000000 +/- 1.08e-37]
quadrature     = [2.35040238728760291376476370119 +/- 1.21e-30]
e - 1/e        = [2.35040238728760291376476370119 +/- 1.21e-30]
```

The enclosure of the quadrature sum is rigorous for the sum, but it
does not include the error of the quadrature rule. For an integral
with a proved error bound, use `acb_calc_integrate`, which bounds the
quadrature error using the analyticity of the integrand
[@Joh2018a].
