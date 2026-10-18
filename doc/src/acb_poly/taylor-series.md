---
title: Taylor series of complex functions
summary: Computing derivatives of complex special functions with acb_poly power series, illustrated with the Riemann zeta function, the Stieltjes constants and the derivative of a polylogarithm with respect to its order.
order: 1
---

[](acb_poly) supports the same power series idiom as [](arb_poly)
(see [](arb_poly/power-series)): the series $s_0 + x$ of length 2,
passed through the `_series` functions with a target length $n$, gives
the first $n$ Taylor coefficients of the composite function at $s_0$,
with rigorous error bounds. Besides the elementary functions, the
module provides series versions of the gamma function, the Hurwitz
zeta function, polylogarithms with a series order and the Lambert W
function, and [](acb_hypgeom), [](acb_elliptic) and [](acb_dirichlet)
add many more (`acb_hypgeom_erf_series`, `acb_dirichlet_hardy_z_series`,
`acb_elliptic_k_series`, ...).

## The zeta function {#zeta}

`acb_poly_zeta_series(z, s, a, deflate, n, prec)` expands the Hurwitz
zeta function $\zeta(s, a)$ in the series $s$. With `deflate` set, it
expands the entire function $\zeta(s, a) - 1/(s - 1)$ instead, which is
needed at $s = 1$: there, $\zeta(1 + x) - 1/x = \sum_k (-1)^k \gamma_k
x^k / k!$ gives the Stieltjes constants $\gamma_k$. The program also
evaluates $\zeta$ and $\zeta'$ at an approximation of the first
nontrivial zero $\rho_1 = 1/2 + 14.1347\ldots i$, and the derivative
with respect to $s$ of the polylogarithm $\operatorname{Li}_s(1/2)$ at
$s = 2$ with `acb_poly_polylog_series`:

```c example name=acb_poly_taylor_zeta
#include <flint/acb_poly.h>

int main(void)
{
    acb_poly_t s, z;
    acb_t a, t;
    slong k, prec = 64;

    acb_poly_init(s);
    acb_poly_init(z);
    acb_init(a);
    acb_init(t);

    /* Stieltjes constants from zeta(1 + x) - 1/x */
    acb_poly_set_coeff_si(s, 0, 1);
    acb_poly_set_coeff_si(s, 1, 1);
    acb_one(a);
    acb_poly_zeta_series(z, s, a, 1, 4, prec);
    for (k = 0; k < 4; k++)
    {
        /* gamma_k = (-1)^k k! [x^k] */
        acb_zero(t);
        arb_fac_ui(acb_realref(t), k, prec);
        acb_mul(t, t, z->coeffs + k, prec);
        if (k % 2 == 1)
            acb_neg(t, t);
        flint_printf("gamma_%wd = ", k);
        acb_printn(t, 15, 0);
        flint_printf("\n");
    }

    /* zeta(s) and zeta'(s) at s = 1/2 + 14.134725141734693790 i */
    acb_poly_zero(s);
    acb_set_d_d(t, 0.5, 0);
    arb_set_str(acb_imagref(t), "14.134725141734693790", prec);
    acb_poly_set_coeff_acb(s, 0, t);
    acb_poly_set_coeff_si(s, 1, 1);
    acb_poly_zeta_series(z, s, a, 0, 2, prec);
    flint_printf("zeta  = ");
    acb_printn(z->coeffs + 0, 10, 0);
    flint_printf("\nzeta' = ");
    acb_printn(z->coeffs + 1, 10, 0);
    flint_printf("\n");

    /* Li_s(1/2) and its derivative with respect to s at s = 2 */
    acb_poly_zero(s);
    acb_poly_set_coeff_si(s, 0, 2);
    acb_poly_set_coeff_si(s, 1, 1);
    acb_set_d(t, 0.5);
    acb_poly_polylog_series(z, s, t, 2, prec);
    flint_printf("Li_2(1/2) = ");
    acb_printn(z->coeffs + 0, 15, 0);
    flint_printf("\nd/ds Li_s(1/2) = ");
    acb_printn(z->coeffs + 1, 15, 0);
    flint_printf("\n");

    acb_poly_clear(s);
    acb_poly_clear(z);
    acb_clear(a);
    acb_clear(t);
    return 0;
}
```
```output
gamma_0 = [0.577215664901533 +/- 1.40e-16]
gamma_1 = [-0.0728158454836767 +/- 2.49e-17]
gamma_2 = [-0.00969036319287232 +/- 1.52e-18]
gamma_3 = [0.00205383442030335 +/- 4.14e-18]
zeta  = [+/- 1.44e-18] + [+/- 2.08e-18]*I
zeta' = [0.7832965119 +/- 3.30e-11] + [0.1246998297 +/- 4.82e-11]*I
Li_2(1/2) = [0.582240526465013 +/- 5.57e-16]
d/ds Li_s(1/2) = [-0.067314759696392 +/- 3.56e-16] + [+/- 2.58e-16]*I
```

The value of $\zeta$ at the approximate zero is a ball that contains
zero, as it must, since the input approximates $\rho_1$ to about 19
digits and $|\zeta'(\rho_1)| \approx 0.79$. The derivative of
$\operatorname{Li}_s(1/2)$ is real, but the computation goes through
complex arithmetic (a combination of Hurwitz zeta functions), which
leaves an imaginary part that is only known to be small; for a real
result, take the real part.

## Efficiency {#efficiency}

The cost of the series of $\zeta(s, a)$ of length $n$ at precision $p$
is dominated by the Euler–Maclaurin sum $\sum_{k<N} (k + a)^{-s-x}$
with $N$ roughly proportional to $p$ (and to $|\operatorname{Im}(s)|$
for large imaginary parts). The power sum is computed term by term in
$O(Nn)$ operations, or by the transposed multipoint evaluation of
`_acb_poly_powsum_series_tree` in $\tilde O(N + n)$ operations when $n$
and $N$ are both large, which makes thousands of derivatives feasible
([](algorithms/special-functions/zeta)). For $a = 1$ the series at the
constant term comes from `acb_dirichlet_zeta_jet`, which uses the
functional equation in the left half-plane and, for the value and first
derivative at very large heights, the Riemann–Siegel formula.
