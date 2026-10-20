---
title: Using acb_hypgeom
summary: How to call the special functions of acb_hypgeom and arb_hypgeom, what their output balls mean, how precision, inexact input, regularized functions, optional outputs and power series arguments work, and what to do when a result is too wide.
order: 1
---

The functions of [`acb_hypgeom`](acb_hypgeom) and
[`arb_hypgeom`](arb_hypgeom) follow the conventions of ball arithmetic
(see [](guide/types/ball-arithmetic)): the output is a ball that is
guaranteed to contain the exact function value, for every point of the
input balls. The algorithms are described in
[](algorithms/special-functions/overview).

## A first program {#first}

```c example name=acb_hypgeom_using
#include <flint/acb.h>
#include <flint/acb_poly.h>
#include <flint/acb_hypgeom.h>

int main(void)
{
    acb_t nu, z, r, a, b;
    acb_poly_t x, y;

    acb_init(nu);
    acb_init(z);
    acb_init(r);
    acb_init(a);
    acb_init(b);
    acb_poly_init(x);
    acb_poly_init(y);

    /* J_0(10) and J_0(10^6) */
    acb_zero(nu);
    acb_set_ui(z, 10);
    acb_hypgeom_bessel_j(r, nu, z, 64);
    flint_printf("J0(10)        = "); acb_printn(r, 18, 0); flint_printf("\n");
    acb_set_ui(z, 1000000);
    acb_hypgeom_bessel_j(r, nu, z, 64);
    flint_printf("J0(10^6)      = "); acb_printn(r, 18, 0); flint_printf("\n");

    /* an inexact argument: the radius propagates */
    acb_set_ui(z, 10);
    mag_set_d(arb_radref(acb_realref(z)), 1e-6);
    acb_hypgeom_bessel_j(r, nu, z, 64);
    flint_printf("J0(10+-1e-6)  = "); acb_printn(r, 18, 0); flint_printf("\n");

    /* regularized 1F1 at a pole of the ordinary function */
    acb_set_d(a, 0.5);
    acb_set_si(b, -1);
    acb_set_d(z, 2.0);
    acb_hypgeom_m(r, a, b, z, 0, 64);
    flint_printf("1F1(1/2,-1,2) = "); acb_printn(r, 18, 0); flint_printf("\n");
    acb_hypgeom_m(r, a, b, z, 1, 64);
    flint_printf("M(1/2,-1,2)   = "); acb_printn(r, 18, 0); flint_printf("\n");

    /* only one of the four Airy functions */
    acb_set_si(z, -5);
    acb_hypgeom_airy(r, NULL, NULL, NULL, z, 64);
    flint_printf("Ai(-5)        = "); acb_printn(r, 18, 0); flint_printf("\n");

    /* Taylor coefficients of erf at 1: erf(1 + t) = c0 + c1 t + c2 t^2 + ... */
    acb_poly_set_coeff_si(x, 0, 1);
    acb_poly_set_coeff_si(x, 1, 1);
    acb_hypgeom_erf_series(y, x, 3, 64);
    for (slong k = 0; k < 3; k++)
    {
        flint_printf("[t^%wd] erf(1+t) = ", k);
        acb_printn(y->coeffs + k, 18, 0);
        flint_printf("\n");
    }

    acb_clear(nu);
    acb_clear(z);
    acb_clear(r);
    acb_clear(a);
    acb_clear(b);
    acb_poly_clear(x);
    acb_poly_clear(y);
    flint_cleanup();
    return 0;
}
```
```output
J0(10)        = [-0.245935764451348335 +/- 2.11e-19]
J0(10^6)      = [0.0003310430137399 +/- 8.69e-17]
J0(10+-1e-6)  = [-0.245936 +/- 4.97e-7]
1F1(1/2,-1,2) = nan
M(1/2,-1,2)   = [8.41930946082163841 +/- 9.04e-18]
Ai(-5)        = [0.350761009024114320 +/- 2.62e-19]
[t^0] erf(1+t) = [0.842700792949714869 +/- 3.69e-19]
[t^1] erf(1+t) = [0.415107497420594703 +/- 4.59e-19]
[t^2] erf(1+t) = [-0.415107497420594703 +/- 4.59e-19]
```

## Precision and accuracy {#precision}

The precision argument (64 bits above) is the precision of the
arithmetic, not a promise about the result. For moderate arguments the
functions estimate the cancellation of the algorithm they use and add
guard bits, so that the relative error of the result is close to
$2^{-\mathit{prec}}$: $J_0(10)$ above has about 60 correct bits. This
fails near the zeros of a function, where the relative error is
unbounded, and the absolute error is then what matters; $J_0(10^6)
\approx 3.3 \cdot 10^{-4}$ is computed with the asymptotic expansion
and an absolute error of about $2^{-53}$.

To get a result with a guaranteed accuracy of $p$ bits, evaluate at a
precision of $p$ plus a few guard bits, check the accuracy with
`acb_rel_accuracy_bits`, and repeat at a higher precision if it is too
low. Precision doubling is the usual strategy; the cost of most
functions is dominated by the last evaluation.

## Inexact input {#inexact}

If the input is a ball of positive radius, the output contains the
image of the whole ball. Most functions evaluate at the midpoint and
add a bound for the propagated error, obtained from a bound for the
derivative on the ball, as in the third line above (the bound for
$|J_0'|$ used there is crude, so the radius is larger than $|J_1(10)|
\cdot 10^{-6}$). The functions also reduce their working precision to
what the input accuracy justifies, so a wide input is cheap. A very
wide input gives a correct but possibly very pessimistic enclosure; to
evaluate a function on a large interval, subdivide it.

When a parameter is supposed to be an exact integer, make sure that it
is an exact ball: many formulas have removable singularities at
integer parameters (for example $U(a,b,z)$, $Y_\nu(z)$ and $K_\nu(z)$
at integer $b$ or $\nu$, and the transformations of ${}_2F_1$ when a
parameter difference is an integer), which the functions handle by
computing a limit only if they can recognize the integer. An inexact
ball containing an integer gives an indeterminate result. For
${}_2F_1$, flags can declare a parameter difference as an exact
integer even if the parameters are inexact (see [](acb_hypgeom/2f1)).

## Indeterminate results and regularization {#nan}

An indeterminate result (printed as `nan`) means that no finite
enclosure was found: the input contains a pole or straddles a branch
cut, the precision is too low for the parameters, or an algorithm is
not implemented for this case. The ordinary function
${}_1F_1(1/2, -1, 2)$ in the program is undefined, because the series
divides by $(b)_k$ which vanishes; the regularized function
$\mathbf{M}(a,b,z) = {}_1F_1(a,b,z)/\Gamma(b)$ is entire in $b$, and
`acb_hypgeom_m` with `regularized` = 1 computes it. Use the
regularized versions whenever a lower parameter can be a nonpositive
integer.

The real functions of `arb_hypgeom` are indeterminate also where the
function value is not real: `arb_hypgeom_bessel_y` at a negative
argument, `arb_hypgeom_lgamma` at a negative argument,
`arb_hypgeom_dilog` for $z > 1$. Use the complex functions there.

## Several outputs {#outputs}

Functions that compute several related functions at once, such as the
Airy functions, the Fresnel integrals, the Coulomb wave functions and
$J_\nu$ together with $Y_\nu$, accept `NULL` for unwanted outputs and
then skip the work for them, as for $\operatorname{Ai}(-5)$ above.

## Power series and derivatives {#series}

The functions with the suffix `_series` evaluate a function at a power
series argument, truncated to a given length. With the argument $x_0 +
t$, the coefficients of the result are the Taylor coefficients
$f^{(k)}(x_0)/k!$; above, $\operatorname{erf}'(1) = 2e^{-1}/\sqrt\pi
\approx 0.4151$ and $\operatorname{erf}''(1)/2 = -2e^{-1}/\sqrt\pi$.
The functions with the suffix `_jet` take a point and return the
Taylor coefficients directly. Power series arguments also propagate
derivatives through compositions.
