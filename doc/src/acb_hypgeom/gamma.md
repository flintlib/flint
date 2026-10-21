---
title: The gamma function of a complex ball
functions: acb_hypgeom_gamma acb_hypgeom_rgamma acb_hypgeom_lgamma
---

`acb_hypgeom_gamma`, `acb_hypgeom_rgamma` and `acb_hypgeom_lgamma`
compute $\Gamma(z)$, $1/\Gamma(z)$ and the principal branch of $\log
\Gamma(z)$ for a complex ball $z$. They are the implementations behind
`acb_gamma`, `acb_rgamma` and `acb_lgamma`. Real arguments are passed
to the real functions (see [](arb_hypgeom/gamma)); otherwise the
functions use a tabulated Taylor series near the real axis and the
Stirling series elsewhere, as described in
[](algorithms/special-functions/gamma).

::: summary algorithms/special-functions/gamma :::

## The principal branch of log-gamma {#lgamma}

The function computed by `acb_hypgeom_lgamma` is not $\log(\Gamma(z))$
with the principal logarithm. It is the analytic continuation of
$\log\Gamma(x)$ from the positive real axis to $\mathbb{C} \setminus (-\infty,
0]$, which is the convention of most computer algebra systems and of
the DLMF [@NIST2012, §5.4]: it satisfies $\log\Gamma(z+1) = \log\Gamma(z)
+ \log z$ exactly (not only modulo $2\pi i$), and its imaginary part
grows like $\operatorname{Im}(z) \log |z|$ instead of staying in
$(-\pi, \pi]$. The two functions differ by $2\pi i k(z)$, where the
integer $k(z)$ changes along curves where $\operatorname{Im}\log\Gamma(z)$
is an odd multiple of $\pi$. On the branch cut, the value is the limit
from the upper half-plane.

Computing this branch requires more than a logarithm of $\Gamma(z)$:

- The Stirling series gives $\log\Gamma(z)$ directly for $z$ of large
  modulus in the right half-plane. Its argument reduction
  $\log\Gamma(z) = \log\Gamma(z+r) - \log\,(z)_r$ needs the
  log-rising factorial $\log\,(z)_r = \sum_{k<r} \log(z+k)$ with the
  correct multiple of $2\pi i$, which `acb_hypgeom_log_rising_ui`
  computes with Hare's algorithm [@Hare1997].
- The reflection formula for the left half-plane becomes $\log\Gamma(z)
  = \log\pi - \log\sin(\pi z) - \log\Gamma(1-z)$, where $\log\sin(\pi
  z)$ must be the analytic branch computed by `acb_log_sin_pi`.
- Near the real axis at moderate precision, the Taylor series of
  $1/\Gamma(1+u)$ is used for $u$ with real part in about $[-1/2,
  1/2]$; the branch of $\log$ applied to it is fixed with a table of
  the lines where $\operatorname{Im}\log\Gamma$ crosses $(k +
  \tfrac12)\pi$.

## Example {#example}

At $z = -10.25 + 0.5i$ the principal branch of $\log\Gamma$ and the
principal logarithm of $\Gamma$ differ by $-10\pi i$:

```c example name=acb_hypgeom_lgamma_branch
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>

int main(void)
{
    acb_t z, y, t;

    acb_init(z);
    acb_init(y);
    acb_init(t);

    acb_set_d_d(z, -10.25, 0.5);

    acb_hypgeom_gamma(y, z, 64);
    flint_printf("gamma(z)      = "); acb_printn(y, 15, 0); flint_printf("\n");

    acb_hypgeom_lgamma(y, z, 64);
    flint_printf("lgamma(z)     = "); acb_printn(y, 15, 0); flint_printf("\n");

    acb_hypgeom_gamma(t, z, 64);
    acb_log(t, t, 64);
    flint_printf("log(gamma(z)) = "); acb_printn(t, 15, 0); flint_printf("\n");

    /* the difference is 2 pi i k */
    acb_sub(t, y, t, 64);
    acb_const_pi(y, 64);
    acb_mul_2exp_si(y, y, 1);
    acb_div(t, t, y, 64);
    flint_printf("difference / (2 pi) = "); acb_printn(t, 15, 0); flint_printf("\n");

    acb_hypgeom_rgamma(y, z, 64);
    flint_printf("rgamma(z)     = "); acb_printn(y, 15, 0); flint_printf("\n");

    acb_clear(z);
    acb_clear(y);
    acb_clear(t);
    flint_cleanup();
    return 0;
}
```
```output
gamma(z)      = [7.08331784482241e-8 +/- 3.66e-23] + [-1.88626361829666e-7 +/- 2.32e-22]*I
lgamma(z)     = [-15.4175381117796 +/- 4.86e-15] + [-32.6274954179361 +/- 2.24e-14]*I
log(gamma(z)) = [-15.4175381117796 +/- 4.86e-15] + [-1.21156888203815 +/- 4.69e-15]*I
difference / (2 pi) = [+/- 2.62e-19] + [-5.00000000000000 +/- 1.09e-18]*I
rgamma(z)     = [1744778.02581992 +/- 4.92e-9] + [4646284.95319217 +/- 4.09e-9]*I
```

$\Gamma(z)$ is small here because $z$ is close to the negative real
axis between two poles; $1/\Gamma(z)$, which is entire, is large.
