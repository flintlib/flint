---
title: Evaluating Jacobi theta functions
functions: acb_modular_theta acb_modular_theta_jet acb_modular_theta_series
---

`acb_modular_theta` computes all four Jacobi theta functions
$\theta_1(z, \tau), \ldots, \theta_4(z, \tau)$ at once, since they share
the powers of $q$ and $w$; `acb_modular_theta_jet` also computes
derivatives with respect to $z$, and `acb_modular_theta_series` composes
with a power series. They are the basis of the Weierstrass functions of
[`acb_elliptic`](acb_elliptic), and the $g = 1$ case of the Riemann theta
functions of [`acb_theta`](acb_theta) is related to them by
$(\theta_{0,0}, \theta_{0,1}, \theta_{1,0}, \theta_{1,1}) = (\theta_3, \theta_4, \theta_2, -\theta_1)$.

## Evaluation {#evaluation}

The evaluation has three steps [@Joh2019]; the mathematics is in
[](algorithms/special-functions/theta-functions#genus-one).

1. **Reduce $\tau$.** `acb_modular_fundamental_domain_approx` finds
   $g \in \mathrm{PSL}_2(\Z)$ such that $\tau' = g\tau$ satisfies
   $|\operatorname{Re}(\tau')| \le 1/2$ and $|\tau'|^2 \ge 63/64$, so that
   $|q'| = e^{-\pi \operatorname{Im}(\tau')} < 0.07$.
   `acb_modular_theta_transform` gives the permutation of the four
   functions, the 8th roots of unity, and whether the factors
   $A = \sqrt{i/(c\tau + d)}$, $B = e^{-\pi i c z^2/(c\tau + d)}$ and the
   new variable $z' = -z/(c\tau + d)$ are needed.
2. **Reduce $z$.** If $|\operatorname{Im}(z')| > 16$, an integer $n$
   close to $\operatorname{Im}(z')/\operatorname{Im}(\tau')$ is
   subtracted: $z' \gets z' - n\tau'$, using the quasi-periodicity
   $\theta_3(z + n\tau, \tau) = e^{-\pi i (n^2 \tau + 2nz)} \theta_3(z, \tau)$
   and its analogues (which add a sign $(-1)^n$ to $\theta_1$ and
   $\theta_4$). This keeps $|w^{\pm 1}|$, $w = e^{\pi i z'}$, from
   overwhelming the decay of the powers of $q'$.
3. **Sum.** `acb_modular_theta_sum` sums the four series with
   $N = O(\sqrt{p})$ terms at precision $p$, computing the powers of $q'$
   with an addition sequence and the needed even powers of $w^{\pm 1}$,
   and adds a rigorous bound for the tail. The result is multiplied by
   the factors from step 1.

The cost is thus $O(\sqrt{p})$ multiplications at precision $p$ plus a
few exponentials, independently of how close $\tau$ is to the real
axis; but the magnitude of the result can be huge or tiny when it is
(as in the example below), and precision is lost in proportion to
$\log |B|$ when $z$ is far from the real line. The `notransform` versions
skip steps 1 and 2, which is only advisable when
$\operatorname{Im}(\tau)$ is already large.

For `acb_modular_theta_jet` the factor $B$ is a power series in $x$
(for $z \to z + x$), computed with `_acb_poly_exp_pi_i_series`, and the
coefficients of the transformed functions are scaled by powers of
$-1/(c\tau + d)$. If $\tau$ is already reduced and $|\operatorname{Im}(z)| \le 16$,
the `notransform` version is called directly.

When $\operatorname{Re}(\tau)$ is an exact integer and $z$ is exactly
real or exactly imaginary, `acb_modular_theta` gives the outputs that are
known to be real or imaginary an exact zero imaginary or real part.

## Derivatives {#derivatives}

The $k$-th entry of the output of `acb_modular_theta_jet` is
$\theta_i^{(k)}(z, \tau) / k!$, where the derivative is taken with respect
to the variable $z$ of this module (for which the period is 1). For
example $\theta_1'(0, \tau) = \pi \theta_2(0, \tau) \theta_3(0, \tau) \theta_4(0, \tau)$,
with the factor $\pi$ that does not appear when the variable is $\pi z$,
as in [@NIST2012, 20.4.6]. Derivatives with respect to $\tau$ follow from
the heat equation $4\pi i \, \partial_\tau \theta_j = \partial_z^2 \theta_j$.

## Example {#example}

The program evaluates the theta functions at
$\tau = 1/2 + 0.001 i$, close to the real line, where
$\theta_3(0.25, \tau)$ is of order $10^{-84}$, then checks Jacobi's
identity $\theta_3^4 = \theta_2^4 + \theta_4^4$ at $z = 0$ and the
derivative formula above.

```c example name=acb_modular_theta
#include <flint/acb.h>
#include <flint/acb_modular.h>

int main(void)
{
    acb_t z, tau, t1, t2, t3, t4, u;
    acb_ptr j1, j2, j3, j4;
    slong prec = 128;

    acb_init(z); acb_init(tau); acb_init(u);
    acb_init(t1); acb_init(t2); acb_init(t3); acb_init(t4);
    j1 = _acb_vec_init(2); j2 = _acb_vec_init(2);
    j3 = _acb_vec_init(2); j4 = _acb_vec_init(2);

    acb_set_d_d(tau, 0.5, 0.001);
    acb_set_d_d(z, 0.25, 0.0);
    acb_modular_theta(t1, t2, t3, t4, z, tau, prec);
    flint_printf("theta_1 = "); acb_printn(t1, 20, 0); flint_printf("\n");
    flint_printf("theta_3 = "); acb_printn(t3, 20, 0); flint_printf("\n");

    acb_zero(z);
    acb_modular_theta(t1, t2, t3, t4, z, tau, prec);
    acb_pow_ui(t2, t2, 4, prec);
    acb_pow_ui(t3, t3, 4, prec);
    acb_pow_ui(t4, t4, 4, prec);
    acb_add(u, t2, t4, prec);
    acb_sub(u, u, t3, prec);
    flint_printf("theta_2^4 + theta_4^4 - theta_3^4 = ");
    acb_printn(u, 5, 0); flint_printf("\n");

    acb_set_d_d(tau, 0.1, 0.8);
    acb_modular_theta_jet(j1, j2, j3, j4, z, tau, 2, prec);
    acb_mul(u, j2, j3, prec);
    acb_mul(u, u, j4, prec);
    acb_const_pi(t1, prec);
    acb_mul(u, u, t1, prec);
    flint_printf("theta_1'(0)                = ");
    acb_printn(j1 + 1, 20, 0); flint_printf("\n");
    flint_printf("pi theta_2 theta_3 theta_4 = ");
    acb_printn(u, 20, 0); flint_printf("\n");

    acb_clear(z); acb_clear(tau); acb_clear(u);
    acb_clear(t1); acb_clear(t2); acb_clear(t3); acb_clear(t4);
    _acb_vec_clear(j1, 2); _acb_vec_clear(j2, 2);
    _acb_vec_clear(j3, 2); _acb_vec_clear(j4, 2);
    return 0;
}
```
```output
theta_1 = [20.658574377159641560 +/- 1.28e-19] + [8.5570616863128383887 +/- 1.67e-20]*I
theta_3 = [1.6845266452348015023e-84 +/- 1.73e-104] + [+/- 2.39e-119]*I
theta_2^4 + theta_4^4 - theta_3^4 = [+/- 3.11e-32] + [+/- 2.57e-32]*I
theta_1'(0)                = [3.2914967126912256704 +/- 1.45e-20] + [0.22014791420183894914 +/- 2.43e-21]*I
pi theta_2 theta_3 theta_4 = [3.2914967126912256704 +/- 1.45e-20] + [0.22014791420183894914 +/- 2.43e-21]*I
```
