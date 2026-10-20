---
title: acb_elliptic.h
header: acb_elliptic.h
description: Complete and incomplete elliptic integrals, Carlson symmetric integrals and Weierstrass elliptic functions of complex variables.
category: Special functions
order: 3
---

This module evaluates elliptic integrals and elliptic (doubly periodic)
functions of complex variables with rigorous error bounds:

- the complete elliptic integrals $K(m)$, $E(m)$ and $\Pi(n, m)$,
  computed with the arithmetic-geometric mean;
- the incomplete Legendre integrals $F(\phi, m)$, $E(\phi, m)$ and
  $\Pi(n, \phi, m)$, reduced to Carlson's symmetric integrals;
- the Carlson symmetric integrals $R_F$, $R_G$, $R_J$ and the special
  cases $R_C$ and $R_D$, computed by the duplication algorithm with a
  high-order Taylor expansion;
- the Weierstrass functions $\wp$, $\wp'$, $\wp^{-1}$, $\zeta$ and
  $\sigma$ for the lattice $\Z + \tau\Z$, and the lattice invariants and
  roots, computed with Jacobi theta functions.

The modular forms and theta functions on which the Weierstrass functions
are built are in [`acb_modular`](acb_modular) (genus 1) and
[`acb_theta`](acb_theta) (any genus). The algorithms are described in
[](algorithms/special-functions/elliptic-modular), following
[@Joh2019].

::: summary algorithms/special-functions/elliptic-modular :::

## Conventions {#conventions}

- The parameter $m$ of the Legendre integrals is the *parameter*, related
  to the modulus $k$ and the modular angle $\alpha$ by
  $m = k^2 = \sin^2 \alpha$ (as in [@NIST2012, §19.1]).
- In the Weierstrass functions, the lattice is
  $\Lambda = \Z + \tau \Z$ with $\operatorname{Im}(\tau) > 0$, that is,
  the half-periods are $\omega_1 = 1/2$ and $\omega_3 = \tau/2$. A
  general lattice $2\omega_1 \Z + 2\omega_3 \Z$ is handled by the
  homogeneity relations, for instance
  $\wp(z; 2\omega_1, 2\omega_3) = (2\omega_1)^{-2} \wp(z/(2\omega_1), \omega_3/\omega_1)$.
- A nonzero flag `times_pi` replaces the amplitude $\phi$ by $\pi \phi$,
  which makes the quasi-period 1 instead of $\pi$ and allows exact input
  at multiples of $\pi/2$.
- The `flags` arguments of the Carlson integrals are reserved for future
  use; pass 0.
- Outputs may alias inputs. A function returns indeterminate balls (or
  balls containing the result with infinite radius) when the input is not
  finite, lies on or straddles a singularity, or the algorithm cannot
  separate the input from a branch cut.

::: warning
Incomplete elliptic integrals have a complicated branch structure when
extended to complex variables. For some functions in this module, branch
cuts are artifacts of the evaluation algorithm rather than having a
natural mathematical justification, so the values may differ from other
systems on some cuts, and there may be points where a function is
mathematically well defined but the implementation fails to produce a
finite result because of an artificial internal singularity.
:::

## Complete elliptic integrals {#complete}

::: function
void acb_elliptic_k(acb_t k, const acb_t m, slong prec)

Sets `k` to the complete elliptic integral of the first kind
$K(m) = \int_0^{\pi/2} (1 - m \sin^2 t)^{-1/2} \, dt$, computed as
$K(m) = \pi / (2 M(\sqrt{1-m}))$ where $M(z) = \operatorname{agm}(1, z)$
is computed by `acb_agm1`. The branch cut is $m \in [1, \infty)$.

::: details
$$
K(m) = \int_0^{\pi/2} \frac{dt}{\sqrt{1-m \sin^2 t}}
= \int_0^1 \frac{dt}{\sqrt{1-t^2}\sqrt{1-mt^2}}.
$$

The AGM converges quadratically, so the cost is $O(\log p)$ square roots
at precision $p$; see [](algorithms/special-functions/elliptic-modular).
:::
:::

::: function
void acb_elliptic_k_jet(acb_ptr w, const acb_t m, slong len, slong prec)

Sets `w` to the Taylor coefficients of $K(m + x)$ at $x = 0$, truncated
to length `len`. The first two coefficients come from
$(M, M')$ evaluated by `acb_agm1_cpx`, the others from the linear
recurrence implied by the hypergeometric differential equation of $K$.

::: details
Since $K(m) = \frac{\pi}{2} \, {}_2F_1(\tfrac12, \tfrac12; 1; m)$, the
coefficients $c_k$ of $K(m + x)$ satisfy

$$
4(k+1)(k+2) \, m(m-1) \, c_{k+2} = -(2k+1)^2 c_k - 4(k+1)^2 (2m-1) c_{k+1},
$$

At $m = 0$ it uses
$c_k = \frac{\pi}{2} \left(\frac{(1/2)_k}{k!}\right)^2$ directly.
:::
:::

::: function
void _acb_elliptic_k_series(acb_ptr res, acb_srcptr m, slong mlen, slong len, slong prec)
void acb_elliptic_k_series(acb_poly_t res, const acb_poly_t m, slong len, slong prec)

Set `res` to $K(m(x))$ truncated to length `len`, for a power series
$m(x)$, by composing the output of `acb_elliptic_k_jet` with the
nonconstant part of $m$. The underscore version requires `len` $\ge 1$
and reads $\min(\text{mlen}, \text{len})$ coefficients of `m`.
:::

::: function
void acb_elliptic_e(acb_t res, const acb_t m, slong prec)

Sets `res` to the complete elliptic integral of the second kind
$E(m) = \int_0^{\pi/2} (1 - m \sin^2 t)^{1/2} \, dt$, computed from
$K(m)$ and $K'(m)$ by $E(m) = (1-m)(2m K'(m) + K(m))$, where $K'$ is the
derivative (not the complementary integral).

::: details
$$
E(m) = \int_0^{\pi/2} \sqrt{1-m \sin^2 t} \, dt =
\int_0^1 \frac{\sqrt{1-mt^2}}{\sqrt{1-t^2}} \, dt.
$$

The special values $E(0) = \pi/2$ and $E(1) = 1$ are handled directly.
:::
:::

::: function
void acb_elliptic_pi(acb_t r, const acb_t n, const acb_t m, slong prec)

Sets `r` to the complete elliptic integral of the third kind
$\Pi(n, m) = \int_0^{\pi/2} (1 - n \sin^2 t)^{-1} (1 - m \sin^2 t)^{-1/2} \, dt$.
For $n = 0$ this is $K(m)$ and for $m = 0$ it is
$\frac{\pi}{2}(1-n)^{-1/2}$; otherwise it is computed as
$R_F(0, 1-m, 1) + \frac{n}{3} R_J(0, 1-m, 1, 1-n)$, the same reduction
as for the incomplete integral, which is slower than the AGM methods for
$K$ and $E$.

::: details
$$
\Pi(n, m) = \int_0^{\pi/2}
\frac{dt}{(1-n \sin^2 t) \sqrt{1-m \sin^2 t}} =
\int_0^1
\frac{dt}{(1-nt^2) \sqrt{1-t^2} \sqrt{1-mt^2}}.
$$
:::
:::

## Legendre incomplete elliptic integrals {#legendre}

The incomplete integrals are computed on the strip
$|\operatorname{Re}(\phi)| \le \pi/2$ by the reductions to Carlson's
integrals of [@NIST2012, §19.25(i)], with $s = \sin \phi$ and
$c = \cos \phi$:

$$
F(\phi, m) = s R_F(c^2, 1 - ms^2, 1), \quad
E(\phi, m) = s R_F(c^2, 1-ms^2, 1) - \tfrac{m}{3} s^3 R_D(c^2, 1-ms^2, 1),
$$
$$
\Pi(n, \phi, m) = s R_F(c^2, 1-ms^2, 1) + \tfrac{n}{3} s^3 R_J(c^2, 1-ms^2, 1, 1-ns^2).
$$

Outside the strip the functions are extended by quasi-periodicity:
writing $\phi = \phi' + d \pi$ with $d \in \Z$ and
$|\operatorname{Re}(\phi')| \le \pi/2$, the integral over a full period
is added $d$ times. If $\operatorname{Re}(\phi)$ is a ball that contains
a point $\pi/2 + d \pi$, both neighbouring values of $d$ are tried and the
union of the results is returned; if it is wider than $\pi/2$, the result
is indeterminate.

::: function
void acb_elliptic_f(acb_t res, const acb_t phi, const acb_t m, int times_pi, slong prec)

Sets `res` to the incomplete elliptic integral of the first kind
$F(\phi, m) = \int_0^\phi (1 - m \sin^2 t)^{-1/2} \, dt$, extended by
$F(\phi + d\pi, m) = 2 d K(m) + F(\phi, m)$. With `times_pi` nonzero,
$\phi$ is replaced by $\pi \phi$. In particular $F(\pi/2, m) = K(m)$.

::: details
$$
F(\phi,m) = \int_0^{\phi} \frac{dt}{\sqrt{1-m \sin^2 t}}
= \int_0^{\sin \phi}
\frac{dt}{\sqrt{1-t^2}\sqrt{1-mt^2}}.
$$

The case $m = 0$ returns $\phi$, and the case of an exact half-integer
$\phi$ with `times_pi` returns $2\phi K(m)$.
:::
:::

::: function
void acb_elliptic_e_inc(acb_t res, const acb_t phi, const acb_t m, int times_pi, slong prec)

Sets `res` to the incomplete elliptic integral of the second kind
$E(\phi, m) = \int_0^\phi (1 - m \sin^2 t)^{1/2} \, dt$, extended by
$E(\phi + d\pi, m) = 2 d E(m) + E(\phi, m)$. With `times_pi` nonzero,
$\phi$ is replaced by $\pi \phi$. In particular $E(\pi/2, m) = E(m)$.

::: details
$$
E(\phi,m) = \int_0^{\phi} \sqrt{1-m \sin^2 t} \, dt =
\int_0^{\sin \phi}
\frac{\sqrt{1-mt^2}}{\sqrt{1-t^2}} \, dt.
$$
:::
:::

::: function
void acb_elliptic_pi_inc(acb_t res, const acb_t n, const acb_t phi, const acb_t m, int times_pi, slong prec)

Sets `res` to the incomplete elliptic integral of the third kind
$\Pi(n, \phi, m) = \int_0^\phi (1 - n \sin^2 t)^{-1} (1 - m \sin^2 t)^{-1/2} \, dt$,
extended by $\Pi(n, \phi + d\pi, m) = 2 d \Pi(n, m) + \Pi(n, \phi, m)$.
With `times_pi` nonzero, $\phi$ is replaced by $\pi \phi$. In
particular $\Pi(n, \pi/2, m) = \Pi(n, m)$ and $\Pi(0, \phi, m) = F(\phi, m)$.

::: details
$$
\Pi(n, \phi, m) = \int_0^{\phi}
\frac{dt}{(1-n \sin^2 t) \sqrt{1-m \sin^2 t}} =
\int_0^{\sin \phi}
\frac{dt}{(1-nt^2) \sqrt{1-t^2} \sqrt{1-mt^2}}.
$$

With `times_pi` and a non-real $\phi$, the function multiplies $\phi$ by
$\pi$ and proceeds as without the flag, because the exact argument
reduction can otherwise land exactly on a branch cut of $R_J$.
:::
:::

## Carlson symmetric integrals {#carlson}

Carlson's symmetric forms are the preferred form of incomplete elliptic
integrals: they are symmetric in their arguments, homogeneous, and can be
computed by a duplication algorithm [@Car1995]. Every elliptic integral
can be expressed in terms of $R_F$, $R_G$ and $R_J$, of which
$R_C(x, y) = R_F(x, y, y)$ and $R_D(x, y, z) = R_J(x, y, z, z)$ are
special cases. [](acb_elliptic/carlson) explains the domains of
validity and the choice of method.

::: function
void acb_elliptic_rf(acb_t res, const acb_t x, const acb_t y, const acb_t z, int flags, slong prec)

Sets `res` to the Carlson symmetric integral of the first kind
$R_F(x, y, z) = \frac12 \int_0^\infty ((t+x)(t+y)(t+z))^{-1/2} \, dt$,
by duplication steps followed by a Taylor series in the elementary
symmetric polynomials of the deviations. The result is indeterminate if
more than one argument contains zero.
:::

::: function
void acb_elliptic_rg(acb_t res, const acb_t x, const acb_t y, const acb_t z, int flags, slong prec)

Sets `res` to the Carlson symmetric integral of the second kind $R_G$,
computed from $R_F$ and $R_D$. There are no restrictions on the
arguments.
:::

::: function
void acb_elliptic_rj(acb_t res, const acb_t x, const acb_t y, const acb_t z, const acb_t p, int flags, slong prec)
void acb_elliptic_rj_carlson(acb_t res, const acb_t x, const acb_t y, const acb_t z, const acb_t p, int flags, slong prec)
void acb_elliptic_rj_integration(acb_t res, const acb_t x, const acb_t y, const acb_t z, const acb_t p, int flags, slong prec)

Set `res` to the Carlson symmetric integral of the third kind
$R_J(x, y, z, p) = \frac32 \int_0^\infty (t+p)^{-1} ((t+x)(t+y)(t+z))^{-1/2} \, dt$.
The `carlson` version uses the duplication algorithm, which can give
wrong values for some complex arguments; the `integration` version moves
the arguments to the right half-plane by numerical integration; the
default version uses the duplication algorithm when it is known to be
correct (including the case $p = z$, which computes $R_D(x, y, z)$) and
integration otherwise.
:::

::: function
void acb_elliptic_rc1(acb_t res, const acb_t x, slong prec)

Sets `res` to $R_C(1, 1+x) = \operatorname{atan}(\sqrt{x})/\sqrt{x}$,
a helper for $R_J$, using a short Taylor series for small $|x|$.
:::

## Weierstrass elliptic functions {#weierstrass}

The Weierstrass functions for the lattice $\Z + \tau\Z$ are computed
from the Jacobi theta functions $\theta_i(z, \tau)$ of
[`acb_modular`](acb_modular), which moves $\tau$ to the fundamental
domain and reduces $z$ before summing the theta series. Because the
Weierstrass functions are invariant (or covariant) under the action of
the modular group, the cost is $O(\log p)$ multiplications plus the
theta series, which needs $O(\sqrt{p})$ terms at precision $p$. The
reference is [@NIST2012, Ch. 23]. When $z$ is real and
$\operatorname{Re}(\tau)$ is a multiple of $1/2$, the functions are real
and the imaginary part of the output is set to exactly zero.

::: function
void acb_elliptic_p(acb_t r, const acb_t z, const acb_t tau, slong prec)

Sets `r` to the Weierstrass elliptic function
$\wp(z, \tau) = z^{-2} + \sum_{\omega \in \Lambda \setminus 0} ((z + \omega)^{-2} - \omega^{-2})$,
computed by

$$
\wp(z, \tau) = \pi^2 \theta_2^2(0,\tau) \theta_3^2(0,\tau)
\frac{\theta_4^2(z,\tau)}{\theta_1^2(z,\tau)} -
\frac{\pi^2}{3} \left[ \theta_2^4(0,\tau) + \theta_3^4(0,\tau)\right].
$$
:::

::: function
void acb_elliptic_p_prime(acb_t r, const acb_t z, const acb_t tau, slong prec)

Sets `r` to the derivative $\wp'(z, \tau)$ with respect to $z$, computed
by $\wp'(z,\tau) = (-2\pi \eta(\tau)^2/\theta_1(z,\tau))^3 \, \theta_2(z,\tau) \theta_3(z,\tau) \theta_4(z,\tau)$,
where $\eta$ is the Dedekind eta function.
:::

::: function
void acb_elliptic_p_jet(acb_ptr r, const acb_t z, const acb_t tau, slong len, slong prec)

Sets `r` to the Taylor coefficients of $\wp(z + x, \tau)$ at $x = 0$,
truncated to length `len`, from the theta function jets of
`acb_modular_theta_jet`. With `len` = 2 this gives
$(\wp(z, \tau), \wp'(z, \tau))$, which generate the field of elliptic
functions with periods 1 and $\tau$.
:::

::: function
void _acb_elliptic_p_series(acb_ptr res, acb_srcptr z, slong zlen, const acb_t tau, slong len, slong prec)
void acb_elliptic_p_series(acb_poly_t res, const acb_poly_t z, const acb_t tau, slong len, slong prec)

Set `res` to $\wp(z(x), \tau)$ truncated to length `len`, for a power
series $z(x)$. The underscore version requires `len` $\ge 1$.
:::

::: function
void acb_elliptic_invariants(acb_t g2, acb_t g3, const acb_t tau, slong prec)

Sets `g2` and `g3` to the lattice invariants $g_2 = 60 G_4(\tau)$ and
$g_3 = 140 G_6(\tau)$, where $G_{2k}$ are the Eisenstein series of
`acb_modular_eisenstein`. They are the coefficients of the differential
equation $\wp'^2 = 4 \wp^3 - g_2 \wp - g_3$.
:::

::: function
void acb_elliptic_roots(acb_t e1, acb_t e2, acb_t e3, const acb_t tau, slong prec)

Sets `e1`, `e2`, `e3` to the lattice roots, the roots of
$4x^3 - g_2 x - g_3$, which are $e_1 = \wp(1/2)$,
$e_2 = \wp((1 + \tau)/2)$ and $e_3 = \wp(\tau/2)$, computed from the
fourth powers of the theta constants [@NIST2012, 23.6.2–23.6.4]. The
outputs are indeterminate if $\operatorname{Im}(\tau)$ is not certainly
positive.

::: details
With $\theta_j = \theta_j(0, \tau)$,

$$
e_1 = \frac{\pi^2}{3}(\theta_3^4 + \theta_4^4), \quad
e_2 = \frac{\pi^2}{3}(\theta_2^4 - \theta_4^4), \quad
e_3 = -\frac{\pi^2}{3}(\theta_2^4 + \theta_3^4).
$$

The roots are real when $\operatorname{Re}(\tau)$ is an integer, and
$e_1$ is real when $\operatorname{Re}(\tau)$ is a half-integer; in these
cases the imaginary parts are set to exactly zero.
:::
:::

::: function
void acb_elliptic_inv_p(acb_t res, const acb_t z, const acb_t tau, slong prec)

Sets `res` to the inverse Weierstrass function, a solution $w$ of
$\wp(w, \tau) = z$, given by the elliptic integral

$$
\wp^{-1}(z, \tau) = \frac12 \int_z^\infty \frac{dt}{\sqrt{(t-e_1)(t-e_2)(t-e_3)}}
= R_F(z - e_1, z - e_2, z - e_3).
$$
:::

::: function
void acb_elliptic_zeta(acb_t res, const acb_t z, const acb_t tau, slong prec)

Sets `res` to the Weierstrass zeta function
$\zeta(z, \tau) = z^{-1} + \sum_{\omega \ne 0} ((z - \omega)^{-1} + \omega^{-1} + z \omega^{-2})$,
which satisfies $\zeta' = -\wp$ and is quasi-periodic:
$\zeta(z + 1) = \zeta(z) + 2\eta_1$, $\zeta(z + \tau) = \zeta(z) + 2\eta_3$
with $\eta_1 = \zeta(1/2)$ and $\eta_3 = \zeta(\tau/2)$. It is computed
as $\zeta(z) = 2 \eta_1 z + \theta_1'(z, \tau)/\theta_1(z, \tau)$, where
$\eta_1 = -\theta_1'''(0, \tau) / (6 \theta_1'(0, \tau))$.
:::

::: function
void acb_elliptic_sigma(acb_t res, const acb_t z, const acb_t tau, slong prec)

Sets `res` to the Weierstrass sigma function
$\sigma(z, \tau) = z \prod_{\omega \ne 0} (1 - z/\omega) \exp(z/\omega + z^2/(2\omega^2))$,
which satisfies $\sigma'/\sigma = \zeta$, computed as
$\sigma(z) = e^{\eta_1 z^2} \theta_1(z, \tau) / \theta_1'(0, \tau)$.

::: details
It is quasi-periodic with
$\sigma(z + 1) = -e^{2 \eta_1 (z + 1/2)} \sigma(z)$ and
$\sigma(z + \tau) = -e^{2 \eta_3 (z + \tau/2)} \sigma(z)$.
:::
:::
