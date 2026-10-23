---
title: Bessel, Airy and Coulomb wave functions
summary: How the Bessel functions J, Y, I, K, the Airy functions and the Coulomb wave functions are reduced to the confluent hypergeometric functions 0F1 and U, how the method and the working precision are chosen, and the error bounds of the asymptotic expansions of the Airy functions, with the connection formulas that make them valid in the whole plane.
order: 46
---

The Bessel, Airy and Coulomb wave functions are confluent
hypergeometric functions: each of them is a combination of
${}_0F_1$ or ${}_1F_1$ series near the origin and of $U$ functions with
exponential prefactors near infinity. FLINT evaluates them with the
general methods of [](algorithms/special-functions/hypergeometric),
plus a few specializations: the cancellation estimate
that decides the working precision, propagated error bounds based on
derivative bounds, dedicated error bounds for the Airy functions, and
interval Newton iteration for the zeros of the Airy functions.

## Bessel functions {#bessel}

For all complex $\nu$ and $z$,
$$
J_{\nu}(z) = \frac{1}{\Gamma(\nu+1)} \left(\frac{z}{2}\right)^{\nu}
{}_0F_1\left(\nu+1, -\frac{z^2}{4}\right), \qquad
I_{\nu}(z) = \frac{1}{\Gamma(\nu+1)} \left(\frac{z}{2}\right)^{\nu}
{}_0F_1\left(\nu+1, \frac{z^2}{4}\right),
$$
with the regularized ${}_0F_1$ for negative integer $\nu$, and the
asymptotic expansions come from
$$
J_{\nu}(z) = \frac{z^{\nu}}{2^{\nu} e^{iz} \Gamma(\nu+1)}
{}_1F_1(\nu+\tfrac{1}{2}, 2\nu+1, 2iz) = A_{+} B_{+} + A_{-} B_{-},
$$
$$
A_{\pm} = (\pm iz)^{-1/2-\nu} z^{\nu} (2 \pi)^{-1/2}, \qquad
B_{\pm} = e^{\pm i z} U^{*}(\nu+\tfrac{1}{2}, 2\nu+1, \mp 2iz),
$$
$$
K_{\nu}(z) = \left(\frac{2z}{\pi}\right)^{-1/2} e^{-z}
U^{*}(\nu+\tfrac{1}{2}, 2\nu+1, 2z).
$$
The prefactors $A_\pm$ simplify to $(\pm i)^{n}(2\pi z)^{-1/2}$ when
$\nu + \tfrac12 = n$ is an integer and to $\exp(\mp i
[(2\nu+1)/4]\pi)(2\pi z)^{-1/2}$ when $\operatorname{Re}(z) > 0$; the
general formula is used otherwise. The functions of the second kind
follow from
$$
Y_{\nu}(z) = \frac{\cos(\nu \pi) J_{\nu}(z) - J_{-\nu}(z)}{\sin(\nu \pi)},
\qquad
K_{\nu}(z) = \frac{\pi}{2 \sin(\pi \nu)} \left[
\left(\frac{z}{2}\right)^{-\nu} {}_0\mathbf{F}_1\left(1-\nu, \frac{z^2}{4}\right) -
\left(\frac{z}{2}\right)^{\nu} {}_0\mathbf{F}_1\left(1+\nu, \frac{z^2}{4}\right)
\right]
$$
for non-integer $\nu$. At integer $\nu = n$, $K_n$ is computed as a
limit with power series in $\nu$ (`acb_hypgeom_bessel_k_0f1_series`),
and $Y_n(z) = -\frac{2}{\pi}(i^n K_n(iz) + [\log(iz) - \log z] J_n(z))$.

**Choice of method.** `acb_hypgeom_bessel_j` uses the convergent series
for $|z| < 8$, and the asymptotic expansion for $|z| \ge 2^{64}$ or
$2|z| \ge p$; the asymptotic expansion of $U^*$ reaches about
$(1 + \log 2)/\log 2 \cdot |z| \approx 2.44|z|$ bits for small $\nu$,
and the factor 2 is a safety margin. Between these regions the
convergent series is used with a working precision raised by the
estimated cancellation. The terms of the ${}_0F_1$ series for
$J_\nu(x + yi)$ grow to about $e^{|z|}$ while the function has size
about $e^{|y|}$, so about $(|z| - |y|)/\log 2$ bits are lost; no
precision is added when the terms do not cancel (purely imaginary $z$
with real $\nu > 0$ or integer $\nu$, or $\nu > |z|^2/4$). The modified
functions $I_\nu$ and $K_\nu$ use the series for $|z| < 16$ and for
$2|z| < p$, and the asymptotic expansion otherwise.

**Propagated error.** For inexact $z$ and integer $\nu$, $J_\nu$ is
evaluated at the midpoint and the error is bounded with
`acb_hypgeom_bessel_j_deriv_bound`, which uses for example
$|J_n^{(d)}(z)| \le e^{|\operatorname{Im} z|}$; this avoids the large
overestimation that ball arithmetic produces when a wide ball is pushed
through a series with cancellation. When the resulting enclosure is
still poor, it is intersected with the direct bound for $|J_\nu(z)|$.

**Large orders and integration.** For real $\nu$ and $z$, if the result
has less than half the target accuracy (typically because $\nu$ is
large compared to $z$ and the precision), `acb_hypgeom_bessel_i` and
`acb_hypgeom_bessel_k` retry with rigorous numerical integration of an
integral representation (`arb_hypgeom_bessel_i_integration`,
`arb_hypgeom_bessel_k_integration`).

**Power series.** The Taylor expansion of $J_\nu$ at a point $z_0 \ne
0$ follows from $J_\nu(z_0)$, $J'_\nu(z_0) = \frac12(J_{\nu-1}(z_0) -
J_{\nu+1}(z_0))$ and a five-term recurrence for the higher coefficients
derived from the Bessel differential equation $z^2 y'' + z y' + (z^2 -
\nu^2) y = 0$. This is implemented generically in `gr_bessel_j_jet` and
`_gr_poly_bessel_j_series`, which `arb_hypgeom_bessel_j_series` and
`acb_hypgeom_bessel_j_series` call; the recurrence divides by $z_0^2$,
so $z_0$ must be bounded away from 0 unless it is exactly 0 and $\nu$
is an integer.

## Airy functions {#airy}

The Airy functions are entire solutions of $y'' = zy$. Near the origin
FLINT sums
$$
\operatorname{Ai}(z) = c_1 f(z) - c_2 g(z), \qquad
\operatorname{Bi}(z) = \sqrt3\,(c_1 f(z) + c_2 g(z)),
$$
$$
f(z) = {}_0F_1(;\tfrac23; \tfrac{z^3}{9}), \qquad
g(z) = z\, {}_0F_1(;\tfrac43; \tfrac{z^3}{9}),
$$
with $c_1 = \operatorname{Ai}(0) = 3^{-2/3}/\Gamma(2/3)$ and $c_2 =
-\operatorname{Ai}'(0) = 3^{-1/3}/\Gamma(1/3)$ (cached constants), by
rectangular splitting in powers of $z^3$ with word-size coefficients
(`acb_hypgeom_airy_direct`). With $w = z^3/9$ the terms are bounded by
$3|w|^n/((n-1)!)^2$, which gives the tail bound. Expressing the Airy
functions through Bessel functions of order $\pm 1/3$ is avoided,
since such formulas hold only in sectors (the Stokes phenomenon).

### Asymptotic expansions {#airy-asymp}

The error bounds are based on Olver [@NIST2012, §9.7]. For $|\arg z| <
\pi$ and $\zeta = \frac23 z^{3/2}$,
$$
\operatorname{Ai}(z) = \frac{e^{-\zeta}}{2 \sqrt{\pi} z^{1/4}}
\left[S_n(\zeta) + R_n(z)\right], \qquad
\operatorname{Ai}'(z) = -\frac{z^{1/4} e^{-\zeta}}{2 \sqrt{\pi}}
\left[S'_n(\zeta) + R'_n(z)\right],
$$
$$
S_n(\zeta) = \sum_{k=0}^{n-1} (-1)^k \frac{u(k)}{\zeta^k}, \quad
S'_n(\zeta) = \sum_{k=0}^{n-1} (-1)^k \frac{v(k)}{\zeta^k}, \quad
u(k) = \frac{(1/6)_k (5/6)_k}{2^k k!}, \quad
v(k) = \frac{6k+1}{1-6k} u(k).
$$
For $n \ge 1$ the remainders are bounded by
$$
|R_n(z)|  \le C |u(n)| |\zeta|^{-n}, \qquad |R'_n(z)| \le C |v(n)|
|\zeta|^{-n},
$$
$$
C = \begin{cases}
2 \exp(7 / (36 |\zeta|)) & |\arg(z)| \le \pi/3, \\
2 \chi(n) \exp(7 \pi / (72 |\zeta|)) & \pi/3 \le |\arg(z)| \le 2\pi/3, \\
4 \chi(n) \exp(7 \pi / (36 |\operatorname{Re}(\zeta)|))
|\cos(\arg(\zeta))|^{-n} & 2\pi/3 \le |\arg(z)| < \pi,
\end{cases}
$$
with $\chi(n) = \sqrt\pi\,\Gamma(\frac n2 + 1)/\Gamma(\frac n2 +
\frac12)$.

For $\operatorname{Bi}$ when $z$ is roughly in the right half-plane, the
connection formulas
$$
\operatorname{Bi}(z) = -i (2 w \operatorname{Ai}(z w^{-2}) -
\operatorname{Ai}(z)), \qquad
\operatorname{Bi}(z) = +i (2 w^{-1} \operatorname{Ai}(z w^{2}) -
\operatorname{Ai}(z)), \qquad w = e^{\pi i/3},
$$
are combined with the expansions, which after collecting the roots of
unity gives
$$
\operatorname{Bi}(z) = \frac{1}{2 \sqrt{\pi} z^{1/4}} [2X \pm iY], \qquad
X = e^{\zeta} [S_n(-\zeta) + R_n(z w^{\mp 2})], \quad
Y = e^{-\zeta} [S_n(\zeta) + R_n(z)],
$$
where the upper signs are valid for $-\pi/3 < \arg(z) < \pi$ and the
lower signs for $-\pi < \arg(z) < \pi/3$; the derivative is analogous.
In the left half-plane the formulas
$$
\operatorname{Ai}(z) = e^{\pi i/3} \operatorname{Ai}(z_1) + e^{-\pi
i/3} \operatorname{Ai}(z_2), \qquad
\operatorname{Bi}(z) = e^{-\pi i/6} \operatorname{Ai}(z_1) + e^{\pi
i/6} \operatorname{Ai}(z_2),
$$
with $z_1 = -z e^{\pi i/3}$ and $z_2 = -z e^{-\pi i/3}$, reduce to
arguments with $|\arg(z_1)|, |\arg(z_2)| < \pi$ when $|\arg(-z)| <
2\pi/3$, and collecting roots of unity gives, with $\zeta = \frac23
(-z)^{3/2}$,
$$
\operatorname{Ai}(z) = A_1 [S_n(i \zeta)  + R_n(z_1)] + A_2 [S_n(-i
\zeta) + R_n(z_2)], \qquad
\operatorname{Bi}(z) = -i A_1 [S_n(i \zeta)  + R_n(z_1)] + i A_2
[S_n(-i \zeta) + R_n(z_2)],
$$
$$
A_1 = \frac{\exp(-i (\zeta - \pi/4))}{2 \sqrt{\pi} (-z)^{1/4}}, \qquad
A_2 = \frac{\exp(+i (\zeta - \pi/4))}{2 \sqrt{\pi} (-z)^{1/4}}.
$$
The differentiated formulas are analogous. `acb_hypgeom_airy_asymp`
implements these formulas, and `acb_hypgeom_airy_bound` evaluates them
with $n = 1$ in `mag_t` arithmetic to get cheap upper bounds for the
four functions on a ball.

### Choice of method {#airy-choice}

`acb_hypgeom_airy` uses the asymptotic expansion when $|z| \ge 4$ and
the terms decrease below $2^{-p}$ (the number of terms comes from an
estimate of $\log_2 |u(k)/\zeta^k|$), and the series otherwise. For the
series it estimates in double precision the magnitude of the largest
term, about $\exp(2|z^3/9|^{1/2})$, and the magnitude of the functions
from the leading asymptotic term, and adds the difference to the
working precision (this ignores the extra cancellation near the zeros).
The series is evaluated at the midpoint of $z$; the error propagated
from the radius $r$ is bounded by $r$ times the bounds of
`acb_hypgeom_airy_bound` for the derivatives, using
$\operatorname{Ai}''(z) = z \operatorname{Ai}(z)$.

### Zeros {#airy-zeros}

`arb_hypgeom_airy_zero` computes the $n$-th real zeros $a_n, a'_n,
b_n, b'_n$ of $\operatorname{Ai}, \operatorname{Ai}', \operatorname{Bi},
\operatorname{Bi}'$ from the asymptotic expansions [@NIST2012, §9.9]
$a_n \sim -T(\tfrac38\pi(4n-1))$, $b_n \sim -T(\tfrac38\pi(4n-3))$,
$a'_n \sim -U(\tfrac38\pi(4n-3))$, $b'_n \sim -U(\tfrac38\pi(4n-1))$,
where $T(t) = t^{2/3}(1 + \frac{5}{48}t^{-2} - \cdots)$ and $U(t) =
t^{2/3}(1 - \frac{7}{48}t^{-2} + \cdots)$. Pittaluga and Sacripante
[@PS1991] proved that two consecutive truncations of the expansions for
$a_n$ and $b_n$ bracket the zero; for $a'_n$ and $b'_n$ the truncated
expansion gives a starting interval that is verified by an interval
Newton step. For $n \le 10$ tabulated 48-bit enclosures are used,
since the expansions are poor for the first zero. The enclosure is
then refined by interval Newton steps $x \mapsto m -
f(m)/f'(X)$ with doubling precision, using $f'' = xf$ for the
derivatives.

## Coulomb wave functions {#coulomb}

The Coulomb wave functions $F_\ell(\eta, z)$, $G_\ell(\eta, z)$ and
$H^\pm_\ell(\eta,z) = G_\ell \pm i F_\ell$ solve $y'' + (1 - 2\eta/z -
\ell(\ell+1)/z^2) y = 0$. With the normalization constant
$$
C_\ell(\eta) = \frac{2^\ell e^{-\pi\eta/2} |\Gamma(\ell + 1 + i\eta)|}{\Gamma(2\ell+2)}
$$
(computed through $\log\Gamma(1+\ell \pm i\eta)$ so that it is analytic
in $\ell$ and $\eta$), the regular solution is
$$
F_\ell(\eta, z) = C_\ell(\eta)\, z^{\ell+1} e^{\mp iz}\,
{}_1F_1(1 + \ell \mp i\eta, 2\ell+2, \pm 2iz),
$$
and all four functions are expressed through $U(1+\ell \pm i\eta,
2\ell+2, \mp 2iz)$ with exponential and power prefactors. The formulas,
with the branch cuts placed on the negative real $z$ axis, follow the
connection formulas collected in [@Fungrim] (see also [@DYF1999;
@Gas2018; @Mic2007] and [@NIST2012, chapter 33]); different formulas
are used in the upper and lower half-planes and on the real line so
that each function is computed without unnecessary cancellation, and
$H^\pm$ are computed directly rather than as $G \pm iF$, which can
cancel. The $U$ functions are evaluated with the asymptotic expansion
when $|z|$ is large enough and with `acb_hypgeom_u` otherwise. The jets
(`acb_hypgeom_coulomb_jet`) use $F$ and $G$ at $\ell$ and $\ell + 1$ for
the first derivative and the differential equation for the higher
coefficients.

## Use in FLINT {#use}

- Bessel functions: `acb_hypgeom_bessel_j`, `acb_hypgeom_bessel_y`,
  `acb_hypgeom_bessel_jy`, `acb_hypgeom_bessel_i`,
  `acb_hypgeom_bessel_k` and their scaled versions, with the methods
  `acb_hypgeom_bessel_j_0f1`, `acb_hypgeom_bessel_j_asymp`,
  `acb_hypgeom_bessel_i_0f1`, `acb_hypgeom_bessel_i_asymp`,
  `acb_hypgeom_bessel_k_0f1`, `acb_hypgeom_bessel_k_asymp`; the real
  wrappers in [`arb_hypgeom`](arb_hypgeom); the generic
  `gr_bessel_j_jet` and `_gr_poly_bessel_j_series`.
- ${}_0F_1$: `acb_hypgeom_0f1` uses the Bessel connection formulas for
  its asymptotic expansion.
- Airy functions: `acb_hypgeom_airy`, `acb_hypgeom_airy_direct`,
  `acb_hypgeom_airy_asymp`, `acb_hypgeom_airy_bound`,
  `acb_hypgeom_airy_jet`, the series versions and
  `arb_hypgeom_airy_zero`.
- Coulomb wave functions: `acb_hypgeom_coulomb`,
  `acb_hypgeom_coulomb_jet`, `acb_hypgeom_coulomb_series` and the real
  versions.
