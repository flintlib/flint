---
title: Elliptic integrals, elliptic functions and modular forms
summary: The complex arithmetic-geometric mean with Dupont's error bound, complete and incomplete elliptic integrals through the AGM and Carlson's duplication algorithm with high-order Taylor series, Weierstrass functions via theta functions, and modular forms by reduction to the fundamental domain and q-series with short addition sequences.
complexity: $O(\log p)$ square roots at precision $p$ for the AGM and complete elliptic integrals; $O(p^{0.6})$ duplication steps plus a Taylor series of $O(p^{0.4})$ terms for Carlson's integrals; $O(\sqrt{p})$ multiplications for modular forms and theta functions on the fundamental domain.
order: 60
---

Elliptic integrals and elliptic functions reduce to a few primitives
that converge quickly: the arithmetic-geometric mean, which gives the
complete elliptic integrals; Carlson's duplication algorithm, which gives
the incomplete integrals; and theta functions, which give the
Weierstrass elliptic functions and the modular forms after a modular
transformation has moved $\tau$ to the fundamental domain. The sections
below describe these methods as implemented in [`acb`](acb) (the AGM),
[`acb_elliptic`](acb_elliptic) and [`acb_modular`](acb_modular), following
[@Joh2019]; the summation of theta series is described in
[](algorithms/special-functions/theta-functions).

## The arithmetic-geometric mean {#agm}

For $a, b > 0$, the arithmetic-geometric mean $\operatorname{agm}(a, b)$
is the common limit of $a_{n+1} = (a_n + b_n)/2$,
$b_{n+1} = \sqrt{a_n b_n}$, with $a_0 = a$, $b_0 = b$. The convergence is
quadratic: $|a_{n+1} - b_{n+1}| \le |a_n - b_n|^2/(8 b_{n+1})$, so
$O(\log p)$ steps give $p$ bits once $a_n$ and $b_n$ agree to a few bits.
With complex variables, it is convenient to use the univariate function
$M(z) = \operatorname{agm}(1, z)$; then $\operatorname{agm}(a,b) = a M(b/a)$
(`acb_agm`), and $M$ is analytic on $\mathbb{C} \setminus (-\infty, 0]$.

**Functional equation.** If the real part of $z$ is not nonnegative,
`acb_agm1` first applies $M(z) = \frac{z+1}{2} M(u)$ with
$u = 2\sqrt{z}/(z+1)$, which has nonnegative real part (rounding errors
that make the ball for $u$ cross the imaginary axis only inflate the
result). For the derivative this becomes
$M'(z) = [M(u) - (z-1) M'(u) / ((1+z) \sqrt{z})] / 2$.

**Square roots.** In the iteration the correct square root is the
continuous branch $\sqrt{a_n}\sqrt{b_n}$, which is computed as
$\sqrt{a_n b_n}$, $i\sqrt{-a_n b_n}$ or $-i\sqrt{-a_n b_n}$ when $a_n$
and $b_n$ both have positive real part, both have nonnegative imaginary
part, or both have nonpositive imaginary part, and as
$\sqrt{a_n}\sqrt{b_n}$ otherwise.

::: lemma Dupont [@Dup2006, p. 87] {#lem-dupont}
If $\operatorname{Re}(z) \ge 0$, the iterates of the AGM of $(1, z)$
satisfy $|M(z) - a_n| \le |a_n - b_n|$ for all $n$.
:::

This gives a rigorous enclosure when the iteration is stopped, and for
real positive variables the lower and upper bounds $b_n \le M \le a_n$ can
be used directly (`arb_agm`). It is slightly faster to stop when $a_n$ and
$b_n$ agree to about $p/10$ bits and finish with a series: with
$w = (a - b)/(a + b)$,

$$
\operatorname{agm}(a,b) = \frac{(a+b) \pi}{4 K(w^2)}, \qquad
\frac{\pi}{4 K(w^2)} = \frac{1}{2} - \frac{w^2}{8} - \frac{5 w^4}{128} - \frac{11 w^6}{512} - \frac{469 w^8}{32768} + \ldots,
$$

valid for $|w| < 1$ and $a, b$ with nonnegative real part, where the
tail is bounded by $\sum_{k \ge 10} |w|^k/64$.

**First derivative.** `acb_agm1_cpx` computes $(M(z), M'(z))$ for exact
$z$ with $|\arg z| \le 3\pi/4$ by a central finite difference. Since
$|M(z)| \le \max(1, |z|)$, Cauchy's integral formula gives
$|M^{(k)}(z)/k!| \le C D^k$ with $C = \max(1, |z| + r)$ and $D = 1/r$ for
$0 < r < |z|$ (FLINT takes $r \approx |z|/4$). For $hD < 1$,

$$
\left|\frac{M(z+h) - M(z-h)}{2h} - M'(z)\right| \le \frac{C D^3 h^2}{1 - D h},
\qquad
\left|\frac{M(z+h) + M(z-h)}{2} - M(z)\right| \le \frac{C D^2 h^2}{1 - D h},
$$

so two evaluations at about 1.5 times the precision give both values,
against twice the precision for a one-sided difference. For inexact $z$,
the midpoint is used and the propagated error is bounded with
$|M'(z+\varepsilon)| \le \max(1, |z|+|\varepsilon|+r)/r$ and
$|M''(z+\varepsilon)| \le 2\max(1, |z|+|\varepsilon|+r)/r^2$, valid if
the disk of radius $|\varepsilon| + r$ around $z$ avoids the negative
real axis (with $r \approx |z|/2$).

**Higher derivatives.** $W(z) = 1/M(z) = \frac{2}{\pi} K(1 - z^2)$ is
D-finite. The coefficients of
$W(z+x) = \sum_{k} c_k x^k$ satisfy $-2 z (z^2-1) c_2 = (3z^2-1) c_1 + z c_0$
and

$$
-(k+2)(k+3) z (z^2-1) c_{k+3} = (k+2)^2 (3z^2-1) c_{k+2} + (3k(k+3)+7)z c_{k+1} + (k+1)^2 c_{k},
$$

which reduces at $z = 1$ to
$-(k+2)^2 c_{k+2} = (3k(k+3)+7) c_{k+1} + (k+1)^2 c_{k}$. Given
$(M(z), M'(z))$, the remaining coefficients follow by the recurrence and
one power series inversion.

## Complete elliptic integrals {#complete}

With $M$ as above,

$$
K(m) = \frac{\pi}{2 M(\sqrt{1-m})}, \qquad
E(m) = (1-m)\bigl(2m K'(m) + K(m)\bigr),
$$

where $K'$ is the derivative; both cost $O(\log p)$ square roots
(`acb_elliptic_k`, `acb_elliptic_e`). The Taylor coefficients of
$K(m + x)$ (`acb_elliptic_k_jet`) satisfy the recurrence of the
hypergeometric equation of
$K = \frac{\pi}{2}\,{}_2F_1(\frac12, \frac12; 1; m)$:

$$
4(k+1)(k+2)\, m(m-1)\, c_{k+2} = -(2k+1)^2 c_k - 4(k+1)^2 (2m-1) c_{k+1}.
$$

The complete integral of the third kind, $\Pi(n, m)$, is computed with
Carlson's integrals,
$\Pi(n, m) = R_F(0, 1-m, 1) + \frac{n}{3} R_J(0, 1-m, 1, 1-n)$.

## Carlson's symmetric integrals {#carlson}

The symmetric integrals $R_F$, $R_J$ (with $R_C$, $R_D$ as special cases)
and $R_G$ are defined on [](acb_elliptic/carlson). All three Legendre
integrals reduce to them on the strip $|\operatorname{Re}\phi| \le \pi/2$
[@NIST2012, §19.25], and are extended outside it by quasi-periodicity.

**Duplication.** For $\lambda = \sqrt{x}\sqrt{y} + \sqrt{y}\sqrt{z} + \sqrt{z}\sqrt{x}$,

$$
R_F(x,y,z) = R_F\left(\frac{x+\lambda}{4}, \frac{y+\lambda}{4}, \frac{z+\lambda}{4}\right),
$$

and the differences between the arguments are divided by about 4 at each
step, while their mean $A = (x+y+z)/3$ is preserved up to the factor
introduced by the step. After $j$ steps, with
$X = 1 - x/A$, $Y = 1 - y/A$, $Z = -X - Y$, the elementary symmetric
polynomials $E_2 = XY - Z^2$ and $E_3 = XYZ$ of the deviations are of
order $4^{-2j}$ and $4^{-3j}$, and

$$
R_F(x,y,z) = A^{-1/2} \sum_{a, b \ge 0} \frac{(-1)^a \, (\frac12)_{a+b}}{a!\, b!\, (2(2a+3b)+1)} E_2^a E_3^b.
$$

Carlson [@Car1995] uses the terms of total degree $2a + 3b \le 7$ and
enough duplication steps to make the remainder negligible, which costs
$O(p)$ square roots at precision $p$. FLINT instead sums the series up
to degree $N \approx 2.1\,p^{0.4}$ (real case) or $2.5\,p^{0.4}$ (complex
case), at most 500, so that only $O(p/N)$ duplication steps are needed
[@Joh2019]: the iteration stops as soon as
$(\max(|x - y|, |x - z|)/|x|)^N < 2^{-p}$. The coefficients are computed
exactly with a common denominator, the polynomial in $E_2$ is evaluated
with precomputed powers and the polynomial in $E_3$ by Horner's rule. The
tail is bounded using that the coefficient of $X^{n_1} Y^{n_2} Z^{n_3}$
with $n_1 + n_2 + n_3 = n$ in the expansion has absolute value at most
$2 (9/8)^n$, so that the truncation error is at most
$2 \sum_{n \ge N} (\tfrac98 \max(|X|, |Y|, |Z|))^n$.

**$R_J$.** The duplication theorem for $R_J$ adds a term
$6 \cdot 4^{-j} R_C(1, 1 + \delta_j)$ (up to scaling) at each step, and
the final Taylor series is in the symmetric polynomials $E_2, \ldots, E_5$
of four deviations [@Car1995]; $R_C(1, 1+x) = \operatorname{atan}(\sqrt{x})/\sqrt{x}$
is computed by `acb_elliptic_rc1`. The principal square roots in the
duplication can cross branch cuts of $R_J$ for complex arguments, so the
algorithm is only used where it is known to be valid, and otherwise
$R_J$ is shifted into the right half-plane by numerical integration; see
[](acb_elliptic/carlson#rj). For $R_D$ the duplication algorithm is valid
everywhere.

**$R_G$** is computed from
$2R_G(x,y,z) = z R_F(x,y,z) - \frac13 (x-z)(y-z) R_D(x,y,z) + \sqrt{xy/z}$.

## Weierstrass elliptic functions {#weierstrass}

For the lattice $\Z + \tau\Z$ with half-periods $\frac12$ and
$\frac\tau2$, the Weierstrass functions are quotients of theta functions
[@NIST2012, §23.6]. With $\theta_j = \theta_j(0, \tau)$ and
$\theta_j(z) = \theta_j(z, \tau)$ in the notation of
[`acb_modular`](acb_modular),

$$
\wp(z, \tau) = \pi^2 \theta_2^2 \theta_3^2 \frac{\theta_4^2(z)}{\theta_1^2(z)} - \frac{\pi^2}{3} (\theta_2^4 + \theta_3^4),
\qquad
\wp'(z, \tau) = \left(\frac{-2\pi \eta(\tau)^2}{\theta_1(z)}\right)^3 \theta_2(z)\theta_3(z)\theta_4(z),
$$

$$
\zeta(z, \tau) = 2\eta_1 z + \frac{\theta_1'(z)}{\theta_1(z)}, \qquad
\sigma(z, \tau) = e^{\eta_1 z^2} \frac{\theta_1(z)}{\theta_1'(0)}, \qquad
\eta_1 = -\frac{\theta_1'''(0)}{6 \theta_1'(0)},
$$

where $\eta(\tau)$ is the Dedekind eta function. The lattice roots are
$e_1 = \frac{\pi^2}{3}(\theta_3^4 + \theta_4^4)$,
$e_2 = \frac{\pi^2}{3}(\theta_2^4 - \theta_4^4)$ and
$e_3 = -\frac{\pi^2}{3}(\theta_2^4 + \theta_3^4)$, the invariants are
$g_2 = 60 G_4$ and $g_3 = 140 G_6$, and the inverse function is the
elliptic integral $\wp^{-1}(z) = R_F(z - e_1, z - e_2, z - e_3)$.

All of these need theta values at one $\tau$, and the theta functions are
evaluated after reducing $\tau$ to the fundamental domain, so the cost is
$O(\sqrt{p})$ multiplications for every $\tau$.

## Modular forms {#modular}

A modular form of weight $k$ satisfies
$f(g\tau) = (c\tau + d)^k f(\tau)$ for
$g = \begin{pmatrix} a & b \\ c & d\end{pmatrix} \in \mathrm{SL}_2(\Z)$,
so it can be evaluated at the image $\tau'$ of $\tau$ in the fundamental
domain $\mathcal{F} = \{|\operatorname{Re}\tau| \le \frac12, |\tau| \ge 1\}$,
where $\operatorname{Im}(\tau') \ge \sqrt{3}/2$ and the $q$-series
converge geometrically with ratio $|q| \le e^{-\pi\sqrt 3} \approx 0.0043$
(for $q = e^{2\pi i \tau}$).

::: algorithm Reduction to the fundamental domain {#alg-reduce}
**Input:** $\tau$ in the upper half-plane, a tolerance $\varepsilon$.

**Output:** $g \in \mathrm{PSL}_2(\Z)$ with $g\tau$ in $\mathcal{F}$ up to
$\varepsilon$.

1. $g \gets 1$.
2. Repeat: let $n$ be the nearest integer to $\operatorname{Re}(\tau)$,
   $\tau \gets \tau - n$ and $g \gets T^{-n} g$; if $|\tau|^2 \ge 1 - \varepsilon$,
   return $g$; otherwise $\tau \gets -1/\tau$ and $g \gets S g$.
:::

Each inversion multiplies $\operatorname{Im}(\tau)$ by $1/|\tau|^2 > 1$,
which guarantees termination. `acb_modular_fundamental_domain_approx` runs the algorithm on
the midpoint in `double` arithmetic when possible and in `arf`
arithmetic at the working precision otherwise, and then applies $g$ to
the ball $\tau$ once, which does not need to be exact: reducing to
$|\tau'|^2 \ge 63/64$ instead of 1 is enough for fast convergence.

**Transformation formulas.** The eta function transforms as
$\eta(g\tau) = \varepsilon(a,b,c,d)\sqrt{c\tau + d}\,\eta(\tau)$ with a
24th root of unity $\varepsilon$ given by Dedekind sums
(`acb_modular_epsilon_arg`), and the four theta functions are permuted
with 8th roots of unity and the factors of
`acb_modular_theta_transform` [@Rad1973].

**$q$-series.** On $\mathcal{F}$, the theta constants and the eta
function are lacunary series,

$$
\eta(\tau) = e^{\pi i \tau/12} \sum_{k \in \Z} (-1)^k q^{k(3k-1)/2}, \qquad
\theta_3(0, \tau) = \sum_{k \in \Z} q_1^{k^2},
$$

(with $q = e^{2\pi i\tau}$ and $q_1 = e^{\pi i\tau}$), which need
$O(\sqrt{p})$ terms. The powers of $q$ with exponents in the sequence of
generalised pentagonal numbers (or squares and pronic numbers) are
computed with an addition sequence in which most terms cost one
multiplication; Enge, Hart and Johansson [@EHJ2016] construct short
addition sequences for these exponents
(`acb_modular_addseq_eta`, `acb_modular_addseq_theta`). For many terms,
rectangular splitting is used instead: the series is written as a
polynomial in $q^m$ for a well-chosen $m$ such that the exponents modulo
$m$ take few values, which replaces most full multiplications by
cheaper additions and scalar operations. Each term is computed at a precision reduced according
to its magnitude (`_acb_modular_mul`).

**Forms computed from theta constants.** With $\theta_j = \theta_j(0,\tau)$
at the reduced point:

$$
j(\tau) = 32 \frac{(\theta_2^8+\theta_3^8+\theta_4^8)^3}{(\theta_2 \theta_3 \theta_4)^8}, \qquad
\Delta(\tau) = \eta(\tau)^{24} = \left(\frac{\theta_2\theta_3\theta_4}{2}\right)^8, \qquad
\lambda(\tau) = \frac{\theta_2^4}{\theta_3^4},
$$

$$
G_4(\tau) = \frac{\pi^4}{90}(\theta_2^8 + \theta_3^8 + \theta_4^8), \qquad
G_6(\tau) = \frac{\pi^6}{945}\bigl(\theta_3^{12} + \theta_4^{12} - 3\theta_2^8(\theta_3^4 + \theta_4^4)\bigr),
$$

and the higher Eisenstein series follow from the recurrence
$(2k+1)(k-3) c_k = 3 \sum_{j=2}^{k-2} c_j c_{k-j}$ for
$c_k = (2k-1) G_{2k}$ [@Apostol1997, Theorem 1.13]. $j$ is invariant
under $\mathrm{SL}_2(\Z)$ and $\lambda$ under $\Gamma(2)$ (for $\lambda$,
the transformation of the theta functions is applied).

**Class polynomials.** The Hilbert class polynomial $H_D$ is the product
of $x - j(\tau_i)$ over the CM points $\tau_i = (-b + \sqrt{D})/(2a)$ of
the reduced forms of discriminant $D$. Its coefficients have about
$\pi\sqrt{|D|}\sum_i 1/a_i$ bits; computing the roots to that precision
and expanding the product with a binary product tree gives the
polynomial with a certificate when every coefficient ball contains a
unique integer [@Enge2009]. See [](acb_modular/class-polynomials).

## Use in FLINT {#use}

- `acb_agm1`, `acb_agm1_cpx`, `acb_agm`, `arb_agm` ([`acb`](acb),
  [`arb`](arb)): the AGM.
- `acb_elliptic_k`, `acb_elliptic_k_jet`, `acb_elliptic_e`: complete
  integrals by the AGM; `acb_elliptic_pi`: via $R_F$ and $R_J$.
- `acb_elliptic_rf`, `acb_elliptic_rj_carlson`, `acb_elliptic_rj`,
  `acb_elliptic_rj_integration`, `acb_elliptic_rg`, `acb_elliptic_rc1`:
  Carlson's integrals; `acb_elliptic_f`, `acb_elliptic_e_inc`,
  `acb_elliptic_pi_inc`: Legendre forms.
- `acb_elliptic_p`, `acb_elliptic_p_prime`, `acb_elliptic_p_jet`,
  `acb_elliptic_zeta`, `acb_elliptic_sigma`, `acb_elliptic_roots`,
  `acb_elliptic_invariants`, `acb_elliptic_inv_p`: Weierstrass
  functions.
- `acb_modular_fundamental_domain_approx`, `acb_modular_theta`,
  `acb_modular_eta`, `acb_modular_j`, `acb_modular_delta`,
  `acb_modular_lambda`, `acb_modular_eisenstein`,
  `acb_modular_hilbert_class_poly`: modular forms.
