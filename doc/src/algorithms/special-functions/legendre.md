---
title: Legendre functions and orthogonal polynomials
summary: Legendre polynomials of large degree by the method of Johansson and Mezzarobba (hypergeometric series at 0 and 1, an asymptotic expansion with explicit error bounds, and a recurrence, chosen by cost and cancellation), Gauss–Legendre nodes and weights, associated Legendre functions of both kinds through 2F1 and connection formulas, spherical harmonics, and the classical orthogonal polynomials and their extensions to non-integer degree.
order: 48
---

The Legendre polynomials $P_n(x)$, and more generally the classical
orthogonal polynomials (Jacobi, Gegenbauer, Chebyshev, Laguerre,
Hermite), are terminating hypergeometric series, for example
$$
P_n(x) = {}_2F_1\left(-n, n+1; 1; \frac{1-x}{2}\right), \qquad
P_n^{(a,b)}(z) = \frac{(a+1)_n}{n!}\,{}_2F_1\left(-n,n+a+b+1;a+1;\frac{1-z}{2}\right).
$$
For integer degree they can be evaluated directly as polynomials, but
for large $n$ this is slow and numerically unstable: the coefficients
grow exponentially while $|P_n(x)| \le 1$ on $[-1,1]$. FLINT has a
dedicated algorithm for $P_n(x)$ with real $x$ and word-size $n$, which
also gives the nodes and weights of Gauss–Legendre quadrature; the
other functions use the general ${}_2F_1$ and ${}_1F_1$ code with a few
direct methods for integer degrees.

## Legendre polynomials {#legendre-poly}

The algorithm of `arb_hypgeom_legendre_p_ui` is that of Johansson and
Mezzarobba [@JM2018]. It computes $P_n(x)$ and $P_n'(x)$ together for
real $x$, and by $P_n(-x) = (-1)^n P_n(x)$ only $x \ge 0$ is needed. Four
expansions are available.

- **Expansion at $x = 1$** (`arb_hypgeom_legendre_p_ui_one`):
  $$
  P_n(x) = \sum_{k=0}^{n} \binom{n}{k}\binom{n+k}{k}
  \left(\frac{x-1}{2}\right)^k.
  $$
  The terms decrease quickly when $|1 - x|$ is small, and the series
  can be truncated after $K < n+1$ terms with a bound for the
  remainder.
- **Expansion at $x = 0$** (`arb_hypgeom_legendre_p_ui_zero`): for
  even $n = 2m$,
  $$
  P_{2m}(x) = (-1)^m \frac{(1/2)_m}{m!}\, {}_2F_1\left(-m, m +
  \tfrac12; \tfrac12; x^2\right),
  $$
  and an analogous formula with $x\,{}_2F_1(-m, m + \frac32; \frac32;
  x^2)$ for odd $n$. The terms decrease when $|x|$ is small.
- **Asymptotic expansion** (`arb_hypgeom_legendre_p_ui_asymp`): with
  $x = \cos\theta$, $0 < \theta < \pi$, the Stieltjes series
  $$
  P_n(\cos\theta) = \frac{4}{\pi} \frac{\Gamma(n+1)}{\Gamma(n+\frac32)}
  \sum_{k=0}^{K-1} \frac{(\frac12)_k^2}{k!\,(n+\frac32)_k}
  \frac{\cos\left((n + k + \tfrac12)\theta - (k+\tfrac12)\tfrac{\pi}{2}\right)}{(2\sin\theta)^{k+1/2}}
  + R_K,
  $$
  summed as the real part of a complex series in $(1 - i\cot\theta)$,
  with the explicit bound $|R_K| \le 0.50795\, \frac{K!\, n!}{2^K
  (n+K)!\, y^K \sqrt{ny}}$, $y = \sin\theta$, which follows from the
  analysis of Bogaert, Michiels and Fostier [@Bog2012]. The number of
  terms needed decreases as $n$ grows; for $n$ much larger than the
  precision a handful of terms suffice, independently of $n$.
- **Recurrence** (`arb_hypgeom_legendre_p_ui_rec`): the three-term
  recurrence $(k+1)P_{k+1} = (2k+1)xP_k - kP_{k-1}$, evaluated in
  fixed-point arithmetic with machine words and a rigorous bound for
  the accumulated rounding error. It is stable on $(-1, 1)$, costs
  $O(n)$ cheap operations, and is the fastest method for moderate $n$
  and precision.

`arb_hypgeom_legendre_p_ui` selects the method from estimates of the
cost of each expansion: the number of terms, the precision lost to
cancellation in the series at 0 and at 1 (whose terms can be much
larger than $P_n(x)$), and the size of the error term of the
asymptotic expansion. The recurrence is used for $n$ below a few
hundred at low precision, unless $x$ is very close to 0 or 1; for
$|x| > 1$ only the expansion at 1 is used (with an early stop if it
converges quickly). The working precision is raised by the estimated
cancellation, and for inexact $x$ the function evaluates at the
midpoint and adds the error propagated from the radius $r$, bounded by
$r \cdot \max|P_n'|$ and $r \cdot \max|P_n''|$ on the ball
(`arb_hypgeom_legendre_p_ui_deriv_bound`, which uses $|P_n'(x)| \le
\min((n+1)^2/2,\ c\sqrt{n}\,(1 - x^2)^{-3/4})$ and $|P_n''(x)| \le
\min((n+1)^4/8,\ 2c\, n^{3/2}(1-x^2)^{-5/4})$ for $|x| \le 1$, with $c
\ge 2^{3/2}/\sqrt\pi$).

## Gauss–Legendre nodes and weights {#roots}

The roots $1 > x_0 > \cdots > x_{n-1} > -1$ of $P_n$ and the weights
$w_k = 2/((1-x_k^2)P_n'(x_k)^2)$ are computed by
`arb_hypgeom_legendre_p_ui_root`: an initial approximation with an
explicit error bound due to Petras [@Pet1999], Newton iteration with a
priori error bounds at low precision, and interval Newton iteration at
doubling precisions above; each step evaluates $P_n$ and $P_n'$ with
the method above, and the doubling of the precision means that the
final evaluation dominates the cost. This is what makes high-degree
Gauss–Legendre quadrature in `acb_calc_integrate` practical at high
precision [@JM2018]. See
[](arb_hypgeom/legendre-roots).

## Legendre functions {#legendre-functions}

For complex degree $n$, order $m$ and argument $z$, the associated
Legendre function of the first kind is
$$
P_n^m(z) = \left(\frac{1+z}{1-z}\right)^{m/2}
\mathbf{F}\left(-n, n+1; 1-m; \frac{1-z}{2}\right) \quad\text{(type 0)},
\qquad
\mathcal{P}_n^m(z) = \left(\frac{z+1}{z-1}\right)^{m/2}
\mathbf{F}\left(-n, n+1; 1-m; \frac{1-z}{2}\right) \quad\text{(type 1)},
$$
with the regularized ${}_2F_1$ (type 0 has cuts on $(-\infty, -1]$ and
$[1, \infty)$, type 1 on $(-\infty, 1]$), and the function of the
second kind is
$$
Q_n^m(z) = \frac{\pi}{2 \sin(\pi m)} \left( \cos(\pi m) P_n^m(z) -
\frac{\Gamma(1+m+n)}{\Gamma(1-m+n)} P_n^{-m}(z)\right)
$$
(and analogously for type 1). `acb_hypgeom_legendre_p` evaluates the
${}_2F_1$ representation with `acb_hypgeom_2f1`, which handles the
transformations needed for large $|z|$. For integers $0 \le m \le n$ it
uses instead the recurrence of `acb_hypgeom_legendre_p_uiui_rec` for
$(1-z^2)^{-m/2}P_n^m(z)$, starting from $P_m^m$ and $P_{m+1}^m$, when
the precision is low compared to $n$, $n - m$ is small, or $z$ is near
1 (where the transformations are singular); real $z$ with $m = 0$ goes
to `arb_hypgeom_legendre_p_ui`. If $n - m > 2p$ the recurrence would
lose all accuracy and the result is indeterminate.

For integer $m$, the defining formula for $Q_n^m$ is a $0 \cdot \infty$
limit, so `acb_hypgeom_legendre_q` uses connection formulas without the
removable singularity. For type 0 it evaluates two ordinary ${}_2F_1$
series in $z^2$ with gamma and trigonometric prefactors [@WQ3a]. For
type 1 it evaluates one regularized ${}_2F_1$ in $1/z^2$ [@WQ3b] when
the ball $z$ is away from the real half-line $[-1, \infty)$ and from the
square $|\operatorname{Re} z|, |\operatorname{Im} z| < 1/4$ (the
formula is invalid on $(-1, 0)$), and otherwise combines the type 0
function $Q_n^m$ and $\mathcal{P}_n^m$ [@WQ3c]. (The formula [@WQ3d] would also be convenient, but it
is incorrect in the lower half-plane.)

The spherical harmonics
$$
Y_n^m(\theta, \phi) = \sqrt{\frac{2n+1}{4\pi} \frac{(n-m)!}{(n+m)!}}
e^{im\phi} P_n^m(\cos\theta)
$$
are polynomials in $\cos\theta$ and $\sin\theta$, which
`acb_hypgeom_spherical_y` evaluates with
`acb_hypgeom_legendre_p_uiui_rec`.

## Other orthogonal polynomials {#orthogonal}

The other classical families are defined for complex degree through
${}_2F_1$ or ${}_1F_1$:

| Function | Definition | Direct method for integer $n \ge 0$ |
|---|---|---|
| Chebyshev $T_n(z)$, $U_n(z)$ | ${}_2F_1(-n, n; \frac12; \frac{1-z}{2})$, $(n+1)\,{}_2F_1(-n, n+2; \frac32; \frac{1-z}{2})$ | `acb_chebyshev_t_ui`, `acb_chebyshev_u_ui` (word-size $n$) |
| Jacobi $P_n^{(a,b)}(z)$ | $\frac{(a+1)_n}{n!}\,{}_2F_1(-n, n+a+b+1; a+1; \frac{1-z}{2})$ | explicit double sum for $n \le p$ and $a < -0.9$ |
| Gegenbauer $C_n^m(z)$ | $\frac{(2m)_n}{n!}\,{}_2F_1(-n, 2m+n; m+\frac12; \frac{1-z}{2})$ | recurrence in $n$ for $n \le p$ and $m < 0$ |
| Laguerre $L_n^m(z)$ | $\frac{(m+1)_n}{n!}\,{}_1F_1(-n; m+1; z)$ | recurrence in $n$ for $n \le p$ and $\lvert m\rvert > n$ |
| Hermite $H_n(z)$ | $2^n\sqrt\pi\left(\frac{{}_1F_1(-\frac n2; \frac12; z^2)}{\Gamma(\frac{1-n}{2})} - \frac{2z\,{}_1F_1(\frac{1-n}{2}; \frac32; z^2)}{\Gamma(-\frac n2)}\right)$ | — |

The direct methods cover the cases where the hypergeometric
representation is singular although the polynomial is defined (the
lower parameter $a + 1$, $m + \frac12$ or $m + 1$ is near a
nonpositive integer). For non-integer degree, the Chebyshev functions
use $T_\nu(\cos\theta) = \cos(\nu\theta)$ and $U_\nu(\cos\theta) =
\sin((\nu+1)\theta)/\sin\theta$ (with $\cosh$ for $z$ outside $[-1,1]$)
away from $z = 1$, and the ${}_2F_1$ series near it.

## Use in FLINT {#use}

- `arb_hypgeom_legendre_p_ui` and its methods
  `arb_hypgeom_legendre_p_ui_zero`, `arb_hypgeom_legendre_p_ui_one`,
  `arb_hypgeom_legendre_p_ui_asymp`, `arb_hypgeom_legendre_p_ui_rec`;
  `arb_hypgeom_legendre_p_ui_root` for quadrature, used by
  `acb_calc_integrate` and `arb_fpwrap`.
- `acb_hypgeom_legendre_p`, `acb_hypgeom_legendre_q`,
  `acb_hypgeom_legendre_p_uiui_rec`, `acb_hypgeom_spherical_y` and the
  real versions.
- `acb_hypgeom_chebyshev_t`, `acb_hypgeom_chebyshev_u`,
  `acb_hypgeom_jacobi_p`, `acb_hypgeom_gegenbauer_c`,
  `acb_hypgeom_laguerre_l`, `acb_hypgeom_hermite_h` and the real
  versions.
