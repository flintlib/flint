---
title: Legendre functions and orthogonal polynomials
summary: Legendre polynomials of large degree by the method of Johansson and Mezzarobba (hypergeometric series at 0 and 1, an asymptotic expansion with explicit error bounds, and a recurrence, chosen by cost and cancellation), Gauss–Legendre nodes and weights by Newton iteration, bit-burst Newton steps and Newton–Taylor steps with certified Taylor bounds, associated Legendre functions of both kinds through 2F1 and connection formulas, spherical harmonics, and the classical orthogonal polynomials and their extensions to non-integer degree.
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
  recurrence $(k+1)P_{k+1} = (2k+1)xP_k - kP_{k-1}$, scaled by its
  denominators: $u_k = k!\,P_k$ satisfies $u_{k+1} = (2k+1)x u_k - k^2
  u_{k-1}$ without division. The values are kept as limb vectors of
  fixed length with a common exponent (the low limb is dropped as they
  grow), each step is a high product `flint_mpn_mulhigh_n` by $x$ and
  two operations by words, and the result is divided by $n!$ at the
  end. The rounding errors are bounded a priori with the
  error-propagation estimate of [@JM2018, Proposition 5], which uses
  only $|P_k(x)| \le 1$. It is stable on $(-1, 1)$, costs $O(n)$
  short-product operations, and is the fastest method for moderate $n$
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
`arb_hypgeom_legendre_p_ui_root`. By the symmetry $x_{n-1-k} = -x_k$,
$w_{n-1-k} = w_k$ only the positive roots, $0 \le k < \lfloor n/2
\rfloor$, are computed; the middle root of odd $n$ is $0$. There are
three methods, chosen by $n$ and the precision $p$:

| condition | method |
|---|---|
| $n < 64$ or $p < 128$ | Newton iteration with evaluations of $P_n$ [@JM2018] |
| $0.3n \le p \le 15n$ | bit-burst Newton steps (`_arb_hypgeom_legendre_p_ui_root_bitburst`) |
| otherwise | Newton–Taylor steps (`_arb_hypgeom_legendre_p_ui_root_taylor`) |

The last two certify their result with an interval Newton step and an
isolating enclosure; if a certification step fails, the first method
is used.

### Initial approximation and isolation {#roots-isolation}

All methods start from the approximation of Petras [@Pet1999]: with
$\phi = \frac{4k+3}{4n+2}\pi$ (that is, $(4k'-1)\pi/(4n+2)$ for the
1-based index $k' = k + 1$),
$$
x_k \approx \cos\psi, \qquad \psi = \phi + \frac{2\cos\phi}{(4n+2)^2
\sin\phi}\left(1 - \frac{11}{(4n+2)^2 \sin^2\phi}\right),
$$
with an error of at most $\phi^2 \cdot 22 y^4 (1 + 2y^2)$, $y =
1/((4k'-1)\pi)$, a bound that is squared by each Newton step.

The bit-burst and Newton–Taylor methods identify the root through this
enclosure at about $64 + 2\log_2 n$ bits: it is checked to lie in $(0,
1)$ and to be disjoint from the enclosures of $x_{k-1}$ and $x_{k+1}$,
so that it contains $x_k$ and no other root. Their final interval
Newton step proves that a candidate interval contains exactly one root
of $P_n$, and containment of that interval in the Petras enclosure
proves that it is $x_k$.

### Newton iteration with full evaluations {#roots-newton}

1. Starting from Petras' approximation, Newton steps at a precision
   about 20% above the target square the a priori error bound until it
   is below $2^{-p}$.
2. When 40 bits plus a padding of $8 + 2\log_2 n$ bits is less than
   half of $p$, step 1 is carried out only at a low starting precision,
   and the enclosure is refined by interval Newton steps $x \mapsto m -
   P_n(m)/P_n'(X)$ at doubling precisions, where $m$ is the midpoint of
   the current enclosure $X$ and $P_n'(X)$ is bounded by $P_n'(m)$ plus
   a bound for $P_n''$ from `arb_hypgeom_legendre_p_ui_deriv_bound`
   times the radius.
3. The weight comes from $P_n'$ at the enclosure of $x_k$.

All evaluations of $P_n$ and $P_n'$ use `arb_hypgeom_legendre_p_ui`;
the last one, at full precision and with a full-precision argument,
dominates the cost.

### Taylor bounds {#roots-taylor-bounds}

The other two methods use the Taylor series $P_n(m + t) = \sum_j c_j
t^j$ at a real point $m$. For $z$ in the disk $|z - m| \le r$, the
Laplace integral
$$
P_n(z) = \frac{1}{\pi} \int_0^\pi \left(z + \sqrt{z^2-1}\cos\varphi\right)^n d\varphi
$$
gives $|P_n(z)| \le M = \max \rho(z)^n$ with $\rho(z) = \max|z \pm
i\sqrt{1-z^2}|$, since the integrand is a convex combination of the
two endpoints; $M$ is bounded on a box enclosing the disk, which must
lie in the half-plane $\operatorname{Re}(1 - z^2) > 0$. The Cauchy
estimates $|c_j| \le M/r^j$ give, for $|t| \le \tau = q r$ with $q <
1/2$,
$$
\sum_{j > R} |c_j| \tau^j \le \frac{M q^{R+1}}{1-q}, \qquad
\sum_{j > R} j |c_j| \tau^{j-1} \le \frac{M}{r}\,\frac{(R+1) q^R}{(1-q)^2}.
$$
The radius is $r = 0.99\min(\sqrt{1-m^2}/n,\ (1-m)/4)$, about the
distance over which $P_n$ changes by a bounded factor: about $1/n$ in
the interior and $(k+1)^2/n^2$ near $\pm 1$.

### Bit-burst Newton steps {#roots-bitburst}

For precision comparable to $n$, the root is computed by the bit-burst
idea for holonomic functions [@CC1990; @vdH1999]: every evaluation
point has only about twice as many bits as the accuracy of the
approximation it represents, so that the full-precision work multiplies
full-precision numbers only by short integers.

1. **Starting point.** $x^{(0)} = A/2^B$, from Petras' approximation
   and three Newton steps in double precision. The first transfer needs
   $|x^{(0)} - x_k|$ well below the radius $r$ of the Taylor bounds.
   With $B = w - 2 - \lceil \log_2 2n \rceil$ for the word size $w$,
   the multipliers of stage 0 fit in a word; this is used when $2B \ge
   w$ and $2^{-B} \le r/64$. Otherwise $B$ is the larger of $62 -
   \lceil \log_2 2n \rceil$ and $\log_2(1/r) + 8$, and for $r <
   2^{-44}$, beyond the accuracy of double precision, $x^{(0)}$ is the
   midpoint of Petras' enclosure, whose error is a small fraction of
   $r$.
2. **Stage 0** computes $P_n$ and $P_{n-1}$ at $x^{(0)}$ to the working
   precision. With $u_j = j!\,2^{Bj} P_j(x^{(0)})$ the three-term
   recurrence becomes
   $$
   u_{j+1} = (2j+1) A\, u_j - j^2 2^{2B} u_{j-1},
   $$
   with no division. With the word-size choice of $B$ and below 130000
   bits, the recurrence runs on two's complement limb windows of fixed
   size (one `mpn_mul_1` and one `mpn_submul_1` per step), dropping the
   low limbs as the values grow; otherwise exact products of blocks of
   steps, each of about the working precision and formed by binary
   splitting, are applied in turn to a ball vector.
3. **Transfer stages.** The Newton correction $-P_n/P_n'$ is rounded
   to $\delta = D/2^E$ with about twice as many bits as the current
   accuracy, and $(P_n, P_n')$ are moved from $x = A/2^B$ to $x +
   \delta$ by the Taylor series. From the Legendre equation
   $(1-x^2)y'' - 2xy' + n(n+1)y = 0$, the terms $a_j = c_j \delta^j$
   satisfy
   $$
   \gamma_j a_{j+2} = \alpha_j a_{j+1} + \beta_j a_j, \quad
   \alpha_j = 2AD(j+1)^2 2^{E-B}, \quad
   \beta_j = D^2 (j(j+1) - n(n+1)), \quad
   \gamma_j = \sigma (j+1)(j+2) 2^{2(E-B)},
   $$
   with $\sigma = 2^{2B} - A^2$. The sums $\sum a_j$ and $\sum j a_j$,
   which give $P_n(x + \delta)$ and $\delta P_n'(x+\delta)$, are
   computed in fixed point, with the terms as integers whose lengths
   decrease with $j$, so that the precision tapers by itself; the
   powers of two are split off from the divisions. When $K \cdot p$
   exceeds $10^{10}$, exact blocks of steps formed by binary splitting
   are applied in ball arithmetic instead. The number of terms $K$ is
   chosen so that the terms beyond it are below the working precision,
   and their sum is bounded by the Taylor bound above. Each stage
   doubles the size of the integers and halves the number of terms.
4. **Final step.** At the last point, within about $2^{-p/2}$ of $x_k$,
   an interval Newton step on the Taylor polynomial of degree 2 plus
   its tail bound encloses $x_k$, and $P_n'(x_k)$ for the weight comes
   from the same expansion.

The products of $2 \times 2$ blocks use `fmpz_mat_mul_fft_small` for
entries above 10000 bits.

**Error bound for the limb recurrence.** When a step drops limbs, the
new $u_{j+1}$ and the new $u_j$ lose less than one unit $2^{wS}$ of
the lowest kept limb each ($S$ the limb exponent of the window); in
terms of $P_i = u_i / (i!\,2^{Bi})$, these are errors in $P_{j+1}$ and,
through the multiplier, in $P_{j+2}$, of at most $(j+2)\, 2^{B + wS} /
((j+1)!\,2^{B(j+1)})$ together. The computed sequence is the exact one
plus, for each error $e$ in $P_i$, the solution $y$ of the recurrence
with $y_i = e$, $y_{i-1} = 0$. With $x = \cos\theta$, $s = \sin\theta$
and
$$
E_m = m\,(y_m^2 - 2x y_m y_{m-1} + y_{m-1}^2) = m (y_{m-1} - x y_m)^2 + m s^2 y_m^2,
$$
the recurrence $(m+1) y_{m+1} = (2m+1) x y_m - m y_{m-1}$ gives
$$
E_{m+1} = \frac{m^2}{m+1} (y_{m-1} - x y_m)^2 + (m+1) s^2 y_m^2 \le \frac{m+1}{m} E_m,
$$
so $E_m \le (m/i) E_i = m e^2$ and $s^2 y_m^2 \le E_m/m \le e^2$:
every error reaches $P_n$ and $P_{n-1}$ multiplied by at most $1/s$.
The bound is accumulated during the recurrence in `mag_t` arithmetic.

**Error bound for the fixed-point transfers.** The products in the
recurrence for $a_j$ are exact and each division by $\gamma_j$
truncates, adding less than one unit. With $\epsilon_j$ the error of
$a_j$ in units,
$$
\epsilon_{j+2} \le \frac{|\alpha_j| \epsilon_{j+1} + |\beta_j| \epsilon_j}{\gamma_j} + 1,
$$
a contracting recurrence since $|\delta|$ is small compared with the
radius of convergence; it is also evaluated in `mag_t` arithmetic,
starting from the radii of $P_n$ and $\delta P_n'$.

### Newton–Taylor steps {#roots-taylor}

For precision much larger or much smaller than $n$:

1. An enclosure $m \pm \tau$ of $x_k$ at a quarter of the precision,
   computed recursively by the method chosen for that precision.
2. One evaluation of $P_n(m)$ and $P_n'(m)$ at the full precision with
   `arb_hypgeom_legendre_p_ui`, and the Taylor coefficients $c_j$ from
   the recurrence
   $$
   (1-m^2)(j+1)(j+2) c_{j+2} = 2m(j+1)^2 c_{j+1} + (j(j+1) - n(n+1)) c_j,
   $$
   with the precision of $c_j$ tapered by $\tau^j$. The degree $R$ is
   the smallest for which the Taylor bound at $|t| \le \tau$ is below
   the working precision.
3. The fixed-point iteration $d \leftarrow -(c_0 + \sum_{j \ge 2} c_j
   d^j)/c_1$ on midpoints, which gains about a quarter of the precision
   per step, and a final interval Newton step on the Taylor polynomial
   of degree $R$ plus its tail bound.

The single full-precision evaluation of $P_n$ dominates the cost when
the precision is much larger than $n$; when it is much smaller, that
evaluation uses the asymptotic expansion.

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
  `arb_hypgeom_legendre_p_ui_root` and its methods
  `_arb_hypgeom_legendre_p_ui_root_bitburst` and
  `_arb_hypgeom_legendre_p_ui_root_taylor` for quadrature, used by
  `acb_calc_integrate` and `arb_fpwrap`.
- `acb_hypgeom_legendre_p`, `acb_hypgeom_legendre_q`,
  `acb_hypgeom_legendre_p_uiui_rec`, `acb_hypgeom_spherical_y` and the
  real versions.
- `acb_hypgeom_chebyshev_t`, `acb_hypgeom_chebyshev_u`,
  `acb_hypgeom_jacobi_p`, `acb_hypgeom_gegenbauer_c`,
  `acb_hypgeom_laguerre_l`, `acb_hypgeom_hermite_h` and the real
  versions.
