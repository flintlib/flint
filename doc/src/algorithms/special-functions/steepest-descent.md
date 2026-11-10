---
title: Bessel and incomplete gamma functions by steepest-descent integration
summary: Rigorous evaluation of the Bessel functions J, Y, I, K, the Hankel functions and the incomplete gamma functions for arbitrary complex parameters by Gauss–Legendre quadrature along steepest-descent contours planned in double precision, with L1 bounds along the same contours for the functions and their derivatives; the treatment of Stokes lines, infinite ladders of saddles and branch cuts; precision tapering; numerically generated saddle-point (Debye) expansions with rigorous error bounds for large parameters; and the automatic choice between integration, convergent series and asymptotic expansions.
order: 47
---

The convergent and asymptotic series of
[](algorithms/special-functions/bessel-airy) and
[](algorithms/special-functions/error-functions) leave a gap: when the
order and the argument are both large and of comparable size, the
convergent series needs many terms and loses many digits to
cancellation, and the asymptotic expansion does not reach the target
precision. Both also behave poorly for wide input balls. FLINT fills
the gap by integrating a contour integral representation numerically
along a contour through the saddle points of the integrand, where the
integrand does not oscillate. The cost is nearly independent of the
parameters: a few hundred evaluations of the integrand at machine
precision. The same contours give bounds for the functions and their
derivatives that are tight in all parameter regimes.

The method has three layers: a contour is *planned* heuristically in
double precision; *bounds* for the integral of the absolute value of
the integrand along the planned contour are computed rigorously, still
in double precision; and the integral is *evaluated* by Gauss–Legendre
quadrature with certified error bounds, in multi-double or
arbitrary-precision ball arithmetic. Since the integrand is entire,
correctness never depends on the quality of the contour, only the
accuracy and the speed. When the parameters are large, the quadrature
can be replaced by a *saddle-point expansion* of the integral along
the same contour ([](#expansions)), which is cheaper, especially at
high precision.

## Integral representations {#representations}

All the functions are integrals of $e^{\varphi(w)}$ with
$$
\varphi(w) = \alpha e^{w} + \beta e^{-w} - \nu w,
$$
between two *ends*: valleys at infinity where $\operatorname{Re}
\varphi \to -\infty$, or a finite starting point. Write $w = x + iy$.
There are three kinds of valleys:

- $V^{+}(k)$: $x \to +\infty$ with $y$ near $\pi - \arg \alpha + 2 \pi
  k$, where $\operatorname{Re}(\alpha e^w) \to -\infty$;
- $V^{-}(k)$: $x \to -\infty$ with $y$ near $\arg(-\beta) + 2 \pi k$
  (if $\beta \ne 0$);
- the *linear valley*: $w \to \infty$ in a direction $e^{i \theta}$
  along which $\operatorname{Re}(\nu w) \to +\infty$ while $\alpha e^w +
  \beta e^{-w}$ stays bounded. This requires $\cos \theta = 0$, or
  $\cos \theta < 0$ if $\beta = 0$.

**Bessel functions.** With $\alpha = -\beta = \zeta / 2$, that is,
$\varphi(w) = \zeta \sinh w - \nu w$, the Schläfli–Sommerfeld integrals
[@NIST2012, §10.9] read
$$
J_\nu(z) = \frac{1}{2 \pi i} \int_{V^{+}(-1)}^{V^{+}(0)} e^{\varphi(w)} dw,
\quad
H^{(1)}_\nu(z) = \frac{1}{\pi i} \int_{V^{-}(0)}^{V^{+}(0)} e^{\varphi(w)} dw,
\quad
H^{(2)}_\nu(z) = -\frac{1}{\pi i} \int_{V^{-}(0)}^{V^{+}(-1)} e^{\varphi(w)} dw
$$
with $\zeta = z$, and $Y_\nu = (H^{(1)}_\nu - H^{(2)}_\nu) / (2i)$. The
modified functions are
$$
I_\nu(z) = e^{\mp i \pi \nu / 2} J_\nu(\pm i z), \qquad
K_\nu(z) = \tfrac{1}{2} e^{-i \pi \nu / 2}
    \int_{V^{-}(0)}^{V^{+}(-1)} e^{\varphi(w)} dw \quad (\zeta = -iz),
$$
or $K_\nu(z) = \tfrac{1}{2} e^{i \pi \nu / 2} \int_{V^{-}(0)}^{V^{+}(0)}$
with $\zeta = i z$; the rotation is chosen so that $\arg \zeta$ stays
away from $\pm \pi$. The valleys are defined with $\arg \zeta \in (-\pi,
\pi]$, so that the formulas give the principal branches, continuous
from above on the cut $(-\infty, 0)$.

**Incomplete gamma functions.** In the variable $w = \log t$, with
$\alpha = -1$, $\beta = 0$, $\nu = -s$, that is, $\varphi(w) = -e^w + s
w$,
$$
\Gamma(s, z) = \int_{\operatorname{Log} z}^{V^{+}(0)} e^{\varphi(w)} dw,
\qquad
\gamma(s, z) = \int_{V}^{\operatorname{Log} z} e^{\varphi(w)} dw,
$$
where $V$ is the linear valley, to the left. For $\operatorname{Re} s >
0$ this is the defining integral of $\gamma$; if $\operatorname{Re} s
\le 0$ and $\operatorname{Im} s \ne 0$, the valley has turned away from
the real axis and the integral is the analytic continuation in $s$.

**Symmetries.** For real $\nu$ and $z > 0$ the two Hankel integrals are
complex conjugates, so that $J_\nu$ and $Y_\nu$ are the imaginary and
real parts of a single integral. $J_\nu$ is obtained this way only for
$z \ge \nu$; for $z < \nu$ it is exponentially smaller than $Y_\nu$ and
is integrated directly. The reflection formulas $K_{-\nu} = K_\nu$,
$H^{(1)}_{-\nu} = e^{i \pi \nu} H^{(1)}_\nu$, $H^{(2)}_{-\nu} = e^{-i
\pi \nu} H^{(2)}_\nu$ are applied when $\operatorname{Re} \nu < 0$, as
are $J_{-n} = (-1)^n J_n$, $Y_{-n} = (-1)^n Y_n$ and $I_{-n} = I_n$ at
integers, where the integrals for $J$ and $I$ would cancel.

**Branch cuts.** The integrals are analytic in the parameters, so for
an input ball the computed enclosures and bounds hold for the analytic
continuation from the midpoint. A ball $z$ which crosses the cut
$(-\infty, 0)$ is therefore split into the parts with $\operatorname{Im}
z \ge 0$ and $\operatorname{Im} z \le 0$, and the union of the two
enclosures is returned.

## Planning the contour {#planning}

**Saddles.** With $u = e^w$, $\varphi'(w) = 0$ is the quadratic
equation $\alpha u^2 - \nu u - \beta = 0$, so the saddles form (at
most) two families $\log u_{1,2} + 2 \pi i k$. For the Bessel functions
they are $\pm \operatorname{acosh}(\nu / \zeta) + 2 \pi i k$, and for
the incomplete gamma functions $\log s + 2 \pi i k$. The roots are
computed without cancellation from $u_1 = (\nu + \sqrt{\nu^2 + 4
\alpha \beta}) / (2 \alpha)$ (with the sign of the square root that
makes the sum largest) and $u_1 u_2 = -\beta / \alpha$. Two saddles
closer than the width $|\varphi'''|^{-1/3}$ of the cubic saddle into
which they merge when $\nu^2 + 4 \alpha \beta = 0$ (the transition
region $\nu \approx \pm \zeta$) are replaced by a single degenerate
saddle with three descent directions; the first step along each of
them is long enough (three times the distance between the two saddles)
to leave the region where the gradient is that of the pair and not of
the cubic. A finite starting point is treated as a saddle with a
single descent direction.

**Descent paths.** From each saddle $w_s$, the paths of steepest
descent of $\operatorname{Re} \varphi$ leave in the directions $\pm
\exp(i (\pi - \arg \varphi''(w_s)) / 2)$. They are followed with the
explicit midpoint rule for $w' = -\overline{\varphi'(w)} /
|\varphi'(w)|$, with a step that resolves the local scale of the
integrand: at most $|\varphi''|^{-1/2}$, at most a fixed decrease of
$\operatorname{Re} \varphi$, and small enough that the cubic term in
the bounds below is $O(1)$. A path ends

- in a *valley*, when the integrand has decreased by the target factor
  and the tail bound along a ray into the valley
  ([](#tails)) is negligible. Deep in a valley several rays can
  qualify; the path is assigned to the valley whose ray is best aligned
  with its direction;
- in a *lower saddle*, when the path passes within the width of that
  saddle. This is the situation on a Stokes line, where one saddle
  connects to another; near a Stokes line it avoids a contour that
  nearly returns along itself.

**Routes.** The saddles and the valleys are the nodes of a graph, and
the descent paths its edges. Any path in the graph between the two
ends is a valid contour (by Cauchy's theorem, since the integrand is
entire and decays in the valleys), and its cost is the sum of the
$L^1$ bounds of its edges. The cheapest route is found with the
Bellman–Ford algorithm. Saddles are traced lazily in order of
increasing $\operatorname{Re} \varphi(w_s)$, among those between the
two ends, until a route is found which the remaining saddles cannot
improve. Two consecutive descent paths from the same saddle are joined
into one polygon passing through the saddle.

### Ladders of saddles {#ladders}

When $\operatorname{Im} \nu \ne 0$, the Bessel-type problem has a
linear valley in the vertical direction $\sigma i$, $\sigma =
-\operatorname{sgn} \operatorname{Im} \nu$. The valleys $V^{\pm}(k + j
\sigma)$, $j \to \infty$, lead to it, separated from each other by
periodic copies of the same saddles, whose contributions decrease
geometrically. When $|\operatorname{Im} \nu|$ is large compared to
$|\zeta|$, a contour between two given valleys can then be forced to
pass through the linear valley, and thus over infinitely many
saddles.

::: theorem Ladders {#thm-ladder}
Let $\alpha, \beta \ne 0$, $\sigma = \pm 1$ with $\sigma
\operatorname{Im} \nu < 0$, and $q = e^{-2 \pi i \sigma \nu}$, so that
$|q| < 1$. Then for $V = V^{+}$ or $V^{-}$ and every $k$,
$$
\int_{V(k)}^{\text{linear valley}} e^{\varphi(w)} dw =
\frac{1}{1 - q} \int_{V(k)}^{V(k + \sigma)} e^{\varphi(w)} dw.
$$
:::

::: proof
Since $\varphi(w + 2 \pi i \sigma) = \varphi(w) - 2 \pi i \sigma \nu$,
translating a contour from $V(k)$ to $V(k + \sigma)$ by $2 \pi i
\sigma j$ gives $\int_{V(k + j \sigma)}^{V(k + (j+1) \sigma)} e^{\varphi}
= q^j \int_{V(k)}^{V(k + \sigma)} e^{\varphi}$. The sum over $j \ge 0$
is a geometric series, and the concatenated contours can be deformed
into a single contour from $V(k)$ into the linear valley: on the
segments joining $V(k + j \sigma)$ at $x = \pm X$ to the line $x = 0$,
$\operatorname{Re} \varphi \le (|\alpha| + |\beta|) e^{X} + X
|\operatorname{Re} \nu| - |\operatorname{Im} \nu| |y|$, which tends to
$-\infty$ as $j \to \infty$ for fixed $X$.
:::

If there is no direct route between two valleys $a$ and $b$, each of
them is therefore connected to the linear valley, by a direct route if
there is one and otherwise by [](#thm-ladder), and $\int_a^b = \int_a^{\text{lin}}
- \int_b^{\text{lin}}$. The pieces coming from the theorem are divided
by $1 - q$ after the quadrature. This is the contour-integral form of
connection formulas such as $K_\nu = \frac{\pi}{2} (I_{-\nu} - I_\nu) /
\sin(\pi \nu)$.

For the incomplete gamma functions the same structure appears when
$|\operatorname{Im} s|$ is large: the valley $V$ is then separated from
$V^{+}(0)$ by infinitely many saddles, whose contributions sum to
$\Gamma(s)$. Here the planner looks for a route for the requested
function and, failing that, for the complementary one; the requested
function then follows from $\Gamma(s, z) + \gamma(s, z) = \Gamma(s)$,
normally without cancellation, since the descent path from
$\operatorname{Log} z$ leads to the valley of the smaller of the two.

### Fixed polygons {#polygons}

Where the landscape is flat (small parameters), the descent flow does
not lead anywhere within the range of double precision, and the
classical contours are used instead: for the Bessel functions, rays
along the centres of the two valleys, starting where $|\alpha| e^{x}
\ge 1$ or $|\beta| e^{-x} \ge 1$ and joined by straight segments; for
$\Gamma(s, z)$, the vertical segment from $\operatorname{Log} z$ to the
real axis followed by the real axis. These are the parameters where
the convergent series work best, so this is only a fallback.

## Bounds {#bounds}

::: theorem {#thm-l1}
Let $C$ be any contour between the two ends, let $g \ge 0$, and let
$\alpha, \beta, \nu$ range over balls with midpoints $\alpha_0,
\beta_0, \nu_0$ and radii $r_\alpha, r_\beta, r_\nu$. Then for all
parameters in the balls
$$
\left| \int_C g(w)\, e^{\varphi(w)} dw \right| \le
\int_C g(w) \exp\left(\operatorname{Re} \varphi_0(w) + r_\alpha e^{x} +
r_\beta e^{-x} + r_\nu |w|\right) |dw|,
$$
provided the right-hand side converges.
:::

::: proof
$\varphi$ is linear in the parameters, so $\operatorname{Re}(\varphi -
\varphi_0) \le |\alpha - \alpha_0| e^{x} + |\beta - \beta_0| e^{-x} +
|\nu - \nu_0| |w|$. With this uniform majorant the integral converges
uniformly on the balls, hence is analytic there and equal to the
continuation of the function from the midpoint.
:::

With $g = 1$ this bounds the function itself. Differentiating under
the integral sign, $\partial_\zeta e^{\varphi} = \sinh(w) e^{\varphi}$
and $\partial_\nu e^{\varphi} = -w e^{\varphi}$, so the weights $g =
(e^{x} + e^{-x}) / 2 \ge |\sinh w|$ and $g = |w|$ give bounds for the
derivatives with respect to the argument and the order of the Bessel
functions, and $g = |w|$ for the derivative of the incomplete gamma
functions with respect to $s$. Along a path of steepest descent the
integrand has constant phase, so for a contour through a single saddle
the bound with $g = 1$ equals the integral; the bounds are tight up to
the cancellation between different saddles, that is, they bound the
envelope of an oscillating function.

**Panels.** The contour is a polygon. On a panel $[w_a, w_b]$ with
midpoint $m$ and $h = (w_b - w_a) / 2$, Taylor's theorem gives for $t
\in [-1, 1]$
$$
\operatorname{Re} \varphi(m + h t) \le \operatorname{Re} \varphi(m) +
\operatorname{Re}(\varphi'(m) h)\, t + \tfrac{1}{2}
\operatorname{Re}(\varphi''(m) h^2)\, t^2 + \tfrac{1}{6} |h|^3 S,
\qquad S = \sup |\varphi'''| \le |\alpha| e^{x_{\max}} + |\beta|
e^{-x_{\min}},
$$
and the integral of the exponential of a quadratic $A + B t + C t^2$
over $[-1, 1]$ is at most $2 e^{M}$, $M$ its maximum, and at most
$\sqrt{\pi / (-C)}\, e^{A - B^2 / (4C)}$ if $C < 0$. The weights and
the terms for the input balls are bounded by their suprema on the
panel.

### Tails {#tails}

Beyond the last vertex $w_e = x_e + i y_e$ of a polygon, the contour
continues along a ray.

- *Into $V^{+}(k)$*, horizontally. With $u = e^{x - x_e} \ge 1$, and
  $c = -\cos(\arg \alpha + y_e) > 0$, $\operatorname{Re}(\alpha e^w) =
  -|\alpha| c\, e^{x_e} u$ and $\operatorname{Re}(\beta e^{-w}) =
  \operatorname{Re}(\beta e^{-w_e}) / u$, so using $1 - 1/u \le \log u$
  and $x - x_e = \log u$, the exponent of [](#thm-l1) is at most $G - K
  (u - 1) + a \log u$, where $G$ is its value at $w_e$, $K = (|\alpha|
  c - r_\alpha) e^{x_e}$ and $a = \max(-\operatorname{Re}(\beta
  e^{-w_e}), 0) + \max(-\operatorname{Re} \nu, 0) + r_\nu$. If $K > 0$,
  the tail is at most $e^{G + K} \int_1^\infty e^{-K u} u^{a - 1} du$, an
  incomplete gamma function which is bounded in closed form. The valley
  index $k$ is read off from $y_e$.
- *Into $V^{-}(k)$*: the same, by the symmetry $\varphi(-w; \alpha,
  \beta, \nu) = \varphi(w; \beta, \alpha, -\nu)$.
- *Into the linear valley*, along $w_e + t e^{i \theta}$ with $\mu =
  \operatorname{Re}(\nu e^{i \theta}) - r_\nu > 0$. Since $\cos \theta
  \le 0$, $|\alpha e^w|$ does not increase (and $|\beta e^{-w}|$ is
  constant if $\beta \ne 0$, where $\cos \theta = 0$), so the exponent
  is at most $G' - \mu t$ with $G' = (|\alpha| + r_\alpha) e^{x_e} +
  (|\beta| + r_\beta) e^{-x_e} - \operatorname{Re}(\nu w_e) + r_\nu
  |w_e|$, and the tail is at most $e^{G'} / \mu$. The direction is
  $e^{i \theta} = \overline{\nu} / |\nu|$ if $\beta = 0$ and
  $\operatorname{Re} \nu \le 0$, and vertical otherwise.

The weights $g$ are carried along: they are at most $g(w_e) u$ on a
horizontal ray, and $|w| \le |w_e| + t$ on any ray.

**Rounding.** All of this is evaluated in double precision. Only
$O(1)$ absolute accuracy of $\operatorname{Re} \varphi$ matters, and
an explicit margin $16 \epsilon (|\alpha| e^{x} + |\beta| e^{-x} +
|\nu| |w| + 1)$ is added to every exponent, with the radii of the
input balls enlarged to cover the rounding of the midpoints to
doubles. This requires $|\alpha|, |\nu| < 2^{44}$. The bounds cost a
few hundred evaluations of `exp`, `sin` and `cos`: some tens of
microseconds.

## Quadrature {#quadrature}

Each piece of the contour (a polygon with tails) is integrated as
$e^{c} \int e^{\varphi(w) - c} dw$, where $c = \varphi(w_s)$ at the
saddle of the piece, rounded to a complex double; $|e^{\varphi - c}|
\le 1$ along the piece up to rounding, and the pieces are combined,
with the factors $e^{c}$, in `acb` arithmetic. The tails are added to
the radius.

::: theorem {#thm-gl}
Let $f$ be analytic with $|f| \le M$ on the closed disk $|w - m| \le
|h| X$, $X > 1$, and let $\rho = X + \sqrt{X^2 - 1}$. The error of the
$n$-point Gauss–Legendre rule for $\int_{m - h}^{m + h} f(w) dw$ is at
most
$$
\frac{64}{15} \frac{M |h|}{(\rho - 1) \rho^{2n - 1}}.
$$
:::

This is the standard bound in terms of the Bernstein ellipse with
parameter $\rho$ [@Tre2008; @Joh2018a], which is contained in the
disk. For $f = e^{\varphi - c}$, $M$ follows from an exact identity:
since $\alpha e^{m + u} + \beta e^{-m - u} = \varphi''(m) \cosh u +
\varphi'''(m) \sinh u$ (with $\varphi''' = \varphi' + \nu$),
$$
\varphi(m + u) - \varphi(m) = \varphi'(m) u + \varphi''(m) (\cosh u -
1) + \varphi'''(m) (\sinh u - u),
$$
so that for $|u| \le R = |h| X$,
$$
\operatorname{Re}(\varphi(m + u) - \varphi(m)) \le |\varphi'(m)| R +
|\varphi''(m)| (\cosh R - 1) + |\varphi'''(m)| (\sinh R - R),
$$
to which $r_\alpha e^{x + R} + r_\beta e^{-x + R} + r_\nu (|m| + R)$
is added for the input balls ([](#thm-l1)).

**Chords.** The vertices of the planned polygon are doubles, and the
rule is applied on *chords* joining two of them. For a chord, the
number of nodes required for the error to be below the tolerance is
minimized over $X$; starting from each vertex, the chord is extended
over the following panels of the polygon for as long as this lowers
its number of nodes per panel (provided the vertices skipped stay
close to the chord). The rules are those with $n = 1, 2, 4, 6, 8, 12,
16, 22, 32, \ldots$ nodes, in steps of about $2^{1/2}$ up to about
$p/2$ nodes at $p$ bits, from the
cache of `acb_calc_integrate` (which is shared; a conversion to
double-based balls is cached in addition). Panels on
which the bound of [](#bounds) is below the tolerance are not
integrated; their bound is added to the radius. A panel on which no
rule is certified is bisected.

**Starting point.** For the incomplete gamma functions, the contour
starts at the double approximation of the midpoint of
$\operatorname{Log} z$, and one more chord joins it to
$\operatorname{Log} z$ itself, a ball. The Gauss–Legendre sum with a
ball endpoint encloses the integrals over all the chords from points
of the ball, and the error bound is taken on a disk containing all
their Bernstein ellipses. Since the chord is very short (for an exact
$z$, of length $2^{-53} |\operatorname{Log} z|$), $X$ can be taken
very large, and one or two nodes suffice.

**Arithmetic.** The nodes are evaluated in ball arithmetic, with the
input balls themselves as parameters, so that the result encloses the
function on the balls. Computing $\varphi(w) - c$ cancels about
$\log_2 (|\alpha| e^{x} + |\beta| e^{-x} + |\nu| |w|)$ bits, which are
added to the working precision. There are two implementations of the
rule on a chord:

- up to about 200 bits, in the vectorized multi-double ball arithmetic
  of [`dfloat`](dfloat) with $N \le 4$ doubles: the exponential and
  trigonometric functions of all the nodes of a chord are computed
  with a few vector operations, and complex numbers are kept as pairs
  of real vectors;
- at higher precision, or if `dfloat` is not supported, in `arb`
  arithmetic, with the nodes of a chord evaluated in parallel when
  several threads are available (`flint_set_num_threads`).

**Tapering.** A node where $|e^{\varphi - c}| = 2^{-k}$ only needs to
be evaluated with $k$ bits less for the same absolute error. With
`arb` arithmetic the precision is chosen per node from a double
precision estimate of $\operatorname{Re}(\varphi - c)$; with `dfloat`,
$N$ is chosen per chord. Since the integrand decays like a Gaussian
along the contour, this saves about a third of the work at high
precision.

**Input balls.** The radius of the result due to the input balls is
about $r_\alpha L_1 + r_\nu L_2$ in terms of the weighted $L^1$ norms
of [](#bounds), which the planner has computed. The tolerance of the
quadrature and the working precision are relaxed accordingly, so that
wide balls are cheaper, not more expensive, than exact inputs.

## Saddle-point expansions {#expansions}

For large parameters, the integrand is concentrated in a small
neighbourhood of each saddle, and the integral over the two descent
paths of a simple saddle $w_s$ can be expanded there. This is the
expansion of Debye type; for the Bessel functions it reproduces
Debye's expansions, but the coefficients are generated numerically for
the given parameters instead of as polynomials, so that no table is
needed, complex parameters need no special treatment, and the error
bound holds for the parameters at hand.

Let $c = \varphi(w_s)$ and $A = \varphi''(w_s)$. By the identity of
[](#quadrature) at $m = w_s$,
$$
\varphi(w_s + u) - c = P(u) = A (\cosh u - 1) + \nu (\sinh u - u).
$$
On the two descent paths $P(u) = -\tau^2$ with $\tau$ real ($\tau > 0$
on the outgoing path), and $u(\tau) = \sum_{k \ge 1} b_k \tau^k$ is
analytic near $\tau = 0$, with $b_1 = \sqrt{-2/A}$. The integral over
the piece is
$$
e^{c} \int e^{-\tau^2} u'(\tau) \, d\tau \sim e^{c} \sum_{k
\text{ odd}} k \, b_k \, \Gamma(k/2).
$$

**Coefficients.** With $D = P'(u) = A \sinh u + \nu (\cosh u - 1)$,
the series satisfies $D u' = -2 \tau$. The power series of $\sinh u$
and $\cosh u$ follow from $(\sinh u)' = u' \cosh u$, $(\cosh u)' =
u' \sinh u$, and comparing the coefficients of $\tau^k$ in $D u'$
gives $b_k$ from $b_1, \ldots, b_{k-1}$: three convolutions per
coefficient, $O(n^2)$ operations for $n$ coefficients, done with
`acb_dot`. Since the terms decrease by a factor of about $2^{-p}$ over
the $n$ coefficients used, the working precision is lowered linearly
along the recurrence.

::: theorem {#thm-saddle}
For $r > 0$ let
$$
G(r) = \frac{2 (\cosh r - 1 - r^2/2)}{r^2} + \left| \frac{\nu}{A}
\right| \frac{2 (\sinh r - r)}{r^2}.
$$
If $G(r) < 1$, then $u(\tau)$ is analytic with $|u(\tau)| < r$ for
$|\tau| < \rho = \sqrt{|A|/2} \; r \sqrt{1 - G(r)}$, and $|b_k| \le r
\rho^{-k}$.
:::

Indeed, $P(u) = (A/2) u^2 (1 + g(u))$ with $|g(u)| \le G(r)$ for $|u|
\le r$, so $s(u) = u \sqrt{1 + g(u)}$ satisfies $|s(u) - u| \le r (1 -
\sqrt{1 - G})$ on $|u| = r$, and by Rouché's theorem takes every value
of modulus less than $r \sqrt{1 - G}$ exactly once in $|u| < r$; $\tau
= \sqrt{-A/2} \, s(u)$, and the bound for $b_k$ is Cauchy's
inequality. The radius $r$ is chosen to maximize $\rho$.

**Error bound.** Let $0 < \tau_0 < \rho$ and $t = \tau_0 / \rho$. On
each of the two sides, the difference between the integral over the
descent path and the $n$ terms consists of three parts:

- the truncation of the series on $[0, \tau_0]$, at most $r \sum_{k >
  n} \min(k \Gamma(k/2) / 2, \tau_0^k) \rho^{-k}$;
- the terms integrated beyond $\tau_0$, at most $r e^{-\tau_0^2}
  \tau_0^{-2} t / (1 - t)^2$ when $\tau_0^2 \ge n - 2$;
- the integral over the rest of the path. The point $w_s + u(\tau_0)$
  is enclosed by the Taylor polynomial with the remainder $r t^{n+1} /
  (1 - t)$, and since the integrand is entire the rest of the path can
  be replaced by the segment from this point to a vertex of the
  planned polygon followed by the rest of the polygon and its tail,
  which are bounded as in [](#bounds).

$\tau_0$ is taken at the first vertex of the polygon where the
integrand is below the tolerance. The first bound decreases like the
terms of the expansion until $k \approx 2 \rho^2$, so the expansion
reaches a precision of about $\rho^2 / \log 2$ bits, and needs fewer
terms the larger $\rho^2 \approx |A|$ is: for $\nu = 10^{10}$, $z =
1.3 \cdot 10^{10}$, 5 coefficients for 53 bits, 24 for 333 bits and
264 for 3333 bits.

**Use.** The expansion is implemented for the Bessel kernel with
exact parameters, where $A = \pm\sqrt{\nu^2 - 4 \alpha^2}$, $e^{w_s} =
(\nu + A) / (2 \alpha)$ and $c = A - \nu w_s$ are computed directly
in ball arithmetic, and for contours consisting of pairs of descent
paths from one valley to another through simple saddles (not for the
half paths that meet at a second saddle on a Stokes line, nor for the
incomplete gamma functions, where the contour starts at a finite
point). As for the quadrature, the tolerance is relative to the
$L^1$ norm of the integrand on the whole contour, so that the
expansion at a subdominant saddle needs fewer terms, and a piece
below the tolerance is only bounded. The expansions replace the
quadrature when their estimated cost (from $n$ and the precision) is
lower; otherwise the contour is integrated as before. With the
expansions, the Bessel functions of exact arguments $|\nu|, |z|
\gtrsim 10^5$ take 10 to 40 percent less time in total at 53 bits
(where most of the remaining time is spent planning the contour), and
are 2 to 8 times faster at 333 bits and, for $|\nu|, |z| \gtrsim
10^8$, at 1000 bits.

## Choice of method {#choice}

`acb_hypgeom_bessel_j`, `acb_hypgeom_bessel_y`, `acb_hypgeom_bessel_i`,
`acb_hypgeom_bessel_k` (and the scaled versions),
`acb_hypgeom_gamma_upper` and `acb_hypgeom_gamma_lower` choose between
the convergent series, the asymptotic expansion and the integration
with the same scheme. The methods are tried at a working precision
that doubles from one round to the next; each result is a rigorous
enclosure, so the best is kept and the first that is accurate enough
is returned.

- Unless the integration is expected to be cheaper, the first attempt
  is the classical implementation
  (`acb_hypgeom_bessel_j_nointegration` and so on) at the precision
  raised by the estimated cancellation of the series, or the
  asymptotic expansion directly where it reaches the precision (for
  $K$, $Y$ and the incomplete gamma functions; the implementations of
  $J$ and $I$ change from the series to the expansion by themselves,
  the series being cheaper up to about $|z| = p/2$). For the Bessel
  functions the expansion reaches $p$ bits when $2.9 |z| \ge p$ and
  $|z| > 1.9 |\nu| + 2$, and the series loses about $1.45 (|z| -
  |\operatorname{Im} z|)$ bits for $J$ and $Y$, $1.45 (|z| -
  |\operatorname{Re} z|)$ for $I$ and $2.9 \operatorname{Re} z$ for
  $K$ (less if the terms grow less than that, for large $|\nu|$). For
  the incomplete gamma functions the number of terms of the expansion
  is estimated from the ratios $|s - k| / |z|$ of its terms, and the
  cancellation is that of the lower series in the better of its two
  forms, plus $\log_2 |\Gamma(s) / \Gamma(s,z)|$ for the upper
  function. With these estimates, small orders and arguments are
  normally settled by this first attempt.
- In each further round, the series and the expansion are tried in the
  order of their estimated cost; the series recovers what remains of
  its cancellation as the precision doubles, and the expansion either
  reaches the precision or does not.
- The integration is tried *first* if its estimated cost is below that
  of the cheaper of the expansion (where that converges) and the
  series at the precision raised by its cancellation, and otherwise as
  soon as it is cheaper than the next round. Its cost is modelled as a
  constant per descent path which depends on the precision only: about
  45 microseconds at 64 bits and 200 at 190 bits, and growing
  quadratically from there; beyond 200 bits and for exact input, the
  cost of the saddle-point expansion ([](#expansions)) is used
  instead where that is lower.
- For inexact input, the accuracy of the series is limited by the
  input balls passing through its cancellation. A result within a few
  bits of the accuracy that the input permits (estimated from the
  condition numbers, with the accuracy of an order $\nu$ counted
  relative to $\max(1, |\nu|)$) is accepted.
- If only $z$ is inexact (and the integration is not expected to be
  cheaper anyway), direct evaluation on the ball is kept only if it is
  within two bits of that accuracy, which is the case for the
  asymptotic expansion well inside its region of convergence and for
  the series without cancellation. Otherwise the function is evaluated
  at the midpoint $m$ of $z$ and the error $r \sup |\partial_z F|$ is
  added, where $r$ is the radius. For the incomplete gamma functions $\partial_z F = \mp
  z^{s-1} e^{-z}$ is evaluated on the ball. A Bessel function $g$
  (or $e^{\sigma z}$ times $I_\nu$ or $K_\nu$ for the scaled functions,
  $\sigma = \mp 1$) satisfies $g'' = a g' - b g$ with
  $a = 2 \sigma - 1/z$ and $b = 1 - \nu^2/z^2$ for $J$, $Y$,
  $b = -\nu^2/z^2 - \sigma/z$ for $I$, $K$. For any constant $c > 0$, the vector $v = (g,
  c g')$ satisfies $v' = A v$ with $\|A\|_F^2 = c^{-2} + c^2 |b|^2 +
  |a|^2$, so that along the segment from $m$ to any point of the ball
  $|g'(z)| \le \sqrt{|g(m)|^2 / c^2 + |g'(m)|^2}\, \exp(r \sup \|A\|_F)$,
  where $g'(m) = \mp g_{\nu+1}(m) + (\nu/m + \sigma) g_\nu(m)$ costs
  one more function value at low precision. The constant $c \approx
  |b|^{-1/2}$ balances the two terms, and is increased where $g'$ is
  small compared to $g$ as far as the exponential factor allows. The
  bound is within a small factor of $|g'(m)|$ unless the ball is very
  wide, so the enclosure is about as tight as possible at little more
  than the cost of an evaluation with exact input. A ball crossing
  the branch cut is excluded.
- Otherwise (inexact $\nu$ or $s$, or the above falling short), the
  alternatives are the integration, which works with the input balls
  directly, and the evaluation at the midpoints of the balls with the
  error $\sup |\partial_z F|\, r_z + \sup |\partial_\nu F|\, r_\nu$
  added from the bounds of [](#bounds). Both give about the same
  radius; the integration is preferred unless it is more than about
  four times as expensive, and its result is final. These cost a few
  times as much as the series at small arguments; inputs given with a
  few tens of bits more than the precision avoid them, since the
  series on the balls is then accurate enough.
- Inputs without any accuracy are passed to the classical
  implementation directly.

The integration applies for $|\nu|, |z| < 2^{44}$. Its accuracy is
limited only by cancellation between the contributions of different
saddles, which is inherent in the function ($J_\nu$ near a zero) or
occurs for $J_\nu$ and $I_\nu$ with $\nu$ close to a negative integer,
where the function is much smaller than its two saddle contributions.

## Use in FLINT {#use}

- Evaluation: `acb_hypgeom_bessel_j_integration`,
  `acb_hypgeom_bessel_y_integration`, `acb_hypgeom_bessel_jy_integration`,
  `acb_hypgeom_bessel_h1_integration`,
  `acb_hypgeom_bessel_h2_integration`,
  `acb_hypgeom_bessel_i_integration`,
  `acb_hypgeom_bessel_k_integration`; these are called by the automatic
  functions `acb_hypgeom_bessel_j`, `acb_hypgeom_bessel_y`,
  `acb_hypgeom_bessel_jy`, `acb_hypgeom_bessel_i`,
  `acb_hypgeom_bessel_k`, `acb_hypgeom_gamma_upper`,
  `acb_hypgeom_gamma_lower` and their real wrappers in
  [`arb_hypgeom`](arb_hypgeom).
- The versions restricted to series and asymptotic expansions:
  `acb_hypgeom_bessel_j_nointegration` and the like,
  `acb_hypgeom_gamma_upper_nointegration`,
  `acb_hypgeom_gamma_lower_nointegration`.
- The implementation is in `src/acb_hypgeom/sd_plan.c` (contours and
  bounds), `sd_integrate.c` (quadrature), `sd_expand.c` (saddle-point
  expansions), `bessel_integration.c` and
  `gamma_incomplete_integration.c` (the two families), with the
  dispatchers in `bessel.c` and `gamma_incomplete.c`.
