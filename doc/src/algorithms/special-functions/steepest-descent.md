---
title: Bessel and incomplete gamma functions by steepest-descent integration
summary: Rigorous evaluation of the Bessel functions J, Y, I, K, the Hankel functions and the incomplete gamma functions for arbitrary complex parameters by Gauss–Legendre quadrature along steepest-descent contours planned in double precision, with L1 bounds along the same contours for the functions and their derivatives; the treatment of Stokes lines, infinite ladders of saddles and branch cuts; numerically generated saddle-point (Debye) expansions with rigorous error bounds for large parameters; and the automatic choice between integration, convergent series and asymptotic expansions.
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
the integrand along parts of the contour are computed rigorously; and
the integral is *evaluated* by Gauss–Legendre quadrature with certified
error bounds. Since the integrand is entire, correctness never depends
on the quality of the contour, only the accuracy and the speed:
everything computed in double precision (the contour, the quadrature
rules, tolerances and working precisions) is a choice among valid
alternatives, and every term of the error bound is computed in ball
arithmetic for the choice made. For large parameters the quadrature can
be replaced by a *saddle-point expansion* along the same contour
([](#expansions)).

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
\pi]$, so that the formulas give the principal branches.

**Incomplete gamma functions.** In the variable $w = \log t$, with
$\alpha = -1$, $\beta = 0$, $\nu = -s$, that is, $\varphi(w) = -e^w + s
w$,
$$
\Gamma(s, z) = \int_{\operatorname{Log} z}^{V^{+}(0)} e^{\varphi(w)} dw,
\qquad
\gamma(s, z) = \int_{V}^{\operatorname{Log} z} e^{\varphi(w)} dw,
$$
where $V$ is the linear valley, to the left. For $\operatorname{Re} s >
0$ this is the defining integral of $\gamma$; otherwise it is the
analytic continuation in $s$.

**Symmetries.** For real $\nu$ and $z > 0$ the two Hankel integrals are
complex conjugates, so that $J_\nu$ and $Y_\nu$ are the imaginary and
real parts of a single integral (for $z \ge \nu$; for $z < \nu$, $J_\nu$
is exponentially smaller than $Y_\nu$ and is integrated directly). The
reflection formulas for $K_{-\nu}$, $H^{(1,2)}_{-\nu}$ and, at integers,
$J_{-n}$, $Y_{-n}$, $I_{-n}$ are applied when $\operatorname{Re} \nu <
0$. The integrals are analytic in the parameters, so for an input ball
the computed enclosures hold for the analytic continuation from the
midpoint; a ball $z$ which crosses the cut $(-\infty, 0)$ is split into
its parts above and below it.

## Planning the contour {#planning}

**Saddles.** With $u = e^w$, $\varphi'(w) = 0$ is the quadratic
equation $\alpha u^2 - \nu u - \beta = 0$, so the saddles form (at
most) two families $\log u_{1,2} + 2 \pi i k$: $\pm
\operatorname{acosh}(\nu / \zeta) + 2 \pi i k$ for the Bessel
functions, and $\log s + 2 \pi i k$ for the incomplete gamma functions.
Two saddles closer than the width $|\varphi'''|^{-1/3}$ of the cubic
saddle into which they merge (the transition region $\nu \approx \pm
\zeta$) are treated as one degenerate saddle with three descent
directions. A finite starting point is a saddle with a single descent
direction.

**Descent paths.** From each saddle $w_s$ the paths of steepest
descent of $\operatorname{Re} \varphi$ are followed numerically, with
steps that resolve the local scale of the integrand. A path ends in a
*valley*, when the integrand has decreased by the target factor and the
tail bound along a ray into the valley ([](#tails)) is negligible, or
in a *lower saddle*, when it passes within the width of that saddle:
this is the situation on a Stokes line, where one saddle connects to
another.

**Routes.** The saddles and the valleys are the nodes of a graph, and
the descent paths its edges. Any path in the graph between the two
ends is a valid contour (the integrand is entire and decays in the
valleys), and its cost is the sum of the $L^1$ bounds of its edges. The
cheapest route is found with the Bellman–Ford algorithm, tracing
saddles lazily in order of increasing $\operatorname{Re}
\varphi(w_s)$ until the remaining ones cannot improve the route.

### Ladders of saddles {#ladders}

When $\operatorname{Im} \nu \ne 0$, the Bessel-type problem has a
linear valley in the vertical direction $\sigma i$, $\sigma =
-\operatorname{sgn} \operatorname{Im} \nu$. The valleys $V^{\pm}(k + j
\sigma)$, $j \to \infty$, lead to it, separated by periodic copies of
the same saddles, whose contributions decrease geometrically. When
$|\operatorname{Im} \nu|$ is large compared to $|\zeta|$, a contour
between two valleys can be forced through the linear valley, and thus
over infinitely many saddles.

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
them is connected to the linear valley, directly or by
[](#thm-ladder), and $\int_a^b = \int_a^{\text{lin}} -
\int_b^{\text{lin}}$. This is the contour-integral form of connection
formulas such as $K_\nu = \frac{\pi}{2} (I_{-\nu} - I_\nu) / \sin(\pi
\nu)$. For the incomplete gamma functions with large $|\operatorname{Im}
s|$ the contributions of the infinitely many saddles sum to
$\Gamma(s)$; the planner then looks for a route for the complementary
function and uses $\Gamma(s, z) + \gamma(s, z) = \Gamma(s)$, normally
without cancellation, since the descent path from $\operatorname{Log}
z$ leads to the valley of the smaller of the two.

Where the landscape is flat (small parameters), the descent flow does
not lead anywhere within the range of double precision, and classical
contours are used instead: rays along the centres of the valleys
joined by segments, or for $\Gamma(s, z)$ a segment to the real axis
followed by the axis. The convergent series work best there, so this is
only a fallback.

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

With $g = 1$ this bounds the function itself. Since $\partial_\zeta
e^{\varphi} = \sinh(w) e^{\varphi}$ and $\partial_\nu e^{\varphi} = -w
e^{\varphi}$, the weights $g = \cosh x \ge |\sinh w|$ and $g = |w|$
bound the derivatives with respect to the argument and the order. Along
a path of steepest descent the integrand has constant phase, so for a
contour through a single saddle the bound with $g = 1$ equals the
integral; the bounds are tight up to the cancellation between
different saddles, that is, they bound the envelope of an oscillating
function.

**Panels.** On a segment $[w_a, w_b]$ of the contour with midpoint $m$
and $h = (w_b - w_a) / 2$, Taylor's theorem gives for $t \in [-1, 1]$
$$
\operatorname{Re} \varphi(m + h t) \le \operatorname{Re} \varphi(m) +
\operatorname{Re}(\varphi'(m) h)\, t + \tfrac{1}{2}
\operatorname{Re}(\varphi''(m) h^2)\, t^2 + \tfrac{1}{6} |h|^3 S,
\qquad S = \sup |\varphi'''|,
$$
and the integral of the exponential of a quadratic over $[-1, 1]$ is
bounded in closed form.

### Tails {#tails}

Beyond the last vertex $w_e$ of a polygon the contour continues along a
ray: horizontally into $V^{\pm}(k)$, or into the linear valley. On a
horizontal ray into $V^{+}(k)$, with $u = e^{x - x_e}$, the exponent of
[](#thm-l1) is at most $G - K(u - 1) + a \log u$, where $G$ is its
value at $w_e$, $K$ the rate at which $\operatorname{Re}(\alpha e^w)$
decreases and $a$ collects the terms which can grow; the tail is then
at most $e^{G + K} \int_1^\infty e^{-Ku} u^{a-1} du$, an incomplete
gamma function bounded in closed form. $V^{-}(k)$ follows by the
symmetry $w \to -w$, $\alpha \leftrightarrow \beta$, $\nu \to -\nu$.
Into the linear valley the exponent decreases at least linearly, and
the tail is at most $e^{G'}/\mu$ with $\mu$ that rate. The weights $g$
are carried along. The bound also verifies that the ray stays in the
valley which the plan requires.

**Arithmetic.** While the contour is planned, the bounds are
*estimated* in double precision; nothing depends on these estimates
except the quality of the contour. Where a bound enters the radius of a
result it is *computed* in ball arithmetic (`sd_bound.c`): $\varphi$ and
its derivatives are evaluated in `acb` arithmetic with the input balls
themselves as parameters, which takes care of the suprema over the
balls, of rounding, and of the cancellation between the terms of
$\varphi$ and the reference value. The polygon between two vertices can
be replaced by the straight segment between them wherever the bound for
the segment is not much larger, which leaves some 15 to 20 segments per
contour.

### Large parameters {#frames}

Doubles resolve the saddles and the values of $\varphi$ for
$|\alpha|, |\nu| < 2^{44}$ or so. Beyond, the widths of the saddles,
$|\varphi''|^{-1/2}$, fall below the spacing of the doubles around
them, and $\varphi$, whose terms are of the size of the parameters,
cannot be evaluated to an absolute accuracy below 1. Everything above
carries over if the points are taken relative to a *base point* $w_0$
close to the saddle: with $w = w_0 + \delta$,
$$
\varphi(w_0 + \delta) - \varphi(w_0)
  = (a + b)(\cosh \delta - 1) + (a - b)(\sinh \delta - \delta) + r \delta,
$$
where $a = \alpha e^{w_0}$, $b = \beta e^{-w_0}$ and $r = \varphi'(w_0)$.
For the $\delta$ of interest the terms on the right are of the size of
the left-hand side, so that it can be evaluated in double precision in
$\delta$ without cancellation. A *frame* consists of an exact base
point $w_0$, of $a$, $b$, $a + b$ and $r$ computed in ball arithmetic
and rounded, and of an exact reference value close to $\varphi(w_0)$.
Each saddle, and the starting point, has a frame; the planner works in
double precision as before with the vertices as offsets $\delta$ in the
frame of their saddle. The rigorous bounds are computed from $\varphi$
at the exact points $w_0 + \delta$, with as many more bits as the
parameters have beyond $2^{44}$, and the quadrature evaluates the
integrand in the frame, from balls for $a$, $b$, $r$ and $\varphi(w_0)$
computed once per frame. Its working precision is therefore that of
small parameters, whatever the size of $\nu$ and $z$; the cost of an
evaluation grows only weakly with the parameters. The midpoints of the
parameters must be doubles: $|\nu|, |z| < 2^{900}$ for the Bessel
functions and $2^{860}$ for the incomplete gamma functions.

## Quadrature {#quadrature}

Each piece of the contour is integrated as $e^{c} \int e^{\varphi(w) -
c} dw$, where $c$ is $\varphi$ at the saddle of the piece, rounded, so
that $|e^{\varphi - c}| \lesssim 1$ along it; the pieces are combined
in `acb` arithmetic and the tails are added to the radius.

::: theorem {#thm-gl}
Let $f$ be analytic with $|f| \le M$ on the closed disk $|w - m| \le
|h| X$, $X > 1$, and let $\rho = X + \sqrt{X^2 - 1}$. The error of the
$n$-point Gauss–Legendre rule for $\int_{m - h}^{m + h} f(w) dw$ is at
most
$$
\frac{64}{15} \frac{M |h|}{(\rho - 1) \rho^{2n - 1}}.
$$
:::

This is the second bound of
[](algorithms/ball-arithmetic/integration#thm-gl), which
`acb_calc_integrate` uses as well [@Joh2018a], for the segment: the
Bernstein ellipse with parameter $\rho$ of the segment is contained in
the disk. (Trefethen [@Tre2008] states the bound for the rule with $n
+ 1$ points.) For $f = e^{\varphi - c}$, $M$ follows from the exact
identity
$$
\varphi(m + u) - \varphi(m) = \varphi'(m) u + \varphi''(m) (\cosh u -
1) + \varphi'''(m) (\sinh u - u)
$$
(with $\varphi''' = \varphi' + \nu$), so that for $|u| \le R = |h| X$,
$$
\operatorname{Re}(\varphi(m + u) - \varphi(m)) \le |\varphi'(m)| R +
|\varphi''(m)| (\cosh R - 1) + |\varphi'''(m)| (\sinh R - R),
$$
with the derivatives evaluated on the input balls.

**Chords.** The rule is applied on *chords* joining vertices of the
planned polygon. For each chord the number of nodes is minimized over
$X$, and a chord is extended over the following panels for as long as
this lowers the number of nodes per panel. The rules, with $n = 1, 2,
4, 6, 8, 12, 16, 22, 32, \ldots$ nodes, are those of the cache of
`acb_calc_integrate`. All of this is decided from double precision
estimates; for the chord, $n$ and $X$ chosen, the bound is then
computed rigorously and added to the radius. The outer parts of each
polygon, where the integrand is below the tolerance, are bounded
instead of integrated. For the incomplete gamma functions one more
chord joins the double approximation of $\operatorname{Log} z$ to the
ball $\operatorname{Log} z$ itself; the sum with a ball end point
encloses the integrals over all chords from points of the ball.

**Arithmetic.** The nodes are evaluated in ball arithmetic with the
input balls as parameters, so that the result encloses the function on
the balls, and with the bits which $\varphi(w) - c$ cancels added to the
working precision: up to about 200 bits in the vectorized multi-double
ball arithmetic of [`dfloat`](dfloat), beyond (or where `dfloat` is not
supported) in `arb` arithmetic, in parallel when several threads are
available. A node where $|e^{\varphi - c}| = 2^{-k}$ needs $k$ bits
less for the same absolute error; this *tapering* saves about a third
of the work at high precision. The radius of the result due to the
input balls is about that of the weighted $L^1$ norms of [](#bounds),
and the tolerance and working precision are relaxed accordingly, so
that wide balls cost less than exact input.

## Saddle-point expansions {#expansions}

For large parameters the integrand is concentrated near each saddle,
and the integral over the two descent paths of a simple saddle $w_s$
can be expanded there. This is the expansion of Debye type, with the
coefficients generated numerically for the given parameters instead of
as polynomials, so that no table is needed, complex parameters need no
special treatment, and the error bound holds for the parameters at
hand.

Let $c = \varphi(w_s)$ and $A = \varphi''(w_s)$. By the identity of
[](#quadrature) at $m = w_s$,
$$
\varphi(w_s + u) - c = P(u) = A (\cosh u - 1) + \nu (\sinh u - u).
$$
On the two descent paths $P(u) = -\tau^2$ with $\tau$ real, and $u(\tau)
= \sum_{k \ge 1} b_k \tau^k$ is analytic near $\tau = 0$, with $b_1 =
\sqrt{-2/A}$. The integral over the piece is
$$
e^{c} \int e^{-\tau^2} u'(\tau) \, d\tau \sim e^{c} \sum_{k
\text{ odd}} k \, b_k \, \Gamma(k/2).
$$
The coefficients follow from $P'(u)\, u' = -2 \tau$ by comparing
coefficients, $O(n^2)$ operations for $n$ of them.

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
inequality.

**Error bound.** For $0 < \tau_0 < \rho$, the difference between the
integral over each descent path and $n$ terms consists of the
truncation of the series on $[0, \tau_0]$, bounded with the
coefficient bounds of [](#thm-saddle); the terms integrated beyond
$\tau_0$; and the integral over the rest of the path, which is replaced
(the integrand being entire) by a segment from an enclosure of $w_s +
u(\tau_0)$ to a vertex of the planned polygon followed by the rest of
the polygon, bounded as in [](#bounds). The expansion reaches a
precision of about $\rho^2 / \log 2 \approx |A|/\log 2$ bits, and
needs fewer terms the larger $|A|$ is.

It is implemented for the Bessel kernel with exact parameters and
contours made of pairs of descent paths through simple saddles, and
replaces the quadrature where its estimated cost is lower, which is
the case mainly at higher precision (2 to 8 times faster at 333 bits
for $|\nu|, |z| \gtrsim 10^5$).

## Choice of method {#choice}

`acb_hypgeom_bessel_j`, `acb_hypgeom_bessel_y`, `acb_hypgeom_bessel_i`,
`acb_hypgeom_bessel_k` (and the scaled versions),
`acb_hypgeom_gamma_upper` and `acb_hypgeom_gamma_lower` choose between
the convergent series, the asymptotic expansion and the integration by
the same scheme. For each method the number of guard bits it needs and
its cost are estimated in double precision, and the cheapest is tried
first at the precision raised by its guard bits. Each result is a
rigorous enclosure, so the estimates only affect the running time: the
best result is kept, and the first which is accurate enough is
returned.

**Guard bits.** With exact input, a method which loses $c$ bits to
cancellation gives full accuracy when evaluated with $c$ more bits.

- A convergent series loses the size of its largest term relative to
  the function. The largest term follows from the ratios of the terms,
  and the size of the function from the leading-order (Debye)
  approximations, $\log |J_\nu(z)| \approx \operatorname{Re} \xi$
  beyond the turning point and $|\operatorname{Re} \xi|$ in the
  oscillatory region, with $\xi = \sqrt{\nu^2 - z^2} - \nu \log((\nu +
  \sqrt{\nu^2 - z^2}) / z)$, and similarly for the other functions. For
  small orders this is the classical $1.45 (|z| - |\operatorname{Im}
  z|)$ bits for $J$ and $Y$; for $\nu > z$ the alternating series for
  $J$ loses twice the growth of its terms, whereas the series for $I$
  loses nothing. For $I$ and $K$ with complex orders the sizes are
  $e^{\pm \operatorname{Re} \eta}$ with $\eta = \sqrt{\nu^2 + z^2} + \nu
  \log(z / (\nu + \sqrt{\nu^2 + z^2}))$ in the sector where this holds
  [DLMF, 10.41(iv)](https://dlmf.nist.gov/10.41#iv).
- The asymptotic expansions of the Bessel functions have terms which
  grow before they decrease when $|\nu|^2$ is not small compared to
  $|z|$; the growth is lost to cancellation, and the expansion is used
  only where its smallest term, times the factor by which the error
  bound exceeds it, reaches the precision. For a half-integer order the
  expansions terminate and are exact, and are used for any $z \ne 0$
  where they are cheapest.
- For the incomplete gamma functions, the series $\gamma(s,z) = s^{-1}
  z^s e^{-z} \sum_k z^k / (s+1)_k$ loses nothing for $|z| < |s|$ and
  about $1.45 (|z| - \operatorname{Re} z)$ bits otherwise, less where
  the end point dominates the integral; the upper function as
  $\Gamma(s) - \gamma(s,z)$ loses $\log_2 |\Gamma(s) / \Gamma(s,z)|$
  bits more. For a positive integer $s$ the expansion of $\Gamma(s, z)$
  terminates, $\Gamma(s,z) = (s-1)!\, e^{-z} \sum_{j<s} z^j / j!$, and
  loses nothing for $s \le |z|$. For $s = 1, 2, 3$ the closed forms are
  returned at once.

**Costs.** The numbers of terms follow from the same ratios, and the
cost per term from timings; with exact short $\nu$ and $z$ the series
are summed by binary splitting, which costs less. The cost of the
integration is modelled per descent path as a function of the
precision only (the saddle-point expansion where that is cheaper).
These are costs with the quadrature nodes at hand: the nodes are
computed on the first call at each precision (about a second at 3333
bits) and cached, and a tenth of this is added to the cost of the
quadrature, so that the series are preferred for a single evaluation
and the integration where it is much faster.

**Order of the attempts.** In the easiest cases the estimates cost as
much as the function, and a first attempt with a few guard bits comes
before them (small $|z|$, or large $|z|$ for $K$ and the upper
incomplete gamma function). Otherwise the cheapest method comes first,
and a result which falls short by a known number of bits is repeated
with that many more; else the precision is doubled. The integration is
tried before any attempt which is estimated to cost more. On surveys
of exact inputs with $|\nu|$ and $|z|$ up to 30 times the precision,
from 333 to 30000 bits, the automatic functions take a median of about
the time of the best single method run with exactly the guard bits it
needs, and at most about four times.

**Inexact input.** A result within a few bits of the accuracy which the
input allows is accepted (that accuracy is estimated from condition
numbers). For input which is as accurate as the precision, up to 1000
bits, a few more bits may be missing in a first result on the input
balls, since the alternatives cost several times as much there and the
caller can raise the precision. A ball $z$ which meets a branch cut
from both sides is left to the series, since no method gives more than
the jump of the function. Otherwise the routes are, in order of cost:

- *A first result which is already tight.* Close to a zero of the
  function the accuracy which the input allows is lower than its
  estimate. A first result on the balls is final if its radius is
  within three bits of $r\, |\partial_z F|$ or $r\, |\partial_\nu F|$,
  the image of the input to first order, computed from one more value
  at low precision.
- *Only $z$ inexact.* The function is evaluated at the midpoint $m$ of
  $z$, and $r \sup |\partial_z F|$ is added. For the incomplete gamma
  functions $\partial_z F = \mp z^{s-1} e^{-z}$ is evaluated on the
  ball. A Bessel function $g$ satisfies $g'' = a g' - b g$ (with $a =
  -1/z$, $b = 1 - \nu^2/z^2$ for $J$ and $Y$, and the analogues for $I$,
  $K$ and the scaled functions), so that $v = (g, c g')$ satisfies $v' =
  A v$ with $\|A\|_F^2 = c^{-2} + c^2 |b|^2 + |a|^2$, and along the
  segment from $m$ to any point of the ball
  $$
  |g'(z)| \le \sqrt{|g(m)|^2 / c^2 + |g'(m)|^2}\, \exp(r \sup \|A\|_F),
  $$
  where $g'(m) = \mp g_{\nu+1}(m) + (\nu/m) g_\nu(m)$ costs one more
  value at low precision (for $Y$ and $K$ at integer orders, from $J$
  or $I$ and the Wronskians). With $c \approx |b|^{-1/2}$ the bound is
  within a small factor of $|g'(m)|$ unless the ball is very wide.
- *Monotone functions.* For real arguments where the function is
  monotone in each inexact variable, the range over the input box is
  the interval between the values at two corners, the optimal
  enclosure, at the cost of two evaluations with exact input. This
  holds for $I_\nu(x)$, $K_\nu(x)$, $e^x K_\nu(x)$, and the incomplete
  gamma functions with $x > 0$, including $Q$ and $P$ (for $s > 0$),
  for which it also avoids the dependency between $\Gamma(s,x)$ and
  $\Gamma(s)$. It is used for input with at most 64 accurate bits;
  narrower balls at high precision are better served by the next
  route.
- *An inexact order.* The functions are entire in the order, so with a
  bound $M$ for $|F|$ on the disk $|\nu' - m| \le R$, Cauchy's
  estimates give
  $$|F(\nu, z) - F(m, z)| \le r\, |\partial_\nu F(m, z)| + \frac{r^2 M}{(R - r)^2}, \qquad
  |\partial_\nu F(m, z)| \le \frac{|F(m + h, z) - F(m, z)|}{h} + \frac{h M}{(R - h)^2}.$$
  With $r \ll h \ll R$ the bound is $r\, |\partial_\nu F|$ up to terms
  of higher order, so that a crude $M$ is enough: the bounds of the
  integral representations, such as $|K_\nu(z)| \le K_{|\operatorname{Re}
  \nu|}(\operatorname{Re} z)$, carried to orders with a negative real
  part by the recurrences. The second value needs $\log_2(h/r)$ fewer
  bits than the first.
- *Otherwise* the integration, which works with the input balls
  directly, or the evaluation at the midpoints with the derivative
  bounds of [](#bounds), whichever is estimated to cost less; both give
  about the same radius.

The integration applies for $|\nu|, |z| < 2^{900}$ ($2^{860}$ for the
incomplete gamma functions). Its accuracy is limited only by
cancellation between the contributions of different saddles, which is
inherent in the function ($J_\nu$ near a zero) or occurs for $J_\nu$
and $I_\nu$ with $\nu$ close to a negative integer.

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
- The profile program `src/acb_hypgeom/profile/p-catalog.c` compares
  the automatic functions with the integrals and with the versions
  without them on a catalog of inputs; see
  [](algorithms/special-functions/hypergeometric#profiling).
- The implementation is in `src/acb_hypgeom/sd_plan.c` (contours and
  estimates, in double precision, and the frames for large
  parameters), `sd_bound.c` (the rigorous bounds),
  `sd_integrate.c` (quadrature), `sd_expand.c` (saddle-point
  expansions), `bessel_integration.c` and
  `gamma_incomplete_integration.c` (the two families), with the
  dispatchers in `bessel.c` and `gamma_incomplete.c`.
