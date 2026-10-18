---
title: Rigorous numerical integration
summary: How acb_calc_integrate computes enclosures of integrals of analytic functions along line segments in the complex plane, combining direct interval enclosures, Gauss–Legendre quadrature with error bounds from bounds of the integrand on Bernstein ellipses, and adaptive bisection; with the stopping rules and the role of the user's evaluation function.
complexity: For an integrand that is analytic in a neighbourhood of the path and bounded there, $O(p)$ evaluations at $p$ bits (Gauss–Legendre degree $O(p)$ on $O(1)$ subintervals); with endpoint singularities or near-singularities, $O(p^2)$ evaluations in the worst case. Nodes and weights cost $O(p^2)$ operations per degree and are cached.
order: 110
---

Let $f$ be a function that is analytic on a neighbourhood of the
segment $[a, b] \subset \C$ (except possibly at isolated points where it
is integrable) and that can be evaluated in ball arithmetic. The goal is
a ball containing $I = \int_a^b f(t)\,dt$ whose radius is close to a
requested tolerance. Classical numerical integration estimates the error
from differences between quadrature rules, which can be fooled; in ball
arithmetic, the error must be bounded rigorously, and FLINT does this
with bounds of $|f|$ on complex neighbourhoods of the path, following
Petras [@Pet2002] and [@Joh2018a].

## Direct enclosure {#direct}

If $B$ is a ball (a box in the complex plane) containing the segment
from $a$ to $b$, and $F$ is a ball containing $f(B)$, then $I \in (b -
a) F$. This direct enclosure needs only one evaluation of $f$ on a wide
ball and is valid for any integrable $f$, but it converges only linearly
as the segment shrinks: its radius is about $|b - a| \cdot |f'| \cdot
|b - a|$. It is used for the initial estimate on each subinterval and
to stop at segments where $f$ is so small or the segment so short that
the enclosure already meets the tolerance.

## Gauss–Legendre quadrature with ellipse bounds {#gl}

For $[a, b] = [-1, 1]$, the Gauss–Legendre rule of degree $n$ is
$$Q_n = \sum_{k=1}^{n} w_k f(x_k),$$
where the nodes $x_k$ are the roots of the Legendre polynomial $P_n$ and
$w_k$ are the corresponding weights; it is exact for polynomials of
degree less than $2n$. For $\rho > 1$, let $E_\rho$ be the Bernstein
ellipse with foci $\pm 1$ and semi-axes $(\rho + \rho^{-1})/2$ and
$(\rho - \rho^{-1})/2$.

::: theorem {#thm-gl}
If $f$ is analytic inside $E_\rho$ and $|f(z)| \le M$ there, then
$$|I - Q_n| \le \frac{64 M}{15 (\rho^2 - 1) \rho^{2n}} \le \frac{64 M}{15 (\rho - 1) \rho^{2n - 1}}.$$
:::

::: proof
The first inequality is due to Trefethen [@Tre2008] (see also
[@Tre2013, Chapter 19]); it follows from the Chebyshev expansion of $f$,
whose coefficients are bounded by $2M\rho^{-k}$, and the exactness of the
rule on polynomials of degree less than $2n$. The second follows from
$(\rho^2 - 1)\rho^{2n} = (\rho - 1)(\rho + 1)\rho^{2n} \ge (\rho -
1)\rho^{2n-1}$.
:::

`acb_calc_integrate_gl_auto_deg` maps the segment to $[-1, 1]$ by $t =
m + \delta x$ with $m = (a + b)/2$, $\delta = (b - a)/2$, and tries the
ellipses with $(\rho + \rho^{-1})/2 = X$ for $X = 2, 4, 8, 32, 512,
\ldots$ (that is, $X = 2^{e+1}$ with $e = 0, 1, 2, 4, 8, \ldots$, $e <
p$). For each $X$, it evaluates $f$ once on the box $m + \delta([-X, X]
+ [-Y, Y] i)$ with $Y = \sqrt{X^2 - 1}$, which contains the image of
$E_\rho$, with the flag that asks for a proof of analyticity (the
evaluator must return a non-finite value if $f$ might not be analytic in
the box; see [](acb_calc/integrate)). If the result is finite, $M =
|\delta| \max|f|$ bounds the integrand of the transformed integral, and
the smallest degree $n$ in the list $1, 2, 4, 6, 8, 12, 16, 22, 32, 46,
64, \ldots$ (steps of about $\sqrt 2$, up to the degree limit) for which
the bound of [](#thm-gl) is below the tolerance is recorded. The search
stops at the first non-finite evaluation, since larger boxes contain the
smaller ones, and otherwise picks the smallest degree over all $X$. The
rule is then evaluated in ball arithmetic and the bound added to the
radius. If no degree works, the function reports failure and the caller
bisects.

The nodes and weights are computed with
`arb_hypgeom_legendre_p_ui_root` (an initial isolating interval following
Petras, refined by interval Newton iteration on $P_n$ with doubling
precision) and cached per degree, at a precision that at
least doubles plus 30 bits when more precision is needed, in a
thread-local cache; only half of them are stored, using the symmetry
$x_{n+1-k} = -x_k$. Computing them takes $O(n^2)$ operations, but in an
adaptive integration the same few degrees are used many times.

The choice of ellipses reflects the trade-off in [](#thm-gl): a larger
$\rho$ gives faster convergence, but $M$ grows, and singularities of $f$
limit $\rho$. A pole at distance $d$ from a segment of half-length $h$
allows $\rho \approx 1 + d/h$ for small $d/h$, so that $n \approx
(p \ln 2)/(2 d/h)$ nodes are needed; bisection reduces $h$ and with it
the degree.

## Adaptive subdivision {#adaptive}

::: algorithm Adaptive integration {#alg-adaptive}
Input: $f$, $a$, $b$, a relative goal $g$ (bits), an absolute tolerance
$\tau$, limits on the subdivision depth, the number of evaluations and
the degree, and a precision $p$.

1. Compute the direct enclosure $v$ on $[a, b]$ and set the working
   tolerance to $\tau' = \max(\tau, 2^{-g} |v|_{\mathrm{low}})$, where
   $|v|_{\mathrm{low}}$ is a lower bound for the absolute value.
2. Keep a stack (or, optionally, a priority queue ordered by radius) of
   subintervals with their direct enclosures. Repeatedly take a
   subinterval $J$:
   - if the radius of its direct enclosure is below $\tau'$, or $J$ is
     so short that its endpoints overlap at precision $p$, or a limit has
     been reached, add the direct enclosure to the sum;
   - otherwise, if Gauss–Legendre quadrature with the tolerance $\tau'$
     succeeds on $J$, add its result to the sum;
   - otherwise, bisect $J$, compute the direct enclosures of both
     halves and push them, the one with the larger radius on top.
3. After each new enclosure $u$, update $\tau' \leftarrow \max(\tau',
   2^{-g} |u|_{\mathrm{low}})$.
4. Return the sum, and a status that indicates whether a limit was
   reached.
:::

The result is always a valid enclosure: every subinterval contributes
either a direct enclosure or a Gauss–Legendre result with a rigorous
error bound, and the subintervals partition $[a, b]$. The tolerance and
the limits only affect the work and the radius. The update of $\tau'$
implements the relative goal: the integral over a subinterval with
large values gives a lower bound for the magnitude of the result only
heuristically (the integrals over different subintervals may cancel),
but this only makes the radius larger than requested in such cases, not
the result wrong. If the direct enclosure of a subinterval is real, the
imaginary part of the Gauss–Legendre result is set to zero, since the
integral is then real.

The default limits are a degree of at most $0.5 \min(g, p) + 60$, at
most $1000p + p^2$ evaluations and a depth of at most $2p$ (the stack
size; with the depth-first strategy, this bounds the nesting of
bisections). Bisection handles singularities at or near the path: near
an integrable endpoint singularity, the intervals shrink geometrically
towards the singularity, and each level costs one Gauss–Legendre
evaluation of low degree on the subinterval adjacent to the one
containing the singularity, so about $p$ levels and $O(p)$ evaluations
per level are needed. Functions that are not analytic, such as $|t|$ or
$\lfloor t \rfloor$, can be integrated if the evaluator reports
non-analyticity on boxes that contain a discontinuity, so that the
integration falls back to direct enclosures on the subintervals that
contain such a point. Since the direct enclosure on an interval of
length $h$ containing a jump has radius of order $h$ (and of order $h^2$
for a kink), about $p$ (respectively $p/2$) bisection levels are needed
around each such point; it is cheaper to pass these points as endpoints
of separate integrals.

## Use in FLINT {#use}

- `acb_calc_integrate` ([](#alg-adaptive)),
  `acb_calc_integrate_gl_auto_deg` ([](#gl)) and
  `acb_calc_integrate_opt_t` ([](acb_calc/integrate)).
- `acb_calc_integrate_taylor` is an older method that uses Taylor
  expansions with Cauchy bounds (`acb_calc_cauchy_bound`) at uniformly
  spaced points; it needs derivatives of the integrand and a known
  distance to the singularities.
- The Gauss–Legendre nodes come from `arb_hypgeom_legendre_p_ui_root`
  ([](arb_hypgeom)).
- Integral representations of special functions in FLINT are usually
  evaluated by other means (series, recurrences), but
  `acb_calc_integrate` is a convenient way to check them, and the
  example programs in `examples/integrals.c` cover many classical
  integrals.
