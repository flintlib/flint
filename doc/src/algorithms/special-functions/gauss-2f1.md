---
title: The Gauss hypergeometric function
summary: Evaluation of 2F1(a,b;c;z) on the whole complex plane — the linear fractional transformations and how one is chosen, limits at integer parameter differences by power series perturbation, polynomial and special cases, analytic continuation by Taylor series with Cauchy–Kovalevskaya majorant bounds near exp(±πi/3), and numerical integration for real arguments.
order: 47
---

The Gauss hypergeometric function
$$
{}_2F_1(a,b;c;z) = \sum_{k=0}^\infty \frac{(a)_k (b)_k}{(c)_k}
\frac{z^k}{k!}
$$
converges only for $|z| < 1$ and slowly when $|z|$ is close to 1. It is
extended to $\mathbb{C} \setminus [1, \infty)$ by analytic continuation. The
regularized function $\mathbf{F}(a,b;c;z) = {}_2F_1(a,b;c;z)/\Gamma(c)$
is entire in the parameters. FLINT evaluates both with the scheme
of [@Joh2016], implemented in `acb_hypgeom_2f1`: special cases, then a
linear fractional transformation that maps $z$ into a disk where the
series converges quickly, or a short analytic continuation from such a
disk; the limit cases of the transformation formulas are handled with
power series arithmetic. The series themselves are summed and bounded
with the methods of [](algorithms/special-functions/hypergeometric).

## Transformations {#transformations}

The six maps $z \mapsto z, \frac{z}{z-1}, \frac1z, \frac{1}{1-z}, 1-z,
1 - \frac1z$ permute the singular points $0, 1, \infty$. Each gives a
formula for ${}_2F_1$ in terms of series in the new variable $w$
[@NIST2012, §15.8]:

| Index | $w$ | Formula | Singular when |
|---|---|---|---|
| 0 | $z$ | direct series | — |
| 1 | $z/(z-1)$ | $(1-z)^{-a}\,{}_2F_1(a, c-b; c; w)$ (Pfaff) | never |
| 2 | $1/z$ | two series in $w$ with prefactors $(-z)^{-a}$, $(-z)^{-b}$ | $b - a \in \Z$ |
| 3 | $1/(1-z)$ | two series with prefactors $(1-z)^{-a}$, $(1-z)^{-b}$ | $b - a \in \Z$ |
| 4 | $1-z$ | two series, one with prefactor $(1-z)^{c-a-b}$ | $c - a - b \in \Z$ |
| 5 | $1-1/z$ | two series with prefactors $z^{-a}$, $z^{a-c}(1-z)^{c-a-b}$ | $c - a - b \in \Z$ |

For example, transformation 4 reads, in regularized form,
$$
\frac{\sin(\pi(c-a-b))}{\pi} \mathbf{F}(a,b;c;z) =
\frac{\mathbf{F}(a,b;a+b-c+1;1-z)}{\Gamma(c-a)\Gamma(c-b)} -
(1-z)^{c-a-b}\frac{\mathbf{F}(c-a,c-b;c-a-b+1;1-z)}{\Gamma(a)\Gamma(b)}.
$$
The regularized form of every formula is valid for all parameters
except where the factor $\sin(\pi(b-a))$ or $\sin(\pi(c-a-b))$ on the
left vanishes. `acb_hypgeom_2f1_transform` evaluates the right-hand
side with regularized series, divides by the sine and multiplies by
$\Gamma(c)$ for the ordinary function. The Pfaff transformation is
used with the regularized series and multiplied by $\Gamma(c)$
afterwards, because the identity in non-regularized form would read
$0 \cdot \infty$ at nonpositive integer $c$.

**Choice of transformation.** The cost of a series is governed by $|w|$:
for $|w| \le 1/4$ each term gains at least 2 bits. `acb_hypgeom_2f1_choose`
computes $|z|$, $|z/(z-1)|$, $|1/z|$, $|1/(1-z)|$, $|1-z|$ and $|1 -
1/z|$ in double precision and returns:

- 0 (direct) if $|z| \le 1/4$, and 1 (Pfaff) if $|z/(z-1)| \le 1/4$;
- the better of 0 and 1 if either is at most $3/4$;
- otherwise the index with the smallest $|w|$, if it is at most $3/4$;
- otherwise 6, the corner case.

The union of the disks $|w| \le 3/4$ for the six maps covers the plane
except two small regions around $e^{\pm \pi i/3}$, the fixed points of
the group of transformations, where every $|w|$ is close to 1.

## Limits at integer parameter differences {#limits}

If $b - a$ (for transformations 2 and 3) or $c - a - b$ (for 4 and 5) is
an integer, both sides of the transformation formula vanish to first
order. Writing $a = a_0 + \varepsilon$, the function value is the
limit
$$
\mathbf{F}(z) = \lim_{\varepsilon \to 0} \frac{\pi\, R(\varepsilon)}{\sin(\pi
\delta(\varepsilon))} = \frac{\pi R'(0)}{\pi\delta'(0)\cos(\pi\delta(0))},
$$
where $R$ is the right-hand side and $\delta$ the parameter difference.
`acb_hypgeom_2f1_transform_limit` evaluates $R$ with $a$ replaced by
the power series $a + \varepsilon \in \mathbb{C}[[\varepsilon]]/(\varepsilon^2)$,
using `acb_hypgeom_2f1_series_direct`, `acb_poly_rgamma_series` and
`acb_poly_pow_series`, and divides the coefficients of $\varepsilon$ of
the numerator and of $\sin(\pi\delta)$. The power series tail bounds of
[](algorithms/special-functions/hypergeometric#series-tail) make the
derivatives rigorous. For transformations 2 and 3 the roles of $a$ and
$b$ are exchanged if $b - a < 0$, to avoid a series with a large
negative parameter.

The transformation code can only recognize an integer difference that
is proved by the balls, that is, when the difference is an exact
integer. If the parameters are inexact but known to differ by an
integer, the caller can pass the flags `ACB_HYPGEOM_2F1_AB` or
`ACB_HYPGEOM_2F1_ABC` to `acb_hypgeom_2f1`, which then uses the limit
formula unconditionally (see [](acb_hypgeom/2f1)).

## Special cases {#special}

Before choosing a transformation, `acb_hypgeom_2f1` handles:

- **Polynomials.** If $a$ or $b$ is a nonpositive integer $-n$ with $n <
  p$, the series terminates and is summed directly for any $z$; the
  rounding errors of a long polynomial with large terms are compensated
  only by the working precision chosen by the caller. If $c - a$ or $c
  - b$ is a nonpositive integer and $c$ is exact, the Pfaff
  transformation produces a terminating series.
- **Regularized zeros.** If $c$ is a nonpositive integer and $a$ or $b$
  is a nonpositive integer that is at least $c$, $\mathbf{F} = 0$.
  If $a = c$ or $b = c$, $\mathbf{F} = (1-z)^{-b}/\Gamma(c)$ or
  $(1-z)^{-a}/\Gamma(c)$.
- **$z = 1$.** Gauss's summation $\mathbf{F}(a,b;c;1) = \Gamma(c-a-b)/
  (\Gamma(c-a)\Gamma(c-b))$ for $\operatorname{Re}(c-a-b) > 0$; the
  function is infinite or undefined otherwise, and the result is
  indeterminate.

## The corner case {#corner}

Near $z = e^{\pm \pi i/3}$ no transformation gives $|w|$ much below 1.
`acb_hypgeom_2f1_corner` instead continues the solution of the
hypergeometric differential equation analytically. The function $f(x)
= {}_2F_1(a,b;c;z_0+x)$ satisfies
$$
f''(x) = -\frac{(z_0+x)(a+b+1)-c}{(z_0+x)(z_0-1+x)} f'(x) -
\frac{a b}{(z_0+x)(z_0-1+x)} f(x).
$$
Given $f(0)$ and $f'(0)$, the Taylor coefficients $f_{[k]}$ follow
from a three-term recurrence, and evaluating the truncated Taylor
series gives $f$ and $f'$ at a nearby point. The corner case starts
with the direct series at $z_1 = 0.375 \pm 0.625i$ (where $|z_1| \approx
0.73$; $f'(z_1) = \frac{ab}{c}\,{}_2F_1(a+1,b+1;c+1;z_1)$ is also a
direct series), continues to $z_2 = 0.5 \pm 0.8125i$ and from there to
$z$, with the sign of the imaginary parts taken from $z$.

::: theorem Majorant bound {#thm-majorant}
Let $\nu \ge \max(1/|z_0|, 1/|z_0 - 1|)$, $M_0 = 2\nu|ab|$, $M_1 =
\nu|a+b+1| + 2|c|$, $N \ge \max(\sqrt{2M_0}, 2M_1)/\nu$ and $A \ge
\max(|f(0)|, |f'(0)|/(\nu(N+1)))$. Then the Taylor coefficients of $f$
at 0 satisfy
$$
|f_{[k]}| \le A \binom{N+k}{k} \nu^k.
$$
:::

This follows from the Cauchy–Kovalevskaya majorant method, as outlined
by van der Hoeven [@Hoe2001]: the coefficients of the differential
equation, as power series in $x$, are majorized by those of
$$
g''(x) = \frac{N+1}{2} \left( \frac{\nu}{1-\nu x} \right) g'(x) +
\frac{(N+1)N}{2} \left( \frac{\nu}{1-\nu x} \right)^2 g(x),
$$
whose solution $g(x) = A(1 - \nu x)^{-N}$ majorizes $f$ when the
initial values do, and $g_{[k]} = A\binom{N+k}{k}\nu^k$.
`acb_hypgeom_2f1_continuation` computes these quantities in `mag_t`
arithmetic, chooses the number of terms so that the tail of the
majorant series at $|x|\nu < 1$ is below $2^{-p}$, and adds the tail
bound for $f$ and $f'$ to the result. The points must therefore be
well separated from 0 and 1 compared to the step length.

## Numerical integration {#integration}

For large real parameters, the series in all the transformations
converge slowly and cancel. If the result of the steps above has less
than half the requested accuracy and all arguments are real (and below
$2^{60}$), `acb_hypgeom_2f1` evaluates
$$
{}_2F_1(a,b;c;z) = \frac{\Gamma(c)}{\Gamma(b)\Gamma(c-b)} \int_0^1
t^{b-1}(1-t)^{c-b-1}(1-zt)^{-a}\, dt
$$
with `arb_hypgeom_2f1_integration` (for $b \ge 1$, $c - b \ge 1$ and $z
< 1$, possibly after exchanging $a$ and $b$), using the rigorous
integration of `acb_calc_integrate` [@Joh2018a], and returns the better
enclosure.

## Use in FLINT {#use}

- `acb_hypgeom_2f1` and the real `arb_hypgeom_2f1` implement the
  scheme; `acb_hypgeom_2f1_direct`, `acb_hypgeom_2f1_transform`,
  `acb_hypgeom_2f1_transform_limit`, `acb_hypgeom_2f1_corner`,
  `acb_hypgeom_2f1_continuation` and `acb_hypgeom_2f1_choose` expose
  the parts.
- `acb_hypgeom_pfq` calls `acb_hypgeom_2f1` for $(p, q) = (2, 1)$.
- Functions defined through ${}_2F_1$: the incomplete beta function
  `acb_hypgeom_beta_lower`, the Legendre functions
  `acb_hypgeom_legendre_p` and `acb_hypgeom_legendre_q`, the Jacobi,
  Gegenbauer and Chebyshev functions of non-integer degree
  (see [](algorithms/special-functions/legendre)).
- The complete elliptic integrals of [`acb_elliptic`](acb_elliptic) are
  ${}_2F_1$ functions but are computed with the arithmetic-geometric
  mean instead ([](algorithms/special-functions/elliptic-modular)).
