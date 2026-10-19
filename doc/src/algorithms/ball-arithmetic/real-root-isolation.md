---
title: Isolating and refining real roots
summary: Bisection with exclusion and monotonicity tests for isolating the simple real roots of an analytic function given by a ball evaluator, Newton iteration with a rigorous quadratic error bound, the interval Newton operator, and precision doubling, as used by arb_calc, arb_poly and arb_fmpz_poly.
complexity: Isolation: one evaluation of $f$ and $f'$ on an interval per examined subinterval, $O(k \cdot \mathrm{depth})$ subintervals for $k$ roots in favourable cases. Refinement to $p$ bits: $O(p)$ evaluations by bisection, $O(\log p)$ Newton steps at doubling precisions, total cost about twice that of one evaluation at $p$ bits.
order: 70
---

Let $f$ be a real function that is analytic on an interval $[a, b]$, and
suppose that for any subinterval $I$ we can compute balls containing
$f(I) = \{f(t) : t \in I\}$ and $f'(I)$, as well as the value of $f$ at
exact points. The root isolation problem asks for disjoint intervals that
each contain exactly one root of $f$ and together contain all roots in
$[a, b]$. With interval evaluations, this can be solved for simple roots
by bisection, and isolated roots can then be refined rapidly by Newton's
method with rigorous error bounds [@Moo1966; @Neu1990].

## Exclusion and isolation tests {#tests}

::: lemma Exclusion {#lem-exclusion}
If the computed enclosure of $f(I)$ does not contain 0, then $f$ has no
root in $I$.
:::

::: lemma Isolation {#lem-isolation}
Let $I = [u, v]$. If $f(u)$ and $f(v)$ are nonzero with opposite signs
and the enclosure of $f'(I)$ does not contain 0, then $f$ has exactly one
root in $I$, it lies in the interior, and it is simple.
:::

::: proof
$f'$ is continuous and does not vanish on $I$, so it has constant sign
and $f$ is strictly monotonic on $I$; by the intermediate value theorem
it has a root in $(u, v)$, which is unique by monotonicity and simple
since $f' \ne 0$ there.
:::

::: algorithm Isolation by bisection {#alg-bisect}
Input: $f$, $[a, b]$, limits on the depth, the number of evaluations and
the number of roots. Output: a sorted list of subintervals with flags.

1. Compute the signs $\sigma_a$, $\sigma_b$ of $f(a)$ and $f(b)$ (0 if
   the computed value contains zero).
2. Process $[a, b]$ recursively: for a subinterval $I = [u, v]$ with
   known endpoint signs, if [](#lem-exclusion) applies, discard $I$; if
   [](#lem-isolation) applies, output $I$ with flag 1; if a limit has been
   reached, output $I$ with flag 2 (unknown); otherwise compute the sign
   of $f$ at the midpoint $w$ and process $[u, w]$ and $[w, v]$.
:::

Every root of $f$ in $[a, b]$ lies in some output interval, because only
intervals that certainly contain no root are discarded. The intervals
are output from left to right and overlap at most in endpoints. If all
roots are simple and none of them is a bisection point, then once the
intervals are small enough, each contains at most one root and the
enclosures of $f(I)$ and $f'(I)$ are accurate enough for the tests to
succeed, so the algorithm terminates with all roots isolated. A multiple
root, or a root at a bisection point (where the sign is zero), is never
isolated and ends up in intervals of the minimal width with flag 2.

The cost is dominated by the evaluations on intervals. If the
enclosures are tight, roots at distance $\delta$ from each other are
separated after about $\log_2((b - a)/\delta)$ bisection levels; the
dependency problem of interval evaluation (the enclosure of $f(I)$ can be
much wider than the true range) makes the tests succeed only later.
Evaluating $f$ through a Taylor expansion at the midpoint of $I$
(centred forms) reduces the overestimation from $O(|I|)$ to
$O(|I|^2)$.

## Newton iteration with a quadratic error bound {#newton}

::: theorem {#thm-newton}
Let $f$ be twice continuously differentiable on an interval $J$, let $x =
[m - r, m + r] \subseteq J$ contain a root $\xi$ of $f$, and let
$$C \ge \sup_{s, t \in J} \frac{|f''(s)|}{2 |f'(t)|}.$$
Then $m' = m - f(m)/f'(m)$ satisfies $|m' - \xi| \le C r^2$.
:::

::: proof
By Taylor's theorem, $0 = f(\xi) = f(m) + f'(m)(\xi - m) + \tfrac12
f''(\eta)(\xi - m)^2$ for some $\eta$ between $m$ and $\xi$, hence
$$m' - \xi = m - \xi - \frac{f(m)}{f'(m)} = \frac{f''(\eta)}{2 f'(m)} (\xi - m)^2,$$
and $|m' - \xi| \le C |\xi - m|^2 \le C r^2$.
:::

`arb_calc_newton_step` and `_arb_poly_newton_step` compute $m'$ in ball
arithmetic at the working precision (the rounding errors are part of the
resulting ball) and add $C r^2$ to the radius. The step is accepted if
the new ball is contained in $J$ (so that the next step can use the same
$C$) and has a smaller radius. Once $C r < 1$, the radius is roughly
squared at each step, so the number of correct bits doubles.
`arb_calc_refine_root_newton` and `_arb_poly_newton_refine_root` exploit
this by running the steps at precisions that also double: the last step
at $p + O(\log C)$ bits, the one before at about half of that, and so
on, down to the accuracy of the starting ball. The total cost is then
about twice that of the last step.

The factor $C$ is computed by evaluating $f'(J)$ and $f''(J)$ on the
whole interval $J$ (`arb_calc_newton_conv_factor`,
`_arb_poly_newton_convergence_factor`). It is finite only if $f'$ is
bounded away from zero on $J$, which requires $J$ to be a narrow
interval around a simple root. (In FLINT 3.6 and
earlier `arb_calc_newton_conv_factor` returned only $C/2$, for which
[](#thm-newton) does not apply.)

## The interval Newton operator {#interval-newton}

A variant needs no bound for $f''$. For an interval $X$ with midpoint
$m$, the interval Newton operator is
$$N(X) = m - \frac{f(m)}{f'(X)},$$
where $f'(X)$ is an enclosure of the range of $f'$ on $X$.

::: theorem Interval Newton {#thm-interval-newton}
If $0 \notin f'(X)$, then every root of $f$ in $X$ lies in $N(X)$. If
moreover $N(X)$ is contained in the interior of $X$, then $f$ has
exactly one root in $X$.
:::

::: proof
For a root $\xi \in X$, the mean value theorem gives $f(m) = f(m) -
f(\xi) = f'(\eta)(m - \xi)$ with $\eta \in X$, so $\xi = m -
f(m)/f'(\eta) \in N(X)$.

For existence, let $X = [u, v]$ and suppose that $f'(X) \subset (0,
\infty)$ (the negative case is symmetric). By the mean value theorem,
$f(v) = f(m) + d (v - m)$ with $d = f'(\eta) \in f'(X)$. Since $m -
f(m)/d \in N(X)$ lies in the interior of $X$, we have $m - f(m)/d < v$,
that is $-f(m) < d (v - m)$, so $f(v) > 0$. In the same way, $f(u) =
f(m) - d'(m - u)$ with $d' \in f'(X)$ and $m - f(m)/d' > u$ give $f(u) <
0$. By the intermediate value theorem, $f$ has a root in $(u, v)$, which
is unique since $f$ is strictly increasing on $X$ [@Moo1966; @Neu1990].
:::

When the root is well isolated, $N(X)$ has width $O(|X|^2)$, so the
interval Newton method also converges quadratically.
`arb_fmpz_poly_refine_root_arb` combines it with ordinary Newton steps
and bisection: it keeps an exact interval $[\ell, h]$ that contains the
root, tries an interval Newton step (accepted if $N(X)$ lies in the
interior of $X$ and gains at least 10% in accuracy), otherwise an
ordinary Newton step whose result is accepted if a small interval around
it has a sign change at its endpoints, and otherwise bisects. The
working precision follows twice the current accuracy, so the cost is
again dominated by the last steps.

## Use in FLINT {#use}

- [](#alg-bisect): `arb_calc_isolate_roots`, with `arb_calc_partition`
  for the bisection step; `arb_calc_refine_root_bisect` refines by
  bisection.
- [](#thm-newton): `arb_calc_newton_step`, `arb_calc_refine_root_newton`
  and `arb_calc_newton_conv_factor` for functions given by evaluators;
  `_arb_poly_newton_step`, `_arb_poly_newton_refine_root` and
  `_arb_poly_newton_convergence_factor` for polynomials.
- [](#thm-interval-newton): `arb_fmpz_poly_refine_root_arb`, used by
  `arb_fmpz_poly_real_roots` after exact isolation by Descartes' rule
  ([](algorithms/polynomials/real-roots-integers)).
- The zeros of the Riemann zeta function and of Dirichlet $L$-functions
  on the critical line are isolated by sign changes of the Hardy
  $Z$-function and refined in [](acb_dirichlet).
