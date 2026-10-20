---
title: Counting and isolating real roots over Z
summary: Root bounds, Descartes' rule of signs and the Vincent-Collins-Akritas bisection algorithm, Sturm sequences computed as subresultant or primitive remainder sequences, and sign-change isolation for real-rooted polynomials, as used by fmpz_poly to count and isolate the real roots of squarefree integer polynomials exactly.
complexity: The Descartes bisection algorithm visits a number of nodes polynomial in the degree $n$ and the coefficient size, each costing a Taylor shift of length $n$; a Sturm sequence costs $O(n^2)$ operations on its coefficients to compute and $O(n)$ operations per evaluation.
order: 85
---

`fmpz_poly` counts the real roots of a squarefree $f \in \Z[x]$ of
degree $n$ exactly and isolates them, that is, finds disjoint intervals
with rational (here dyadic) endpoints that each contain exactly one root. Numerical
refinement of isolated roots and complex root isolation with ball
arithmetic are in `arb_fmpz_poly`.

## Root bounds {#bounds}

All isolation methods start from a bound for the roots.

::: proposition Cauchy and Fujiwara bounds {#prop-bounds}
Every complex root $\alpha$ of $f = \sum_{i=0}^n a_i x^i$ satisfies
$|\alpha| < 1 + \max_{i < n} |a_i / a_n|$ (Cauchy) and

$$
|\alpha| \le 2 \max\left( \left|\frac{a_{n-1}}{a_n}\right|,
\left|\frac{a_{n-2}}{a_n}\right|^{1/2}, \ldots,
\left|\frac{a_1}{a_n}\right|^{1/(n-1)},
\left|\frac{a_0}{2 a_n}\right|^{1/n} \right)
$$

(Fujiwara [@Fuj1916]).
:::

::: proof
For Cauchy's bound, let $M = \max_{i<n} |a_i/a_n|$ and $|\alpha| \ge 1
+ M$; then $|\sum_{i < n} a_i \alpha^i| \le M |a_n| (|\alpha|^n -
1)/(|\alpha| - 1) < |a_n| |\alpha|^n$, so $f(\alpha) \ne 0$. For
Fujiwara's bound, let $B$ be the right-hand side and $|\alpha| > B$, so
$|a_{n-k}| < |a_n| (|\alpha|/2)^k$ for $1 \le k < n$ and $|a_0| < 2
|a_n| (|\alpha|/2)^n$. Then $|\sum_{i < n} a_i \alpha^i| < |a_n|
|\alpha|^n (\sum_{k=1}^{n-1} 2^{-k} + 2^{1-n}) = |a_n| |\alpha|^n$.
:::

`fmpz_poly_bound_roots` evaluates Fujiwara's bound with integer roots
(see the warning there about the last term). For positive roots, the
*local-max* bound of Akritas, Strzeboński and Vigklas [@ASV2008],
which pairs each negative coefficient with a positive coefficient of
higher degree, is usually much sharper;
`fmpz_poly_positive_root_upper_bound_2exp` computes it as a power of
two, following Tsigaridas' SLV library.

## Descartes' rule of signs {#descartes}

::: theorem Descartes' rule of signs {#thm-descartes}
Let $V(f)$ be the number of sign changes in the sequence of nonzero
coefficients of $f$, and $P(f)$ the number of positive roots of $f$
counted with multiplicity. Then $P(f) \le V(f)$ and $V(f) - P(f)$ is
even. In particular, $V(f) = 0$ implies $P(f) = 0$, and $V(f) = 1$
implies $P(f) = 1$.
:::

The roots of $f$ in $(0, 1)$ are the positive roots of the
transformed polynomial $T(f) = (x + 1)^n f(1/(x + 1))$, so $V(T(f))$
bounds the number of roots in $(0, 1)$ with the same parity.
`_fmpz_poly_descartes_bound_0_1` computes $V(T(f))$ by a Taylor shift
of the reversal, row by row, in fixed-size limb arithmetic, and stops
as soon as the count exceeds a given bound, which is usually after a
few rows. The rule is exact when all roots are real, and it becomes
exact on small intervals by the one- and two-circle theorems of
Obreschkoff [@Obr1963], as used by Collins and Akritas [@CA1976]: for
the interval $(a, b)$, if the open disc with diameter $[a, b]$ contains
no root of $f$, then $V = 0$; if the union of the two open discs
bounded by the circumcircles of the equilateral triangles with base
$[a, b]$ contains exactly one root, then $V = 1$.

### The Vincent–Collins–Akritas algorithm {#vca}

::: algorithm Descartes bisection (VCA) {#alg-vca}
**Input:** squarefree $f$ with all positive roots in $(0, 1)$.

**Output:** isolating intervals $(c 2^{-k}, (c + 1) 2^{-k})$ and exact
roots $c 2^{-k}$ of $f$ in $(0, 1)$.

1. Start with the interval $(0, 1)$, that is, $c = 0$, $k = 0$ and
   $g = f$, where $g(x) = 2^{kn} f((x + c)/2^k)$ represents $f$ on
   $(c 2^{-k}, (c + 1) 2^{-k})$ mapped to $(0, 1)$.
2. If $g(0) = 0$, record the exact root $c 2^{-k}$ and divide $g$ by $x$.
3. Let $v = V(T(g))$. If $v = 0$, the interval has no root; if $v = 1$,
   it is an isolating interval. Otherwise bisect: replace $g$ by
   $2^n g(x/2)$ and $(c, k)$ by $(2c, k + 1)$, and go to 2.
4. After an interval with $v \le 1$, move to the next interval at the
   same or a higher level: replace $c$ by $c + 1$ and $g$ by $g(x + 1)$,
   then go up while $c$ is even, halving $c$ and replacing $g$ by
   $g(2x)$ (divided by the power of 2 of its content); stop when back
   at level 0.
:::

By the two-circle theorem, an interval whose length is small compared
to the distances between the roots gives $V \le 1$, so the bisection
terminates, and the size of the recursion tree is bounded by a
polynomial in $n$ and the bit size $\tau$ of the coefficients. The implementation is
`_fmpz_poly_isolate_real_roots_0_1_vca`; it walks the tree depth-first
with the in-place transformations of step 4, so that only one
polynomial is stored, and scales by powers of two with
`_fmpz_poly_scale_2exp`, which also removes common powers of 2. The
positive roots of a general $f$ are first mapped into $(0, 1)$ by
`_fmpz_poly_scale_positive_roots_0_1`, using the local-max bound, and
the negative roots are the positive roots of $f(-x)$. See also
Rouillier and Zimmermann [@RZ2004] for a memory-efficient formulation.

## Sturm sequences {#sturm}

::: theorem Sturm [@Stu1835] {#thm-sturm}
Let $f$ be squarefree, $p_0 = f$, $p_1 = f'$ and $p_{k+1} = -(p_{k-1}
\bmod p_k)$ until a constant $p_s$. For $a < b$ not roots of $f$, the
number of roots of $f$ in $(a, b)$ is $W(a) - W(b)$, where $W(t)$ is
the number of sign changes in $p_0(t), \ldots, p_s(t)$.
:::

The theorem remains true if each $p_k$ is multiplied by a positive
constant, so the sequence can be computed over $\Z$ by
pseudo-division, with the signs adjusted. `_fmpz_poly_num_real_roots_sturm`
uses the subresultant pseudo-remainder sequence
([](algorithms/polynomials/gcd-integers#subresultant)), whose
coefficients grow linearly along the sequence, and counts the sign
changes at $-\infty$, 0 and $+\infty$ from the leading and constant
coefficients. `_fmpz_poly_num_real_roots_sturm_bounded` uses the
primitive sequence, dividing each remainder by its content, which is
much cheaper when the remainders stay small, as for classical
orthogonal polynomials or sparse polynomials whose degrees drop
quickly, and gives up as soon as a remainder exceeds a size limit.
`fmpz_poly_signature` counts the real roots with a Sturm sequence
(Algorithm 4.1.11 of [@Coh1996]).

For isolation, the primitive sequence is stored through the
pseudo-quotients: with $\ell_k^{d} p_{k-1} = Q_k p_k + G_k p_{k+1}$,

$$
p_{k+1} = \frac{\ell_k^{d} p_{k-1} - Q_k p_k}{G_k},
$$

so the whole sequence can be evaluated at a point with $O(n)$
operations, in ball arithmetic with increasing precision and exact
arithmetic as a last resort. Bisection with these evaluations isolates
the roots without the large coefficients that the Taylor shifts of
VCA produce for polynomials with many close roots.

## Real-rooted polynomials {#signs}

If $f$ has only real roots, Descartes' rule is exact: the number of
positive roots is $V(f)$. Then it suffices to find $V(f)$ disjoint
intervals with a sign change of $f$ at their endpoints, each of which
must contain exactly one root. `fmpz_poly_isolate_real_roots` detects
the likely real-rooted case cheaply, by checking Newton's inequalities
$a_k^2 k (n - k) \ge a_{k-1} a_{k+1} (k + 1)(n - k + 1)$ (necessary for
real-rootedness) and that $V(f(x)) + V(f(-x)) = n$, and, when the
roots spread over many orders of magnitude, as for Eulerian
polynomials, it evaluates $f$ at dyadic points $2^L, \ldots, 2^K$ in
ball arithmetic and bisects until the number of sign changes equals
$V$. The result is certified by the count; if the evaluation budget is
exhausted, the function falls back to the other methods.

## Use in FLINT {#use}

- `fmpz_poly_num_real_roots` uses closed formulas for degree at most 4
  (the sign of the discriminant, and for quartics further invariants),
  then Descartes' rule on $f(x)$ and $f(-x)$, then the subresultant
  Sturm sequence for short polynomials with small coefficients, the
  bounded primitive Sturm sequence, and finally VCA.
- `fmpz_poly_num_real_roots_0_1` and `fmpz_poly_has_real_root` start
  with Descartes' rule and use the same fallbacks.
- `fmpz_poly_isolate_real_roots` and `fmpz_poly_isolate_positive_roots`
  try Sturm bisection for large degree, then sign changes for
  real-rooted polynomials, then VCA.
- `fmpz_poly_bound_roots`, `fmpz_poly_positive_root_upper_bound_2exp`
  and `_fmpz_poly_descartes_bound_0_1` provide the bounds.
- `arb_fmpz_poly_real_roots` and `arb_fmpz_poly_complex_roots` compute
  enclosures of the roots to any precision; `fmpz_poly_roots_fmpq`
  finds the rational roots exactly.

::: functions fmpz_poly_num_real_roots fmpz_poly_isolate_real_roots fmpz_poly_num_real_roots_sturm _fmpz_poly_num_real_roots_sturm_bounded _fmpz_poly_isolate_real_roots_0_1_vca _fmpz_poly_descartes_bound_0_1 fmpz_poly_positive_root_upper_bound_2exp fmpz_poly_signature :::
