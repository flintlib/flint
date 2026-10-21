---
title: Certified roots of polynomials
summary: How acb_poly_find_roots computes all complex roots of a polynomial by the Weierstrass–Durand–Kerner iteration from Newton polygon starting values and certifies them with inclusion disks; how real roots are recognised; root bounds; and how arb_fmpz_poly combines this with exact real root isolation.
complexity: $O(n^2)$ operations at $p$ bits per iteration for degree $n$, with $O(n)$ iterations typically; validation $O(n^2)$ operations; isolation needs a precision of the order of $n$ times the logarithm of the inverse root separation in the worst case.
order: 60
---

Given a polynomial $f(x) = \sum_{k=0}^{n} a_k x^k$ of degree $n$ with
complex ball coefficients, `acb_poly_find_roots` computes balls $B_1,
\ldots, B_n$ and a number $s$ such that each of $B_1, \ldots, B_s$
contains exactly one root of every polynomial in the coefficient balls,
and each of the remaining balls contains at least one root. If $s = n$,
all roots are isolated. The computation separates the numerical part,
which approximates the roots with floating-point arithmetic and no error
control, from the certification, which is a short computation in ball
arithmetic based on the following classical bound.

## Inclusion disks {#inclusion}

::: theorem Inclusion disk {#thm-inclusion}
Let $f$ be a polynomial of degree $n \ge 1$, $m \in \C$ with $f'(m) \ne
0$, and $\rho = n |f(m)/f'(m)|$. Then the closed disk of radius $\rho$
around $m$ contains a root of $f$.
:::

::: proof
Let $\zeta_1, \ldots, \zeta_n$ be the roots of $f$. If $f(m) = 0$ there
is nothing to prove. Otherwise $f'(m)/f(m) = \sum_i 1/(m - \zeta_i)$. If
$|m - \zeta_i| > \rho$ for all $i$, then
$$\left|\frac{f'(m)}{f(m)}\right| \le \sum_{i=1}^{n} \frac{1}{|m - \zeta_i|} < \frac{n}{\rho} = \left|\frac{f'(m)}{f(m)}\right|,$$
a contradiction [@Kob2010].
:::

::: corollary Isolation {#cor-isolation}
If disks $D_1, \ldots, D_n$ around points $m_1, \ldots, m_n$ with the
radii of [](#thm-inclusion) are pairwise disjoint, then each contains
exactly one root of $f$, and these are all the roots.
:::

::: proof
Each disk contains at least one root, the disks are disjoint, and $f$
has exactly $n$ roots counted with multiplicity.
:::

In ball arithmetic, `_acb_poly_root_inclusion` evaluates $f(m)$ and
$f'(m)$ at the exact point $m$ (the midpoint of an approximation) with
the coefficient balls, takes an upper bound $u \ge |f(m)|$ and an upper
bound $v \ge 1/|f'(m)|$ over all polynomials in the balls, and returns
the box around the disk of radius $n u v$. Since $u v$ bounds
$|f(m)/f'(m)|$ for every polynomial in the coefficient balls, the box
contains a root of each of them, and [](#cor-isolation) applies to each
of them. If $f'(m)$ contains zero, $v$ is infinite and the box is the
whole plane. The evaluation uses Horner's scheme with the error of each
partial result bounded by a disk (the larger of the real and imaginary
radii would overestimate it by a factor $\sqrt 2$ per step, which
matters for high degree).

`_acb_poly_validate_roots` computes these boxes for all approximations,
determines which boxes overlap no other box, moves them to the front
and returns their number $s$; by [](#thm-inclusion) and the argument of
[](#cor-isolation), each of these $s$ boxes contains exactly one root
if $s = n$, and in general each of them contains at least one root and
the boxes that overlap no other one are isolated from each other, but
two roots may share an overlapping group.

The radius $n|f(m)/f'(m)|$ is about $n$ times the distance of $m$ to
the nearest root when that root is simple and well separated, so the
enclosure is only $\log_2 n$ bits less accurate than the approximation.
At a root of multiplicity $k$, $|f/f'| \approx |m - \zeta|/k$, but $m$ is
only accurate to about $p/k$ bits, so the disks of the $k$
approximations overlap at every precision: multiple roots can never be
isolated, and clusters of roots need a precision that resolves them.

## The Weierstrass–Durand–Kerner iteration {#wdk}

Given approximations $z_1, \ldots, z_n$ of the roots of $f$, the
iteration
$$z_k \leftarrow z_k - \frac{f(z_k)}{a_n \prod_{j \ne k} (z_k - z_j)}$$
of Weierstrass [@Wei1891], Durand [@Dur1960] and Kerner [@Ker1966] is
Newton's method for the system of equations that expresses the
coefficients of $a_n \prod_k (x - z_k)$ in terms of the $z_k$ (Vieta's
formulas). It converges quadratically to simple roots from good starting
values, and in practice from almost all starting values on a circle
around the roots. The Aberth iteration [@Abe1973], which converges
cubically, is implemented in `_gr_poly_refine_roots_aberth` but not
used by `acb_poly_find_roots`. Each iteration costs $O(n^2)$
operations.

The starting values are chosen by the Newton polygon technique of Bini
[@Bin1996]. Let $k_1 < k_2 < \cdots < k_q$ be the indices of the
vertices of the upper convex hull of the points $(k, \log_2 |a_k|)$.
For each edge, $m_i = k_i - k_{i-1}$ points are placed on the circle of
radius $|a_{k_{i-1}}/a_{k_i}|^{1/m_i}$, which is where the polynomial
has $m_i$ roots of approximately that modulus when the coefficients
dominate as the polygon indicates. The angles are spread with an
irrational offset so that the configuration has no symmetry that the
iteration would preserve.

## Implementation {#implementation}

`_acb_poly_find_roots` proceeds as follows for degree $n$:

1. If the leading coefficient contains zero, the roots are unbounded:
   all outputs are set to $[\pm\infty] + [\pm\infty]i$ and 0 is
   returned. Degree 1 is solved directly.
2. If $p \le 53$, the polynomial is first rescaled so that its
   coefficients fit in double precision: with the line $y = \alpha k +
   \beta$ that best fits the points $(k, \log_2 |a_k|)$ in the minimax
   sense (computed from the upper convex hull and the endpoints, a
   rotating calipers argument), it substitutes $x \to 2^{-\alpha} x$ and
   divides by $2^{\beta}$. `cd_poly_find_roots` then runs the
   iteration in double precision, with steps limited in size (between
   $2^{-10}$ and $2^{-3}$, adapted to progress) and roots of modulus
   above 1 updated through the reversed polynomial to avoid overflow. If
   the last relative correction is below $1/8$, the result goes to
   validation.
3. Otherwise the iteration runs at $p$ bits with floating-point complex
   numbers (`nfloat_complex` if available, `acf` otherwise) through the
   generic `_gr_poly_refine_roots_wdk`, for at most $2n + \sqrt{p}$
   iterations by default. When the largest correction, relative to the
   largest root, falls below $2^{-p/2}$ (respectively $2^{-p/3}$,
   $2^{-p/4}$), the iteration is stopped two (three, four) iterations
   later, since quadratic convergence then reaches full precision.
4. The approximations are validated by `_acb_poly_validate_roots`.

The working precision must exceed the number of bits lost in evaluating
$f$ near the roots, which for a root $\zeta$ is about $\log_2$ of the
condition number $\sum_k |a_k| |\zeta|^k / |f'(\zeta)|$. Users therefore
call the function in a loop with doubling precision, reusing the
previous approximations as starting values. For polynomials with exact
coefficients, the radius of the inclusion disks shrinks proportionally
to $2^{-p}$ once the precision suffices, and every squarefree
polynomial is eventually isolated.

## Real roots {#real}

For a real polynomial, the roots in the upper and the lower half-plane
come in conjugate pairs, and a root whose isolating box meets the real
axis might still be nonreal. `acb_poly_validate_real_roots` proves
that it is real: let $r$ be the number of boxes that meet the real axis.
The test requires $r \equiv n \pmod 2$ and checks, at the midpoints
between the consecutive candidate intervals (sorted by midpoint), that
$f$ has signs that alternate starting from the sign at $-\infty$ and
ending with the sign at $+\infty$. Then $f$ changes sign within each
candidate interval, which by the intermediate value theorem contains a
real root; since each box contains exactly one root, that root is real,
and the other boxes contain the nonreal roots.

## Root bounds {#bounds}

`arb_poly_root_bound_fujiwara` and `acb_poly_root_bound_fujiwara` bound
the absolute values of all roots by Fujiwara's bound [@Fuj1916]
$$|\zeta| \le 2 \max\left\{\left|\frac{a_{n-1}}{a_n}\right|, \left|\frac{a_{n-2}}{a_n}\right|^{1/2}, \ldots, \left|\frac{a_1}{a_n}\right|^{1/(n-1)}, \left|\frac{a_0}{2 a_n}\right|^{1/n}\right\},$$
computed with upper bounds for the numerators and a lower bound for
$|a_n|$, so that it holds for every polynomial in the balls. The bound
is within a factor $2n$ of the largest root modulus.

## Integer polynomials {#integer}

For a squarefree integer polynomial, `arb_fmpz_poly_complex_roots`
always succeeds, because it controls the precision itself:

1. It writes $f(x) = x^e g(x^d)$ with the largest $d$, so that the roots
   of $f$ are 0 (if $e = 1$) and the $d$-th roots of the roots of $g$.
2. It isolates the real roots of $g$ exactly by Descartes' rule of
   signs ([](algorithms/polynomials/real-roots-integers)) and refines
   them by interval Newton iteration (`arb_fmpz_poly_refine_root_arb`,
   see [](algorithms/ball-arithmetic/real-root-isolation)).
3. If the real roots are numerous or clustered, it divides them out in
   ball arithmetic, increasing the precision until the quotient (whose
   roots are the nonreal roots) is accurate to the working precision, so
   that the Weierstrass–Durand–Kerner iteration does not have to
   resolve the clusters.
4. It calls `acb_poly_find_roots` at 53, 128, 256, ... bits until all
   roots are isolated and accurate to the target, proves which roots are
   real with `acb_poly_validate_real_roots` when necessary, and takes
   $d$-th roots.

The output lists the real roots in increasing order, then the nonreal
roots in conjugate pairs; the second member of each pair is the
conjugate of the first, so that the pairs are exactly symmetric.

## Use in FLINT {#use}

- `_acb_poly_root_inclusion` ([](#thm-inclusion)),
  `_acb_poly_validate_roots`, `_acb_poly_find_roots` and
  `acb_poly_find_roots`, `_acb_poly_find_roots_double` and the
  `cd_poly` functions for double precision, `acb_poly_validate_real_roots`,
  and the root bounds `arb_poly_root_bound_fujiwara` and
  `acb_poly_root_bound_fujiwara`.
- The generic iterations `_gr_poly_refine_roots_wdk` and
  `_gr_poly_refine_roots_aberth` in [](gr_poly).
- `arb_fmpz_poly_complex_roots` and `arb_fmpz_poly_real_roots`
  ([](arb_fmpz_poly/roots)); `qqbar` uses them to compute algebraic
  numbers ([](qqbar)).
- The generic root finding over the complex ball context
  (`gr_poly_roots` with `_gr_acb_poly_roots`) and factorisation into
  linear factors (`gr_generic_poly_factor_roots`) call
  `acb_poly_find_roots` with doubling precision.
