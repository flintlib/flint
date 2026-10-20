---
title: Power series in ball arithmetic
summary: How truncated power series with ball coefficients give rigorous Taylor expansions, why their correctness is automatic, and how the choice between recurrences and Newton iteration affects the width of the enclosures; with the algorithms used by arb_poly, acb_poly and the generic gr_poly code for inversion, exponentials, roots, logarithms, trigonometric functions and composition.
complexity: $O(M(n))$ operations for inversion, division, square roots, exponentials, logarithms and trigonometric functions of a series of length $n$ by Newton iteration; $O(n \cdot m)$ for the recurrences, where $m$ is the length of the input; composition $O(M(n) \log n)$ to $O(\sqrt{n} M(n))$.
order: 50
---

A power series $f = \sum_k f_k x^k$ truncated to length $n$ is stored as
the polynomial $f_0 + \cdots + f_{n-1} x^{n-1}$, and an operation on
series computes the first $n$ coefficients of its result. Since the
first $n$ coefficients of $f + g$, $fg$, $1/f$, $\exp f$, $F(f)$ for an
analytic function $F$, ... depend only on the first $n$ coefficients of
the inputs, truncated series arithmetic is exact arithmetic in the ring
$R[x]/(x^n)$. With ball coefficients, every algorithm that is correct in
exact arithmetic therefore becomes an algorithm that returns
enclosures, simply by performing each ring operation in ball
arithmetic. Rigor is automatic; what the choice of algorithm decides is
the speed and the *width* of the enclosures. The generic algorithms (Newton iteration, middle
products, the formulas for the elementary functions) are described in
[](algorithms/polynomials/newton); the focus here is on what changes
with balls.

## Rigorous Taylor expansions {#taylor}

If $h$ is the series $x_0 + x$, then $F(h) \bmod x^n$ is the Taylor
polynomial $\sum_{k<n} F^{(k)}(x_0) x^k / k!$, and composing series
operations yields the Taylor polynomial of a composite function (the
chain rule is built into series composition). In ball arithmetic $x_0$
may be a ball, and then each computed coefficient contains $F^{(k)}(t)/k!$
for every $t \in x_0$. This is forward-mode automatic differentiation of
arbitrary order with rigorous error bounds, and it is how FLINT computes
derivatives of special functions, Taylor models for integration and root
finding ([](algorithms/ball-arithmetic/real-root-isolation)), and the
Taylor coefficients used by `arb_calc_func_t` and `acb_calc_func_t`
evaluators.

For a transcendental function $F$, the functions in [](arb_poly) and
[](acb_poly) usually compute the series of $F$ at the constant term
$h_0$ of the argument, $F(h_0 + x) = \sum_k F_k x^k$, by a method
specific to $F$, and then compose it with the nonconstant part of $h$:
$$F(h) = \sum_k F_k (h - h_0)^k \bmod x^n,$$
using `_arb_poly_compose_series` (or simpler formulas when $h - h_0$ is
a monomial). The gamma function (Stirling's series at a shifted
argument, or zeta values at integers), the zeta function
(Euler–Maclaurin summation), the arithmetic-geometric mean and the
sinc function at a ball containing 0 are computed this way. The Lambert
W function, the elementary functions and the inverse functions use
Newton iteration or differential equations directly on the series.

## Recurrences and Newton iteration {#newton}

Many series operations can be computed either by a recurrence for the
coefficients or by Newton iteration. For the inverse $q = 1/a$, the
recurrence is
$$q_0 = 1/a_0, \qquad q_k = -q_0 \sum_{j=1}^{k} a_j q_{k-j},$$
and Newton iteration doubles the number of correct coefficients with
$$q \leftarrow q - q(aq - 1),$$
using that the first $m$ coefficients of $aq - 1$ vanish when $q$ is
correct to $m$ terms, so that only a middle product of length $m$ is
needed ([](algorithms/polynomials/middle-product)). For the exponential
$f = \exp h$, the recurrence follows from the differential equation $f'
= h' f$:
$$k f_k = \sum_{j=1}^{k} j h_j f_{k-j},$$
and Newton iteration solves $\log f = h$, updating $f$ and $1/f$
together as in [@HZ2004]. The recurrences cost $O(nm)$ operations for an
input of length $m$; Newton iteration costs $O(M(n))$.

In exact arithmetic, the choice is purely a matter of speed. In ball
arithmetic, the two methods also produce enclosures of different width.
The recurrence for $q_k$ uses all previous $q_j$, and each $q_j$ already
carries a radius that includes the rounding errors and the propagated
errors of all earlier coefficients. Ball arithmetic cannot exploit the
cancellation between these correlated errors: a rounding error in $q_j$
is propagated to $q_k$ with the weight of the coefficient of $x^{k-j}$
in the majorant series $1/(1 - |q_0| \sum_{j \ge 1} |a_j| x^j)$, whose
coefficients can grow exponentially faster than those of $1/a$ (for
example, when the $a_j$ alternate in sign). Newton iteration computes the second half of
the coefficients from a product involving only the first half, through
a much shorter chain of dependencies, and the radius growth is
correspondingly smaller. FLINT therefore uses the recurrences only for
very short series or special shapes:

- `_arb_poly_inv_series` uses the recurrence for all coefficients only
  if $n \le 8$ or the input has length 2, and otherwise only for the
  first 4 coefficients, followed by Newton iteration, although the
  recurrence would be faster up to much larger lengths;
- `_arb_poly_div_series` uses the recurrence for $n \le 10$ (or a
  divisor of length 2) and the product $A \cdot (1/B)$ otherwise;
- `_arb_poly_exp_series` uses the recurrence up to an input length of
  750 at $p \le 256$ bits (or $10^5/(\ln p)^3$), and Newton iteration
  beyond.

The same reasoning applies to the square root (Newton iteration for
$h^{-1/2}$ without divisions, $g \leftarrow g + g(1 - hg^2)/2$, followed
by one Karp–Markstein step [@KM1997] that gives $\sqrt{h}$ from $h^{-1/2}$
at half the length), and to the trigonometric functions: for long
series, `_gr_poly_sin_cos_series` uses a simultaneous Newton iteration
for $\sin h$ and $\cos h$ derived from that for $e^{ih}$, and
`_gr_poly_tan_series` Newton iteration for the inverse of
$\operatorname{atan}$; the inverse functions are integrals of algebraic
expressions in $h$ and $h'$, such as $\operatorname{atan} h =
\operatorname{atan} h_0 + \int h'/(1 + h^2)$, and the logarithm is
$\log h = \log h_0 + \int h'/h$.

## Scaling {#scaling}

The coefficients of the series that arise in practice often decrease
like $\rho^{-k}$, where $\rho$ is the radius of convergence. Ball
arithmetic measures errors relative to the size of each coefficient, so
the relative accuracy of the coefficients of a product is governed by
the largest terms in each convolution sum: if the coefficients of both
factors decay at the same rate, all products $a_j b_{k-j}$ have the same
order of magnitude, and every coefficient of the result is obtained to
nearly full relative precision. The block multiplication algorithm
([](algorithms/ball-arithmetic/polynomial-multiplication)) preserves
this property by rescaling $x$. If the series decay at different rates,
the radii of the small coefficients are dominated by errors of larger
terms, and it can help to rescale the variable explicitly, computing
$F(h(cx))$ for a suitable $c$ and undoing the scaling at the end.

## Composition and reversion {#composition}

Composition $f(g) \bmod x^n$ with $g_0 = 0$ and reversion (the series of
the compositional inverse) are computed by the generic functions
`_gr_poly_compose_series` and `_gr_poly_revert_series`
([](algorithms/polynomials/composition)). Over the ball contexts they
choose Horner's scheme for short input, divide and conquer when the
composition is only slightly longer than $n$, the baby-step giant-step
algorithm of Brent and Kung [@BrentKung1978] with $O(\sqrt{n} M(n) +
n^{(\omega+1)/2})$ operations below a precision-dependent cutoff
$\max(100, 15000/\sqrt{\max(64, p)})$, and the $O(M(n) \log n)$
algorithm of Kinoshita and Li [@KL2024] above it. Reversion uses the
fast Lagrange inversion of [@Joh2015b] below five times that cutoff and
Newton iteration on the composition above.

## Use in FLINT {#use}

- Power series functions of [](arb_poly) and [](acb_poly) implemented in
  the modules themselves: `_arb_poly_inv_series`,
  `_arb_poly_div_series`, `_arb_poly_exp_series` (with
  `_arb_poly_exp_series_basecase`), `_arb_poly_sqrt_series`,
  `_arb_poly_pow_arb_series`, `_arb_poly_sinh_cosh_series`,
  `_arb_poly_sinc_series`, `_arb_poly_lambertw_series`,
  `_arb_poly_gamma_series` and its relatives, `_arb_poly_zeta_series`,
  and the complex versions.
- Wrappers of [](gr_poly) over the ball contexts: the logarithm, the
  inverse trigonometric and hyperbolic functions, sine and cosine (for
  length at least 3), tangent and cotangent, the reciprocal square root
  (`_gr_poly_rsqrt_series_newton`), binary powering, composition,
  reversion, Taylor shifts, derivatives and integrals.
- Middle products for Newton iteration: `_arb_poly_mulmid`,
  `_acb_poly_mulmid`.
- Users of the Taylor expansion idiom: [](arb_calc), [](acb_calc), the
  hypergeometric functions of [](arb_hypgeom) and [](acb_hypgeom) with
  series arguments, and `acb_dirichlet_zeta_jet`.
