---
title: Algebraic numbers by minimal polynomial and isolating ball
summary: The canonical representation of an algebraic number by its minimal polynomial over Z and a complex ball isolating one root, its certification by the interval Newton method, refinement, and exact decisions of equality, realness and signs, as implemented in qqbar.
complexity: Equality, sign and comparison tests cost one or a few interval Newton refinements with $O(d)$ ball operations each at a precision proportional to the separation of the numbers; creating the numbers costs a factorisation over $\Z$ and a complex root isolation.
order: 10
---

An algebraic number $x \in \overline{\mathbb{Q}}$ is a root of a
nonzero polynomial with integer coefficients. Among these polynomials
there is a unique one, the *minimal polynomial* $P$ of $x$, that is
irreducible in $\Z[t]$, primitive (its coefficients have no common
factor) and has a positive leading coefficient; every other integer
polynomial vanishing at $x$ is a multiple of $P$. The degree of $P$ is
the *degree* of $x$, and the roots of $P$ are the *conjugates* of $x$.
Since $P$ alone does not distinguish $x$ from its conjugates, [](qqbar)
stores a pair

$$
x \;\leftrightarrow\; (P, B),
$$

where $B$ is a complex ball, a rectangle with real and imaginary parts
given by [](arb) balls, that contains $x$ and no other root of $P$ (an
*isolating* ball). This representation, used for example by Loos
[@Loo1983] and Strzebonski [@Str1997], is canonical: equal numbers have
the same $P$, and different balls for the same number are recognised by
refining them. Arithmetic is described in
[](algorithms/exact/qqbar-arithmetic).

## Certifying an isolating ball {#newton}

Isolating balls are certified with the interval Newton method [@Moo1966;
@Kra1969; @Neu1990]. For a ball $B$ with center $m$ and a polynomial $P$,
let $P'(B)$ be a ball containing $\{P'(z) : z \in B\}$ (computed in ball
arithmetic) and define the Newton operator

$$
N(B) = m - \frac{P(m)}{P'(B)}.
$$

::: theorem Interval Newton {#thm-newton}
If $0 \notin P'(B)$ and $N(B)$ is contained in the interior of $B$,
then $B$ contains exactly one root of $P$, and this root lies in
$N(B)$.
:::

::: proof
For $z \in B$, write $P(z) = P(m) + (z - m) s(z)$ with $s(z) =
\int_0^1 P'(m + \tau (z - m)) \, d\tau$. The segment from $m$ to $z$
lies in the convex set $B$, so $s(z)$ is an average of values of $P'$ on
$B$ and lies in the convex hull of $\{P'(w) : w \in B\}$, which is
contained in the ball (a rectangle) $P'(B)$. In particular $s(z) \ne 0$.

*Uniqueness.* If $z_1 \ne z_2$ are roots of $P$ in $B$, then $0 =
P(z_1) - P(z_2) = (z_1 - z_2) \int_0^1 P'(z_2 + \tau(z_1 - z_2)) \,
d\tau$, and the integral lies in $P'(B)$, which would have to contain
0.

*Existence.* The map $T(z) = m - P(m)/s(z)$ is continuous on $B$, since
$s$ is continuous and nonzero, and maps $B$ into $N(B) \subseteq B$. By
Brouwer's fixed point theorem it has a fixed point $z^* = m -
P(m)/s(z^*)$, and then $P(z^*) = P(m) + (z^* - m) s(z^*) = 0$. The
root $z^*$ lies in $T(B) \subseteq N(B)$.
:::

For a real interval and a real polynomial the same argument, with the
mean value theorem on an interval, shows that $N(B) \subset
\operatorname{int} B$ certifies a unique *real* root in $B$. FLINT uses
both forms:

- `_qqbar_validate_existence_uniqueness` applies [](#thm-newton) to an
  arbitrary ball; when the ball is real (or purely imaginary) the step is
  carried out on the real line (or the imaginary axis), which proves that
  the unique root lies there.
- `_qqbar_validate_uniqueness` is used when the ball is already known to
  contain a root, for example a ball computed by a root-finding
  algorithm or an enclosure of $f(x)$ obtained by ball arithmetic from
  an enclosure of $x$. It slightly enlarges the ball before the Newton
  step, which makes the test succeed more easily, and returns the
  (usually much smaller) ball $N(B)$ as the new enclosure.

## Creating algebraic numbers {#creating}

Rational numbers $p/q$ have $P = q t - p$ and an exact enclosure. The
roots of a polynomial $f \in \Z[t]$ are created by factoring $f$ over
$\Z$ (see [](algorithms/polynomials/factoring-integers)) and isolating
the complex roots of each irreducible factor $P$ numerically with
`arb_fmpz_poly_complex_roots`, increasing the precision until all $\deg
P$ balls are certified by [](#thm-newton) (see [](qqbar/roots)). Values
of special functions at algebraic points that are algebraic, such as
$e^{2 \pi i p/q}$ and $\cos(\pi p/q)$, get their minimal polynomials from
formulas: cyclotomic polynomials $\Phi_n$, and the minimal polynomial of
$\cos(2 \pi / n)$ (`fmpz_poly_cos_minpoly`), and their enclosures from
[](arb) evaluations of the function, which are then certified.

Every operation of [](qqbar) produces its result in this form: it
determines the minimal polynomial exactly and an isolating ball for the
right root.

## Refining an enclosure {#refining}

To obtain more precision, `_qqbar_enclosure_raw` iterates the Newton
step $B \gets N(B)$, doubling the working precision each time; by
[](#thm-newton) every $N(B)$ again contains the root, and the
convergence is quadratic once $B$ is small compared to the distance to
the other roots. If a step fails to shrink the ball sufficiently (for
example because the starting ball is too wide for the Newton iteration
to contract), all roots are recomputed from scratch at the higher
precision and the one overlapping the old ball is taken. A real
enclosure is refined with the real version (`arb_fmpz_poly_refine_root_arb`
on the real line).

The *polished* enclosures returned by `qqbar_get_acb` are refined until
both the real and the imaginary part have the requested relative
accuracy; a part that is known to be zero is set to an exact zero, and
a part that equals a small dyadic number $n 2^{-k}$ is detected by
testing $x - n 2^{-k}$ exactly, so that, for example, $\frac{1}{2} +
\frac{\sqrt{3}}{2} i$ has an exact real part $\frac{1}{2}$.
`qqbar_cache_enclosure` stores a refined ball in the number itself,
which helps when the same number is used many times.

## Deciding predicates {#predicates}

All predicates are decided exactly by combining refinement with simple
algebraic facts.

::: lemma Equality {#lem-equal}
Let $(P_1, B_1)$ and $(P_2, B_2)$ represent $x_1$ and $x_2$. Then $x_1 =
x_2$ if and only if $P_1 = P_2$ and some ball containing $B_1 \cup B_2$
contains exactly one root of $P_1$.
:::

::: proof
If $x_1 = x_2$, the minimal polynomials agree by uniqueness, and as the
balls are refined towards the common root their union shrinks to it, so
eventually a ball containing it isolates it. Conversely, a ball
containing $B_1 \cup B_2$ that contains a single root of $P_1 = P_2$
must contain both $x_1$ and $x_2$, so they are that root.
:::

`qqbar_equal` therefore compares the minimal polynomials first (which
decides most cases at once), then the balls, and refines both balls
until they are disjoint (the numbers differ) or the union of the balls
is certified by `_qqbar_validate_uniqueness` (they are equal). Since
distinct numbers have disjoint balls at sufficient precision, this
terminates.

::: lemma Real and imaginary numbers {#lem-real}
Let $(P, B)$ represent $x$ with $\deg P \ge 2$.

1. $x$ is real if and only if the real interval $B \cap \mathbb{R}$
   contains a root of $P$.
2. If $\operatorname{Re}(x) = 0$, then $P$ is even, that is, $P(-t) =
   P(t)$.
:::

::: proof
1. If $x$ is real, $x \in B \cap \mathbb{R}$. Conversely, a root of $P$
   in $B \cap \mathbb{R} \subseteq B$ is the unique root in $B$, which is
   $x$.
2. If $x = iy$ with $y \ne 0$ real, then $-x = \overline{x}$ is also a
   root of the real polynomial $P$. So $P(t)$ and $P(-t)$ have the common
   root $x$; since $P$ is irreducible, $P(-t) = c P(t)$, and comparing
   leading coefficients gives $c = \pm 1$. If $c = -1$, then $P(0) = 0$,
   so $t$ divides $P$, which contradicts $\deg P \ge 2$ and
   irreducibility.
:::

`qqbar_sgn_im` reads the sign of $\operatorname{Im}(x)$ off the ball if
the imaginary part of the ball excludes zero or is exactly zero.
Otherwise it alternates refinement with an attempt to certify, with the
real version of the Newton step, that the real interval obtained by
setting the imaginary part of the ball to zero contains a root; by
[](#lem-real) this proves that $x$ is real. If $x$ is not real, the
refinement eventually separates the ball from the real axis. In the same
way `qqbar_sgn_re` decides the sign of the real part, attempting the
certification on the imaginary axis only when $P$ is even.
$\operatorname{Re}(x) = 0$ and $x \ne 0$ is impossible for odd degree,
and a real number of degree at least 2 is never zero.

The comparisons `qqbar_cmp_re`, `qqbar_cmp_im` and `qqbar_cmpabs`
compare the corresponding parts of the balls, refining them while they
overlap. Overlapping balls may represent equal parts of distinct
numbers (for example complex conjugates have equal real parts), so if
refinement does not separate the parts quickly, the difference is
computed exactly with [](qqbar) arithmetic (as $x - y$, or as $x -
\overline{y}$ for conjugates) and its sign is determined as above. Exact
tests on the minimal polynomial (degree, evenness, cyclotomicity) are
used first wherever possible: `qqbar_is_root_of_unity` only checks
whether $P$ is a cyclotomic polynomial with `fmpz_poly_is_cyclotomic`,
and then reads off the index from a low-precision value of $\arg x$.

## Size of the representation {#size}

The minimal polynomial of a number of degree $d$ has $d + 1$
coefficients, and for numbers obtained by arithmetic their size can grow
quickly: a sum $x + y$ of numbers of degrees $d_1$, $d_2$ and heights
$H_1$, $H_2$ generically has degree $d_1 d_2$ and height roughly
$H_1^{d_2} H_2^{d_1}$. The enclosure needs a precision related to the
separation of the roots of $P$, which by Mahler's bound [@Mah1964] is at worst
exponentially small in $d \log(d H)$, although typically it is
moderate. The representation is therefore best suited to individual
numbers of moderate degree; long computations in a fixed number field
are better done in that field ([](nf_elem), or [](ca), which uses
`qqbar_t` only to describe the generator).

## Use in FLINT {#use}

- `qqbar_equal`, `qqbar_sgn_re`, `qqbar_sgn_im`, `qqbar_is_real`,
  `qqbar_cmp_re`, `qqbar_cmp_im`, `qqbar_cmpabs` and
  `qqbar_cmp_root_order` decide predicates as described above.
- `_qqbar_validate_uniqueness`, `_qqbar_validate_existence_uniqueness`,
  `_qqbar_enclosure_raw` and `qqbar_cache_enclosure` implement the
  interval Newton certification and refinement; `qqbar_get_acb` and
  `qqbar_get_arb` polish enclosures.
- `qqbar_roots_fmpz_poly`, `qqbar_real_roots_fmpz_poly`,
  `qqbar_root_of_unity`, `qqbar_cos_pi` and related functions create
  numbers from polynomials and from formulas.
- The generic contexts `gr_ctx_init_real_qqbar` and
  `gr_ctx_init_complex_qqbar` (see [](gr_domains)) use this
  representation, and [](ca) uses it for the generators of number
  fields and to decide predicates for algebraic numbers.
