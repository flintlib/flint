---
title: The dilogarithm
summary: Evaluation of the dilogarithm Li₂(z) on the whole complex plane — the series at 0 after one of four transformations, a series in log z with Bernoulli numbers near exp(±πi/3), Taylor continuation with explicit coefficient bounds, and the bit-burst algorithm at very high precision.
complexity: O(M(p) log² p) at p bits with the bit-burst algorithm; O(p^{1/2}) full multiplications with rectangular splitting at moderate precision.
order: 49
---

The dilogarithm
$$
\operatorname{Li}_2(z) = -\int_0^z \frac{\log(1-t)}{t}\,dt =
\sum_{k=1}^\infty \frac{z^k}{k^2} = z\,{}_3F_2(1,1,1;2,2;z)
$$
is analytic on $\mathbb{C} \setminus [1, \infty)$. `acb_hypgeom_dilog`
evaluates it for complex balls; it is a holonomic function, and the methods below are
instances of the general ones in
[](algorithms/special-functions/hypergeometric). General polylogarithms
$\operatorname{Li}_s(z)$ are computed in [`acb_dirichlet`](acb_dirichlet)
and [`acb_poly`](acb_poly) by other methods
([](algorithms/special-functions/zeta)).

## Transformations {#transformations}

The functional equations
$$
\operatorname{Li}_2(z) = -\operatorname{Li}_2(1/z) - \tfrac12\log^2(-z) - \tfrac{\pi^2}{6}
\quad (z \notin (0,1)),
$$
$$
\operatorname{Li}_2(z) = -\operatorname{Li}_2(1-z) - \log(1-z)\log(z) + \tfrac{\pi^2}{6},
$$
$$
\operatorname{Li}_2(z) = -\operatorname{Li}_2\left(\tfrac{z}{z-1}\right) - \tfrac12\log^2(1-z)
\quad (z \notin (1, \infty)),
$$
$$
\operatorname{Li}_2(z) = \operatorname{Li}_2\left(\tfrac{1}{1-z}\right) +
\log(1-z)\left[\tfrac12\log(1-z) - \log(-z)\right] - \tfrac{\pi^2}{6}
$$
map $z$ to $1/z$, $1 - z$, $z/(z-1)$ and $1/(1-z)$
(`acb_hypgeom_dilog_transform` with `algorithm` 1 to 4; each returns an
indeterminate result if the ball meets the region where the formula is
invalid). As for ${}_2F_1$, the union of the disks $|w| \le r$ over the
transformed arguments $w$ covers the plane except neighbourhoods of the
points $e^{\pm\pi i/3}$, which are fixed points of the transformation
group.

`acb_hypgeom_dilog` computes $|z|^2$, $|1/z|^2$, $|1-z|^2$,
$|z/(z-1)|^2$ and $|1/(1-z)|^2$ in double precision and uses the series
at 0 directly if $|z| \le 1/4$, and otherwise after the transformation
that minimizes $|w|$. If the best $|w|$ is still above $1/4$, it
considers two alternatives:

- below 10000 bits, the series in $\log z$ described next, if its rate
  of convergence $|\log z/(2\pi)|$ (with a penalty factor) is better;
- if the best $|w|$ exceeds $3/4$ (which happens only near
  $e^{\pm\pi i/3}$ when the Bernoulli series was not chosen), expansion
  at one of the points $i$, $(1+i)/2$, $1+i$ (or their conjugates),
  where $\operatorname{Li}_2$ has a closed form, for example
  $\operatorname{Li}_2(i) = -\pi^2/48 + iG$ with Catalan's constant
  $G$, followed by bit-burst continuation to $z$.

The working precision is limited to the accuracy of the input plus 30
bits, and raised by 0.5% plus 5 bits.

## The series at zero {#series}

`acb_hypgeom_dilog_zero_taylor` sums $\sum_{k \ge 1} z^k/k^2$ for small
$|z|$. Since the coefficients are reciprocals of small integers, the
sum is computed by rectangular splitting in powers of $z$, where
several consecutive denominators $k^2$ are combined into one word-size
divisor, or by binary splitting when $z$ has few bits; the tail after
$n$ terms is bounded by $|z|^n/(1-|z|)$. `acb_hypgeom_dilog_zero` uses
this below 40000 bits or for tiny $|z|$, and the bit-burst algorithm
otherwise.

## The series in log z {#bernoulli}

With $w = \log z$, for $|w| < 2\pi$,
$$
\operatorname{Li}_2(e^w) = \frac{\pi^2}{6} + w(1 - \log(-w)) -
\frac{w^2}{4} - \sum_{k \ge 3,\ k \text{ odd}} \frac{B_{k-1}}{(k-1)\,k!} w^k,
$$
which follows from $\operatorname{Li}_2(e^w) = \zeta(2) + w(1 -
\log(-w)) + \sum_{k \ge 2} \zeta(2-k) w^k/k!$ and $\zeta(2-k) =
-B_{k-1}/(k-1)$. For $k \ge 4$ the terms are bounded by $(|w|/(2\pi))^k$,
which gives the truncation bound as a geometric series.
`acb_hypgeom_dilog_bernoulli` evaluates this with Horner's scheme in
$w^2$ and the cached Bernoulli numbers. Near $z = e^{\pm\pi i/3}$, $|w|
\approx \pi/3$ and each term gains about 2.6 bits. Since the Bernoulli
numbers grow in size, the method is used only below 10000 bits.

## Taylor continuation {#continuation}

`acb_hypgeom_dilog_continuation` computes $\operatorname{Li}_2(z) -
\operatorname{Li}_2(a)$ from the Taylor expansion at $a$. Since
$\frac{d}{dz}\operatorname{Li}_2(z) = -\log(1-z)/z$, the Taylor
coefficients $c_k$ satisfy a recurrence of order 2 with polynomial
coefficients in $k$, $a$ and $a - 1$, so the partial sums are products
of $3 \times 3$ matrices, computed with the forward recurrence or by
binary splitting. The coefficients are bounded by
$$
|c_k| \le \frac{1 + |\log(1-a)|}{\min(|a|, |1-a|)^k},
$$
which follows from the symbolic derivatives, or by $|c_k| \le 4/|1-a|^k$
for $|a| \le 1/2$ and $a = \pm i$ (from Cauchy's integral formula); the
tail after $n$ terms is then bounded by a geometric series in $|z -
a|/\min(|a|, |1-a|)$, which must be less than 1. If the segment from
$a$ to $z$ crosses the branch cut, the result is the analytic
continuation along the segment.

## The bit-burst algorithm {#bit-burst}

For a full-precision $z$, the Taylor series at $a$ needs $O(p)$ terms
with full-precision coefficients. `acb_hypgeom_dilog_bitburst`
([](algorithms/special-functions/hypergeometric#bit-burst)) instead
chooses $z_0$ as $z$ truncated to 16 fractional bits and the points
$z_j$ as $z$ truncated to $2^{j+4}$ bits, and sums
$$
\operatorname{Li}_2(z) - \operatorname{Li}_2(z_0) = \sum_{j}
\left(\operatorname{Li}_2(z_{j+1}) - \operatorname{Li}_2(z_j)\right),
$$
each term by Taylor continuation from $z_j$ with binary splitting: the
argument $z_{j+1} - z_j$ is small and short, so each step costs
$O(M(p)\log p)$ and the total $O(M(p)\log^2 p)$. The value
$\operatorname{Li}_2(z_0)$ comes from the series at 0 (or from a
continuation from a special point), and since $z_0$ has only 16 bits
it is cheap. If the path would come close to the branch cut, the
bit-burst step is skipped.

## Use in FLINT {#use}

- `acb_hypgeom_dilog` and the real `arb_hypgeom_dilog`; the methods
  `acb_hypgeom_dilog_zero`, `acb_hypgeom_dilog_zero_taylor`,
  `acb_hypgeom_dilog_transform`, `acb_hypgeom_dilog_bernoulli`,
  `acb_hypgeom_dilog_continuation`, `acb_hypgeom_dilog_bitburst`.
- `acb_polylog` with $s = 2$ (through `_acb_poly_polylog_cpx`) and the
  generic `gr_dilog` for complex ball contexts call `acb_hypgeom_dilog`.
