---
title: The Epstein zeta function
summary: Epstein zeta functions of lattices in any dimension by Crandall's representation, with rigorous truncation from the theta tail bound; the scaled incomplete gamma function; a double-precision ball kernel with certified tables of the incomplete gamma function; Mellin integrals of theta functions by Gauss–Legendre quadrature with Bernstein-ellipse bounds and direct sums for large arguments; and, for rational lattices, sums by shells with exact counts from products of theta series, enumeration or sums of three squares, used for the Madelung constant.
order: 52
---

The Epstein zeta function of the lattice $\Lambda = A \Z^d$ with shifts
$x$ and $y$ is
$$
Z_{\Lambda,\nu}\!\begin{vmatrix} x \\ y \end{vmatrix}
= \sum_{z \in \Lambda}{}' \frac{e^{-2 \pi i\, y \cdot z}}{|x - z|^{\nu}},
$$
the prime omitting $z = x$, continued meromorphically from
$\operatorname{Re}(\nu) > d$ (see [`acb_mzeta`](acb_mzeta) for the
conventions and the regularised function). All algorithms in FLINT
start from Crandall's representation [@BBG2024], which splits $Z$ into
two rapidly converging sums over $\Lambda$ and the dual lattice
$\Lambda^* = A^{-T} \Z^d$ whose terms are incomplete gamma functions.
They differ in how these sums are evaluated: term by term in ball
arithmetic, in double-precision ball arithmetic, as Mellin integrals of
theta functions, or, for rational lattices, by shells of equal norm.
Throughout, $e(t) = e^{2\pi i t}$ and
$$
G_a(t) = t^{-a} \Gamma(a, t) = \int_1^\infty w^{a-1} e^{-t w}\, dw.
$$

## Reductions {#reductions}

The basis is first LLL-reduced [@LLL1982] through its Gram matrix
(`arb_mat_spd_lll_reduce`), which leaves $Z$ unchanged and makes the
Cholesky factors of the quadratic forms well conditioned. Writing
$x = A(u_0 + m)$ and $A^T y = w = w_0 + l$ with $m, l \in \Z^d$ and
$u_0, w_0$ near zero,
$$
Z_{\Lambda,\nu}\!\begin{vmatrix} x \\ y \end{vmatrix}
= e(-w \cdot m)\, Z_{\Lambda,\nu}\!\begin{vmatrix} A u_0 \\ A^{-T} w_0 \end{vmatrix},
$$
by the substitution $z \to z + Am$ and periodicity in $y$. A shift
whose ball contains a lattice point without being exactly that point
gives an indeterminate result, since $Z$ is discontinuous (or singular)
there.

## Crandall's representation {#crandall}

::: theorem [@BBG2024, Thm. 2.11] {#thm-crandall}
For $\lambda > 0$ and $V = |\det A|$,
$$
\begin{aligned}
Z = \Bigl(\frac{\pi}{\lambda^2}\Bigr)^{\nu/2} \Biggl[
&\frac{1}{\Gamma(\nu/2)} \sum_{z \in \Lambda}{}' G_{\nu/2}\Bigl(\frac{\pi |z - x|^2}{\lambda^2}\Bigr) e(-y \cdot z) \\
&+ \frac{\lambda^d}{V\, \Gamma(\nu/2)} \sum_{k \in \Lambda^*}{}' G_{(d - \nu)/2}\bigl(\pi \lambda^2 |k + y|^2\bigr) e(-x \cdot (k + y)) \\
&- \frac{[x \in \Lambda]}{\Gamma(\nu/2 + 1)}
+ \frac{\lambda^d}{V} \frac{2\, [y \in \Lambda^*]}{(\nu - d)\, \Gamma(\nu/2)} \Biggr],
\end{aligned}
$$
where the primes omit $z = x$ and $k = -y$ and $[\cdot]$ is 1 if the
condition holds and 0 otherwise.
:::

::: proof
For $\operatorname{Re}(\nu) > d$, write
$|x - z|^{-\nu} = \frac{(\pi/\lambda^2)^{\nu/2}}{\Gamma(\nu/2)} \int_0^\infty s^{\nu/2 - 1} e^{-\pi s |x - z|^2 / \lambda^2} ds$
and split the integral at $s = 1$. The part over $[1, \infty)$ gives the
first sum. In the part over $[0, 1]$, Poisson summation over $\Lambda$
turns $\sum_z e(-y \cdot z) e^{-\pi s |x - z|^2/\lambda^2}$ into
$\frac{\lambda^d}{V s^{d/2}} \sum_k e^{-\pi \lambda^2 |k + y|^2 / s} e(-x \cdot (k + y))$,
and $s \to 1/s$ gives the second sum. The sum over $z$ in the Poisson
formula includes $z = x$ if $x \in \Lambda$, whose contribution
$\int_0^1 s^{\nu/2 - 1} ds = 2/\nu$ is subtracted, and the term
$k = -y$, if $y \in \Lambda^*$, is $\int_1^\infty s^{(d - \nu)/2 - 1} ds = 2/(\nu - d)$.
Since $G_a(t)$ is entire in $a$ for $t > 0$ and $1/\Gamma$ is entire,
the right-hand side continues $Z$ to all $\nu \ne d$.
:::

FLINT takes $\lambda = V^{1/d}$, which balances the two sums. With the
upper triangular Cholesky factors $C_D$ of $(\pi/\lambda^2) A^T A$ and
$C_R$ of $\pi \lambda^2 (A^T A)^{-1}$, the arguments are
$t = |C_D(n - u_0)|^2$ and $t = |C_R(m + w_0)|^2$ for $n, m \in \Z^d$.
For the regularised function, the term of $y$ itself is replaced by
$G_b(t) - \Gamma(b) t^{-b} = -t^{-b} \gamma(b, t)$, $b = (d - \nu)/2$,
which removes $\hat s_\nu(y)$; at $y = 0$ this is the singular term.

## Truncation {#truncation}

::: lemma {#lem-g-bound}
For $t > v = \max(0, \operatorname{Re}(a) - 1)$,
$|G_a(t)| \le e^{-t} / (t - v)$.
:::

::: proof
$|w^{a-1}| = w^{\operatorname{Re}(a) - 1} \le e^{v (w - 1)}$ for $w \ge 1$
(as $\log w \le w - 1$), so
$|G_a(t)| \le \int_1^\infty e^{v(w-1) - t w} dw = e^{-t}/(t - v)$.
:::

::: theorem Tail of the Crandall sums {#thm-crandall-tail}
For $T \ge \max(4, v + 1)$ and any shift $c$,
$$
\sum_{|C(n + c)|^2 > T} |G_a(|C(n + c)|^2)|
\le K\, T^{(d-1)/2} \frac{e^{-T}}{T - v},
\qquad K = \Bigl(1 + \sqrt{8/\pi}\Bigr) \prod_j \Bigl(1 + \frac{\sqrt{2\pi}}{C_{jj}}\Bigr).
$$
:::

::: proof
By [](#lem-g-bound) each term is at most $e^{-t}/(T - v)$, and the sum
of $e^{-t}$ over the points outside the ellipsoid is bounded by the
theta tail bound [](algorithms/special-functions/theta-functions#thm-tail)
with $R^2 = T$.
:::

$T$ is the smallest value making the bound at most $2^{-p}$ for the
target precision $p$, so $T \approx p \log 2 + O(d \log p)$, and the
lattice points of both ellipsoids $|C(n + c)|^2 \le T$ number
$O(T^{d/2} / \det C)$.

## The scaled incomplete gamma function {#gamma-upper}

The sums need $G_a(t)$ at many real $t > 0$ with one $a$, each to a
relative accuracy of its own (tapered by the size $e^{-t}$ of the term),
so `_acb_mzeta_gamma_upper_scaled` chooses the method and the working
precision itself:

- **Asymptotic series.** $F(t) = e^t G_a(t) = \int_0^\infty (1 + v)^{a-1} e^{-tv} dv$
  and Taylor's theorem with the integral remainder for $(1 + v)^{a-1}$
  give
  $F(t) = \sum_{k < N} (a-1) \cdots (a-k)/t^{k+1} + R_N$ with
  $|R_N| \le |(a-1) \cdots (a-N)| / t^{N+1}$ for $N \ge \operatorname{Re}(a) - 1$
  (then $|(1 + \theta v)^{a - 1 - N}| \le 1$). It is used when the
  smallest term reaches the tolerance, that is for $t$ of the order of
  the precision or larger.
- **Convergent series.** $G_a(t) = \Gamma(a) t^{-a} - e^{-t} \sum_k t^k/(a)_{k+1}$,
  whose cancellation (about $t/\log 2$ bits, more near the poles of
  $\Gamma$) is compensated by raising the working precision by an
  estimate of the loss. For real $a > -1$ with a short dyadic value (the
  usual case, $a = \nu/2$ or $(d - \nu)/2$ for a dyadic $\nu$) the sum
  $\frac1a \sum_k t^k/(a+1)_k$ is computed by rectangular splitting
  (`_arb_hypgeom_gamma_lower_sum_rs_1`), with no division per term; it is
  stopped at an $N$ with $a + 1 + N \ge 2t$, after which the ratios of
  consecutive terms are at most $\tfrac12$, so the tail is at most twice
  the first omitted term.
- **Nonpositive integers** $a = -n$:
  $G_a(t) = t^n \frac{(-1)^n}{n!} \bigl(E_1(t) - e^{-t} \sum_{k<n} (-1)^k k!/t^{k+1}\bigr)$
  with $E_1(t) = -\gamma - \log t + e^{-t} \sum_{k \ge 1} H_k t^k/k!$, a
  series of positive terms.

A ball $t$ is evaluated at its midpoint and the radius is propagated
with $|G_a'(s)| = |G_{a+1}(s)| \le G_{\operatorname{Re}(a)+1}(\underline{t})$
at the lower endpoint $\underline t$, which avoids the amplification of
the radius by about $e^t$ that the series would cause.

## Term-by-term summation {#arb}

`acb_mzeta_epstein_crandall` evaluates both sums of
[](#thm-crandall) in ball arithmetic. The points of the ellipsoids are
enumerated with `acb_theta_eld_set`. If the Gram matrix is exact, the
quadratic form is exact and the points are grouped by its value, so
that $G_a$ is evaluated once per shell; $\pm n$ are paired when the
center is a lattice point. Each value is computed at a precision
tapered by its magnitude $e^{-t}$. The cost is $O(p^{d/2})$ incomplete
gamma values.

## Double-precision balls {#dfloat}

At precisions up to 95 bits and for real $\nu$,
`acb_mzeta_epstein_crandall_dfloat` evaluates the same sums with
[`dfloat`](dfloat) balls of one or two doubles. The lattice points are
enumerated line by line (Fincke–Pohst [@FP1985], the coordinate $n_0$
varying along a line) in outward-rounded double interval arithmetic on
enclosures of $C$ and the shift, so that no point with $t \le T$ is
missed (extra points are harmless), and processed in batches with the
vector operations of `gr`: per line the partial quadratic forms and
phase offsets, per point $t = y_0^2 + \gamma$, $e^{-t}$
(`_gr_vec_exp`), the factor $F(t) = e^t G_a(t)$ from a table, the phase
reduced modulo 1 exactly, $\sin$ and $\cos$ (`_gr_vec_sin_cos`), and two
dot products. Dyadic phases come from an exact table of $e(r/2^q)$
indexed by an integer computation of $r$. Points with $t < 2^{-12}$
(very close to the shift) are evaluated with
`_acb_mzeta_gamma_upper_scaled`.

**Tables of $F$.** $F$ satisfies $t F' = (t - a) F - 1$, so its Taylor
coefficients $c_k$ at $t_0$ obey
$$
t_0 (k+1) c_{k+1} = (t_0 - a - k) c_k + c_{k-1} - [k = 0].
$$
The table holds Taylor polynomials of degree 28 (resp. 15 for one
double) on intervals $|t - t_0| \le \rho$ with centers decreasing
geometrically from $T + 2$ to $2^{-12}$, built once per $a$ and kept in
the cache.

::: lemma Coefficient bounds {#lem-ftab-coeffs}
With $s = \max(0, a - 1)$ and $\tau = t_0 - s > 0$,
$|c_k| \le \tau^{-(k+1)}$. For $a > 0$ and any $t_0 > 0$,
$|c_k| \le B\, (t_0/2)^{-k}$ with
$B = e^{t_0/2} (t_0/2)^{-a} \Gamma(a)$.
:::

::: proof
$c_k = \frac{(-1)^k}{k!} \int_0^\infty v^k (1+v)^{a-1} e^{-t_0 v} dv$ and
$(1 + v)^{a-1} \le e^{s v}$ give the first bound. For the second, $F$ is
analytic for $\operatorname{Re}(t) > 0$ with
$|F(t)| \le \int_0^\infty (1+v)^{a-1} e^{-\operatorname{Re}(t) v} dv = e^{r} r^{-a} \Gamma(a, r) \le e^r r^{-a} \Gamma(a)$,
$r = \operatorname{Re}(t)$, and Cauchy's inequality on the circle of
radius $t_0/2$ applies.
:::

The first bound is used with $\rho = \tau/16$ when $s \le t_0/4$, the
second with $\rho = t_0/32$ otherwise; either gives the truncation error
of the polynomial on the interval. The values $F(t_0)$ at the centers
are propagated downward from one evaluation at the top, by evaluating
the series of the previous center (of degree 60, resp. 32) at the step
$h = t_1 - t_0$.
Downward is the stable direction: the homogeneous solution
$\eta(t) = e^t t^{-a}$ decreases relative to $F$.

::: proposition Propagated error {#prop-ftab-prop}
If the midpoint of $F(t_0)$ is off by $\delta$, the coefficients
computed from it are those of $F$ plus $\delta$ times those of
$g(y) = \eta(t_0 + y)/\eta(t_0) = e^y (1 + y/t_0)^{-a}$, and the
value at the next center is off by at most
$$
|\delta| \Bigl( e^{h} (t_1/t_0)^{-a} + M \frac{q^{K+1}}{1 - q} \Bigr),
\qquad M = e^{t_0/2} 2^{|a|},\ q = 2|h|/t_0 \le \tfrac14,
$$
for a polynomial of degree $K$, plus the truncation error of $F$.
:::

::: proof
The recurrence is linear and $\eta$ solves its homogeneous version, so
the computed series is $F(t_0 + y) - \delta g(y)$ truncated at degree
$K$. The truncation of $g$ is $g(h)$ minus its tail; the coefficients of
$g$ are majorized by those of $e^x (1 - x/t_0)^{-|a|}$, which is bounded
by $M$ on $|x| = t_0/2$, so by Cauchy's inequality the tail is at most
$M q^{K+1}/(1-q)$. So the propagated error follows the homogeneous
solution, $\delta \eta(t)/\eta(t_0)$, and its size relative to $F$
decreases with $t$; bounding $g(h)$ by the majorant
$e^{|h|}(1 - |h|/t_0)^{-|a|}$ instead would multiply the error by more
than $e^{|h|}$ at every step.
:::

Within an interval, $|g(y)| \le e^\rho (1 - \rho/t_0)^{-|a|}$ bounds the
effect of $\delta$, and the polynomial is evaluated by Horner's rule in
plain double or double-double arithmetic with the a priori bound
$\sum_k (k+1) |c_k| \rho^k$ times the unit roundoff of a step
($2^{-52}$ with fused multiply-add, $2^{-102}$ for double-double, with
the terms of degree at least 14 summed in plain double and bounded
separately). Together these give an absolute bound $E$ for the value and
$D$ for the derivative on each interval, from which a ball input gets
its radius.

## Mellin integrals of theta functions {#theta}

With $Q$ the matrix of one of the two quadratic forms, $\sigma$ its
shift, $t_n = (n - \sigma)^T Q (n - \sigma)$ and phases $e(-p \cdot n)$,
each Crandall sum is
$$
S = \sum_n{}' e(-p \cdot n)\, G_a(t_n) = \int_1^\infty w^{a-1}\, \Theta(w)\, dw,
\qquad
\Theta(w) = \sum_n{}' e(-p \cdot n)\, e^{-w t_n},
$$
by the integral form of $G_a$ (the sums converge absolutely for $w \ge 1$),
and the full sum over $n$ is a Riemann theta function,
$$
\sum_n e(-p \cdot n) e^{-w t_n} = e^{-w \sigma^T Q \sigma}\, \theta\bigl(-p - i w Q \sigma/\pi,\; i w Q/\pi\bigr),
$$
evaluated with `acb_theta_one` in quasi-linear time; the excluded
points are subtracted. `acb_mzeta_epstein_theta` computes the integral
by Gauss–Legendre quadrature in $y = \log w$: the integrand
$f(y) = e^{a y} \Theta(e^y)$ is analytic for $|\operatorname{Im} y| < \pi/2$,
where $\operatorname{Re}(w) > 0$.

::: lemma Bounds for the theta series {#lem-phi}
For $r > 0$ and the Cholesky factor $C$ of $Q$,
$$
\Phi(r) = \sum_n{}' e^{-r t_n} \le \min\Bigl(\Phi_0(r),\ e^{-(r - r_1) t_{\min}} \Phi_0(r_1)\Bigr),
\qquad \Phi_0(r) = \prod_j \Bigl(1 + \frac{\sqrt{\pi/r}}{C_{jj}}\Bigr),
$$
for $0 < r_1 \le r$ and a lower bound $t_{\min}$ of the $t_n$ of the
points that are not excluded. If $\operatorname{Re}(w) \ge r$, then
$|\Theta(w)| \le \Phi(r)$.
:::

::: proof
$t_n = \sum_j (C_{jj} x_j + \sum_{k > j} C_{jk} x_k)^2$ with $x = n - \sigma$.
Summing over $n_0$ first, with the other coordinates fixed, the terms
$e^{-r (C_{00} n_0 + \beta)^2}$ form a unimodal sequence, whose sum is
at most its maximum (at most 1) plus the integral
$\sqrt{\pi/r}/C_{00}$; repeating for $n_1, n_2, \ldots$ gives
$\Phi_0(r)$. For the second bound, $e^{-r t} \le e^{-(r - r_1) t_{\min}} e^{-r_1 t}$
for $t \ge t_{\min}$. Finally $|e^{-w t}| = e^{-\operatorname{Re}(w) t}$.
:::

The interval $[0, \log W]$ is covered by intervals of half-length $h$.
On the Bernstein ellipse $E_\rho$ of an interval with center $y_c$ and
semi-minor axis below $\pi/2$, $\operatorname{Re}(w) \ge e^{y_c - h(\rho + 1/\rho)/2} \cos\theta_{\max}$
and $|e^{a y}| \le e^{\operatorname{Re}(a) \operatorname{Re}(y) + |\operatorname{Im}(a)| \theta_{\max}}$,
so [](#lem-phi) bounds $|f|$ by some $M$, and $n$ nodes have error at
most $\frac{64}{15} h M \rho^{-2n}/(\rho^2 - 1)$
([](algorithms/ball-arithmetic/integration#thm-gl)). For the tail,
$\Phi(w) \le e^{-(w - W/2) t_{\min}} \Phi_0(W/2)$ for $w \ge W/2$ and
$w^{s-1} \le W^{s-1} e^{(s-1)(w - W)/W}$, $s \ge \operatorname{Re}(a)$,
give
$$
\int_W^\infty |w^{a-1}|\, \Phi(w)\, dw \le \Phi_0(W/2)\, e^{W t_{\min}/2}\,
\frac{W^{s-1} e^{-W t_{\min}}}{t_{\min} - \max(0, s-1)/W}.
$$
The node counts taper with these bounds from interval to interval, and
the nodes are shared between the intervals and the two sums: a degree
already computed is reused when it is at most $1/8$ above the one
required. Nodes come from `arb_hypgeom_legendre_p_ui_root`, half of
them by symmetry.

**Direct sums for large $w$.** For large $w$ the series $\Theta(w)$ has
few terms of any size. The points with $t_n \le T_0$ are enumerated once
per sum, with their $t_n$ and phases; then
$$
\Theta(w) = \sum_{t_n \le T} e(-p \cdot n)\, e^{-w t_n} + \epsilon,
\qquad |\epsilon| \le e^{-(w-1) T} \Phi_0(1) \quad (w \ge 1),
$$
since $e^{-w t} \le e^{-(w-1) T} e^{-t}$ for $t > T$. A node uses the
direct sum when the $T$ needed for its accuracy is at most $T_0$, and
$T_0$ is chosen so that at most $\min(80, \max(16, 80000/p))\, 2^d$
points are stored, about where one real exponential per point stops
being cheaper than a theta value. This replaces the reduction and the
quasi-linear algorithm of `acb_theta_one` on most nodes of the upper
intervals.

**Working precision.** The value $w_k f(y_k)$ of a node is needed to an
absolute accuracy of about $2^{-p}/n$. With $\ell = \log_2 |w_k e^{a y_k}|$,
$\Theta(w_k)$ is needed to absolute accuracy $2^{-(p + \ell + \log_2 n)}$. A
direct sum has relative accuracy and needs about $p + \ell + \log_2 \Phi(w_k)$
bits, which decreases rapidly for large $w_k$; `acb_theta_one` has
absolute accuracy (relative to $e^{-w \sigma^T Q \sigma} \theta$), so a
theta value needs $p + \ell$ bits.

The cost is about $p$ integrand values per sum, independent of the
$O(p^{d/2})$ lattice points of the truncated Crandall sums, so the
method wins at high precision: measured crossovers against
term-by-term summation are about 700 bits in dimension 3, 400 bits in
dimension 4 and 220 bits in dimension 5; in dimension 2 the two break
even at about 4000 bits.

## Sums by shells {#shells}

For an integral positive definite Gram matrix $F$ (the lattice Gram
matrix being $F/L$) and rational $u$, $v$, the function
$Z = \sum_n{}' e(-v \cdot n)\, F(n - u)^{-s}$, $s = \nu/2$, is computed
by `acb_mzeta_epstein_gram_shell`. The same Mellin splitting as in
[](#thm-crandall), at a free parameter $\kappa$, gives
$$
Z = \frac{1}{\Gamma(s)} \left[\kappa^s S_D + \frac{\pi^{d/2} \kappa^{s - d/2}}{\sqrt{\det F}}\, e(-u \cdot v)\, S_R \right] + \text{singular terms},
$$
$$
S_D = \sum_n{}' e(-v \cdot n)\, G_s(\kappa F(n - u)), \qquad
S_R = \sum_m{}' e(-u \cdot m)\, G_{d/2 - s}\bigl(\pi^2 F^{-1}(m + v)/\kappa\bigr),
$$
with singular terms $-\kappa^s/\Gamma(s+1)$ if $u \in \Z^d$ and
$\pi^{d/2} (\det F)^{-1/2} \kappa^{s - d/2} / ((s - d/2) \Gamma(s))$ if
$v \in \Z^d$.

**Shells.** On each side the values of the integral quadratic form at
the shifted points, $V(n) = (S n - s)^T F (S n - s)$ with $s = S\sigma$
for the common denominator $S$ of the shift $\sigma$, are $g M + r_0$ with integers $M \ge 0$,
$g$ the gcd of all differences $V(n) - V(0)$ (computed from the values
at $\pm e_i$ and the cross terms $2 S^2 F_{ij}$). The points are counted
by shell $M$ and by the residue $N \ell \cdot n \bmod N$ of the linear
character, so that each side is
$\sum_M \sum_r C_r(M)\, e(r/N)\, G(\gamma (g M + r_0))$: one incomplete
gamma value per shell instead of one per point.

**Choice of $\kappa$.** $\kappa$ is chosen so that the arguments on one
side are $\lambda (g M + r_0)$ with a dyadic $\lambda$. There the values
$G_a(\lambda V)$ at exact rationals come from binary splitting (the
asymptotic series or the series at 0,
`_arb_hypgeom_gamma_upper_fmpq_inf_bsplit_sum` and
`_arb_hypgeom_gamma_lower_fmpq_0_bsplit_sum`, or the convergent series
for $\Gamma(-n, z)$). The arguments on the other side are irrational.
If its index is a positive integer $j$ (the Coulomb potential
$\nu = d - 2$, and even $\nu$), then
$G_j(t) = (j-1)!\, e^{-t} \sum_{i<j} t^{i-j}/i!$, and that side is a sum
of $j$ power series in $q = e^{-\gamma g}$ with rational coefficients
$C_r(M)/(gM + r_0)^{j-i}$. These are summed by rectangular splitting
with one shared fixed-point table of powers of $q$, read once per block
for the sums over all residues and powers; the denominators are
cleared in chunks whose integer multipliers fit in one or two limbs,
each chunk is divided once, and the divisions of several chunks are
interleaved [@GraMol2010]. Otherwise each value on that side is
computed below 16384 bits as in [](#gamma-upper), which is faster there,
and above by a bit-burst step, whose cost grows more slowly with the
precision: binary splitting at a dyadic $z_0$ of
about $2\sqrt{p}$ bits near $z$, then the Taylor series in
$\varepsilon = z - z_0$ with coefficients $g_k$ from the differential
equation $t G'' + (1 + a + t) G' + a G = 0$,
$$
g_{k+2} = -\frac{(k+1)(k+1+a+z_0)\, g_{k+1} + (k+a)\, g_k}{z_0 (k+1)(k+2)},
$$
and the bound $|g_k| \le \Gamma(k+a)/(k!\, z_0^{k+a})$ for $k + a > 0$
from the integral form of $G_a$, so that the tail after degree $K$ is
at most its first term over $1 - q$, $q = \max(1, \frac{K+a}{K+1}) |\varepsilon|/z_0$.
$\lambda$ balances the costs: about $T/(\lambda g_B)$ incomplete gamma
values on the binary splitting side against about
$T \lambda \Pi/(\pi^2 g_O)$ terms on the other.

Both sides are truncated with [](#thm-crandall-tail). The sums are split
into tasks for `flint_parallel_do` and accumulated with an
order-independent fixed-point accumulator: every contribution is
rounded to the grid $2^{-E}$ on its own and added as an exact integer,
with the radius in whole ulps, so the result does not depend on the
number of threads.

**Counts.** The integral form is first divided by the content of its
entries. The counts $C_r(M)$ come from interchangeable methods, each of
which recognises the forms it applies to and estimates its cost; the
cheapest applicable one is used (estimates independent of the number of
threads, so that the choice, and hence the result, is too), and the
sums read the counts block by block.

- **Products of theta series.** $\Z^d$ is split into mutually
  $F$-orthogonal sublattices $L_b$ of rank 2 (and one of rank 1 for odd
  $d$), with bases from an LLL-reduced basis and the integral kernels of
  the previous blocks, so that $\sum L_b$ has finite index. On a coset
  $v + \sum L_b$, with $y = S v - s$, $c_b = B_b F y$ and
  $G_b = B_b F B_b^T$, the points $n = v + \sum_b B_b^T u_b$ have
  $$
  V(n) = \sum_b (S u_b + t_b)^T G_b (S u_b + t_b), \qquad t_b = G_b^{-1} c_b,
  $$
  and $\mathrm{Den} \cdot (S u + t)^T G (S u + t)$ is an integer for the
  common denominator $\mathrm{Den}$ of the $c_b^T G_b^{-1} c_b$. So the
  generating function $\sum_n x^{\mathrm{Den}\, V(n)} \zeta^{N \ell \cdot n}$
  of a coset is the product of the series of the blocks, each
  enumerated directly in linear time (a rank-2 block has about as many
  points as exponents). The $\lceil d/2 \rceil - 1$ products run in
  $\Z[x][\zeta]/(\zeta^N - 1)$ by Kronecker substitution with
  `nmod_poly` arithmetic modulo a prime just above a bound for every
  count. The exponents of a block on a coset lie in one class modulo a
  common step $q$, so the series are taken in $x^q$; among the choices of
  the first block the one with the least total series length is used.
- **Enumeration.** Fincke–Pohst enumeration with exact integer values,
  for any form.
- **Sums of three squares.** For $F = I_3$ and characters of order at
  most 2, the counts are streamed from $r_3$ without a table of all
  shells: $r_3(M)$ for integral shifts, and for the shift
  $(\tfrac12, \tfrac12, \tfrac12)$ with trivial character
  $r_3(8k + 3) = 8 [x^k] T(x)^3$, $T(x) = \sum_{a \ge 0} x^{a(a+1)/2}$.
  The coefficients $t_2(j) = [x^j] T(x)^2 = \sum_{d \mid 4j+1} \chi_4(d)$
  (from $8j + 2 = (2a+1)^2 + (2b+1)^2$ and
  $r_2(n) = 4 \sum_{d \mid n} \chi_4(d)$) come from a segmented divisor
  sieve over $d \le \sqrt{4j + 1}$, since $d$ and $(4j+1)/d$ have the same
  character, and fit in 2 bytes; $r_3(8k+3)$ then needs one sparse
  convolution with $T$ per block, or one FFT for short sums.

## The Madelung constant {#madelung}

The Madelung constant of sodium chloride is
$M = -Z_{\Z^3, 1}\bigl(0; (\tfrac12, \tfrac12, \tfrac12)\bigr)$, computed
by the shell algorithm. With $\kappa = \lambda$ the splitting is
Crandall's expansion [@Cra1999, eq. (1.2)]
$$
-M = \sum_{N \ge 1} (-1)^N r_3(N) \frac{\operatorname{erfc}(\sqrt{\lambda N})}{\sqrt{N}}
+ \frac{4}{\pi} \sum_{M \equiv 3 \bmod 8} r_3(M) \frac{e^{-\pi^2 M/(4\lambda)}}{M}
- 2 \sqrt{\lambda/\pi}.
$$
The first sum has $O(p/\lambda)$ terms, each an incomplete gamma value at
an exact integer by binary splitting; the second is a power series in
$q = e^{-2\pi^2/\lambda}$ with $O(\lambda p)$ terms, whose coefficients
come from the $t_2$ sieve. Both cost about $p^2$ up to logarithmic
factors, and $\lambda$ balances them.

## Use in FLINT {#use}

- `acb_mzeta_epstein` and `acb_mzeta_epstein_generic` choose between
  the algorithms below.
- `acb_mzeta_epstein_crandall`: term-by-term summation.
- `acb_mzeta_epstein_crandall_dfloat`: double-precision balls.
- `acb_mzeta_epstein_theta`: Mellin integrals of theta functions.
- `acb_mzeta_epstein_shell`, `acb_mzeta_epstein_gram_shell`,
  `acb_mzeta_madelung_nacl`: sums by shells.
- `_acb_mzeta_gamma_upper_scaled`: the scaled incomplete gamma function.
