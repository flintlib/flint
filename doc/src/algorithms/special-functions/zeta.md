---
title: Zeta functions, polylogarithms and L-functions
summary: Euler–Maclaurin summation for the Hurwitz zeta function and its derivatives with a rigorous remainder bound, power sums, Taylor expansion in the parameter, the Riemann zeta function at integers, polylogarithms, and the methods used for Dirichlet L-functions.
complexity: $O(p + |t|)$ terms for $\zeta(s, a)$ by Euler–Maclaurin summation at precision $p$ and height $t$; $O(\sqrt{t})$ by the Riemann–Siegel formula; $O(p/\log p)$ prime powers for the Euler product at large $\operatorname{Re}(s)$.
order: 50
---

The Hurwitz zeta function

$$
\zeta(s, a) = \sum_{k=0}^\infty \frac{1}{(k+a)^s}, \qquad \operatorname{Re}(s) > 1,
$$

and its analytic continuation to $s \ne 1$ are the core of a family of
functions: the Riemann zeta function $\zeta(s) = \zeta(s, 1)$, the
Dirichlet L-functions, which are linear combinations of Hurwitz zeta
values, and the polylogarithms, which are linear combinations of two
Hurwitz zeta values with complex parameter. FLINT computes $\zeta(s, a)$
and its derivatives with respect to $s$ by Euler–Maclaurin summation
with the rigorous error bounds of [@Joh2013], and uses specialised
methods where they are faster: exact formulas and the Euler product at
integers, the Riemann–Siegel formula high on the critical line
([](algorithms/special-functions/zeta-zeros)), and a Taylor expansion in
$a$ when many values with the same $s$ are needed.

## Overview of methods {#methods}

| Function | Method |
|---|---|
| `acb_dirichlet_zeta`, `acb_zeta` | integers: Bernoulli numbers (negative), `arb_zeta_ui` (positive); $|t| \ge 24 p^{3/2}$: Riemann–Siegel; $\operatorname{Re}(s) < 0$: reflection formula; exact dyadic $s$ at high precision: approximate functional equation; otherwise Euler–Maclaurin |
| `acb_dirichlet_zeta_jet` | Euler–Maclaurin for any length (reflection for $\operatorname{Re}(s) \le -4$); Riemann–Siegel for length $\le 2$ high on the critical line; exact rational $s$ of small height at very high precision: approximate functional equation |
| `acb_dirichlet_hurwitz`, `acb_hurwitz_zeta` | special values ($a = 1$, $s = 0$, negative integers via Bernoulli polynomials); otherwise Euler–Maclaurin (`_acb_poly_zeta_cpx_series`) |
| `arb_zeta_ui` | $n > 0.7p$: asymptotic; even $n$: Bernoulli numbers or Euler product; $\zeta(3)$, $\zeta(5)$: constants; small odd $n$ at huge precision: Borwein's algorithm by binary splitting; large odd $n$: Euler product; otherwise Borwein's algorithm |
| `acb_dirichlet_hurwitz_precomp_eval` | Taylor expansion in $a$ on a grid ([](#hurwitz-taylor)) |
| `acb_polylog` | $|z| < 1/2$: direct series; $s = 2$: `acb_hypgeom_dilog`; otherwise two Hurwitz zeta values ([](#polylog)) |
| `acb_dirichlet_l` | [](#l-functions) |
| `acb_dirichlet_stieltjes` | integral representation [@JB2018], or a Hurwitz zeta series at $s = 1$ (via `acb_dirichlet_zeta_jet` for $a = 1$, so with the approximate functional equation at very high precision) |

## Euler–Maclaurin summation {#euler-maclaurin}

Applying the Euler–Maclaurin formula to $f(t) = (a + t)^{-s}$ on
$[N, \infty)$ gives the following representation, valid for all $s \ne 1$
by analytic continuation.

::: theorem Euler–Maclaurin formula for the Hurwitz zeta function {#thm-em}
Let $a > 0$ be real, $N \ge 0$ and $M \ge 1$ integers, and $s = \sigma + it$
with $\sigma + 2M > 1$, $s \ne 1$. Then

$$
\zeta(s, a) = \sum_{k=0}^{N-1} \frac{1}{(a+k)^s} + \frac{(a+N)^{1-s}}{s-1} + \frac{1}{2 (a+N)^s}
+ \sum_{j=1}^{M} \frac{B_{2j}}{(2j)!} \frac{(s)_{2j-1}}{(a+N)^{s+2j-1}} + R_{N,M}(s, a),
$$

where $(s)_k = s(s+1)\cdots(s+k-1)$ and

$$
|R_{N,M}(s, a)| \le \frac{4 |(s)_{2M}|}{(2\pi)^{2M}} \frac{(a+N)^{1-\sigma-2M}}{\sigma + 2M - 1}.
$$
:::

::: proof
For $\sigma > 1$, the Euler–Maclaurin formula for $f(t) = (a+t)^{-s}$,
whose derivatives are $f^{(m)}(t) = (-1)^m (s)_m (a+t)^{-s-m}$, gives

$$
\sum_{k=N}^{\infty} f(k) = \int_N^\infty f(t)\,dt + \frac{f(N)}{2} - \sum_{j=1}^{M} \frac{B_{2j}}{(2j)!} f^{(2j-1)}(N)
- \int_N^\infty \frac{\tilde B_{2M}(t)}{(2M)!} f^{(2M)}(t)\,dt,
$$

where $\tilde B_{2M}(t) = B_{2M}(t - \lfloor t \rfloor)$ is the periodic
Bernoulli function. The integral of $f$ is $(a+N)^{1-s}/(s-1)$, and the
terms of the sum are those of the theorem. Both sides are analytic in
$s$ for $\sigma + 2M > 1$, $s \ne 1$, since the last integral converges
there, so the formula holds in that domain. Finally
$|\tilde B_{2M}(t)| \le |B_{2M}| = 2 (2M)! \zeta(2M)/(2\pi)^{2M} < 4 (2M)!/(2\pi)^{2M}$,
so the remainder is at most
$4 |(s)_{2M}| (2\pi)^{-2M} \int_N^\infty (a+t)^{-\sigma-2M}\,dt$.
:::

[@Joh2013] extends the bound to complex $a$ (with a factor depending on
$\arg(a + N)$) and to the derivatives with respect to $s$, by treating
$s + x$ as a power series: the remainder of the expansion of
$\zeta(s + x, a)$ to length $d$ is bounded coefficientwise by an explicit
integral of $(a+t)^{-\sigma-2M} \log^k(a+t)$, which is what
`_acb_poly_zeta_em_bound` computes.

**Choice of parameters.** The terms of the Bernoulli sum first decrease
like $|s + 2j|^{2j}/(2\pi(a+N))^{2j}$ and then grow; the bound is small
when $N$ and $M$ are both large compared with $|s|$ and the precision.
`_acb_poly_zeta_em_choose_param` takes
$M = \min(N, p + \min(N/100, 2000))$ for a target of $p$ bits and finds
the smallest $N$ for which the bound is below $2^{-p}$ times an estimate
of $|a^{-s}|$, by doubling and bisection. For $s$ of moderate size this
gives $N, M = O(p)$. The cost is dominated by the $N$ powers
$(a+k)^{-s}$ and the $M$ Bernoulli numbers (cached by
`bernoulli_cache_compute`). The Bernoulli sum is evaluated term by term
(`_acb_poly_zeta_em_tail_naive`), or by binary splitting
(`_acb_poly_zeta_em_tail_bsplit`) for long power series (length at
least 5 and at least $M/10$). For large
$|t|$, the remainder bound requires $N \gtrsim |t|/(2\pi)$, so the cost
grows linearly with the height; this is why the Riemann–Siegel formula,
with $O(\sqrt{t})$ terms, takes over there.

**Power sums.** The main sum $\sum_{k<N} (a+k)^{-(s+x)}$ is computed by
`_acb_poly_zeta_em_sum` with one of several methods:

- for $a = 1$ and length $\le 2$: `acb_dirichlet_powsum_sieved`, which
  computes $k^{-s}$ only for prime $k$ and obtains the others as
  products $(ij)^{-s} = i^{-s} j^{-s}$ from a table, groups the even
  terms as a polynomial in $2^{-s}$ evaluated by Horner's rule, and so
  needs about $N/\log N$ exponentials and $N/2$ multiplications;
- for $a = 1$ and length $\le 4$: `acb_dirichlet_powsum_smooth`, which
  sieves only by the primes 2, 3 and 5 and uses logarithmic memory;
- for longer series or general $a$: a product tree
  (`_acb_poly_powsum_series_tree`), or a direct sum, in parallel when
  several threads are available.

The powers $k^{-s} = \exp(-s \log k)$ are computed with the logarithm of
a nearby integer when possible (`acb_dirichlet_powsum_term`), and the
precision of each term is lowered according to its magnitude.

**Reflection.** For $\operatorname{Re}(s) < 0$ the Riemann zeta function
uses the functional equation
$\zeta(s) = 2 (2\pi)^{s-1} \sin(\pi s/2) \Gamma(1-s) \zeta(1-s)$, also as
power series, since the cost of Euler–Maclaurin summation grows with
$|s|$.

## The zeta function at integers {#integers}

At even positive integers, $\zeta(2n) = (-1)^{n+1} B_{2n} (2\pi)^{2n}/(2 (2n)!)$
is computed from the Bernoulli number when it is cheap (it is cached), and
at negative integers $\zeta(-n) = -B_{n+1}/(n+1)$. For large $n$ the
Euler product

$$
\zeta(n) = \prod_p (1 - p^{-n})^{-1},
$$

truncated at the primes $p < P$, converges quickly: the relative error is
at most about $\zeta(n, P) \approx P^{1-n}/(n-1)$, so $P \approx 2^{p/n}$
primes suffice (`arb_zeta_ui_euler_product`, which bounds the tail by
`mag_hurwitz_zeta_uiui`). For odd $n$ of moderate size, Borwein's
algorithm [@Bor2000]

$$
\zeta(s) = \frac{-1}{d_n (1 - 2^{1-s})} \sum_{k=0}^{n-1} \frac{(-1)^k (d_k - d_n)}{(k+1)^s} + \varepsilon_n(s),
\qquad d_k = n \sum_{i=0}^{k} \frac{(n+i-1)! \, 4^i}{(n-i)! (2i)!},
$$

has an error $\varepsilon_n(s)$ of order $(3 + \sqrt 8)^{-n}$ for real
$s \ge 1/2$, so it needs $n \approx 0.39 p$ terms; `arb_zeta_ui_vec_borwein`
evaluates it for several $s$ at once, and `arb_zeta_ui_borwein_bsplit`
evaluates it by binary splitting at very high precision.

## Taylor expansion in the parameter {#hurwitz-taylor}

To evaluate $\zeta(s, a)$ at many $a$ for a fixed $s$, as needed for the
L-functions of all characters of a modulus,
`acb_dirichlet_hurwitz_precomp_init` precomputes Taylor expansions at
grid points. With $\zeta(s, a) = \sum_{j=0}^{A-1} (a+j)^{-s} + \zeta(s, a + A)$
and the expansion

$$
\zeta(s,a+x) = \sum_{k=0}^{\infty} (-x)^k \frac{(s)_k}{k!} \zeta(s+k,a),
$$

it stores $K$ coefficients at each of the $N$ points
$A + (2i+1)/(2N)$, $0 \le i < N$, which leaves $|x| \le 1/(2N)$ for any
$a \in (0, 1]$. Each evaluation then costs a polynomial evaluation of
length $K$ and $A$ powers. The truncation error is bounded as follows.

::: lemma {#lem-hurwitz-taylor}
Let $a \ge 1$ be real, $\sigma = \operatorname{Re}(s)$, $K + \sigma > 1$,
and let $T(k) = |x|^k \frac{|(s)_k|}{k!} \left[ a^{-\sigma-k} + \frac{a^{1-\sigma-k}}{\sigma+k-1} \right]$.
Then the tail of the series satisfies

$$
\sum_{k=K}^{\infty} |x|^k \frac{|(s)_k|}{k!} \zeta(\sigma+k,a) \le \sum_{k = K}^\infty T(k),
$$

and if

$$
C = \frac{|x|}{a} \left(1 + \frac{1}{K+\sigma+a-1}\right) \left(1 + \frac{|s-1|}{K+1}\right) < 1,
$$

then $\sum_{k \ge K} T(k) \le T(K)/(1 - C)$.
:::

::: proof
For real $\sigma + k > 1$,
$\zeta(\sigma + k, a) \le a^{-\sigma-k} + \int_a^\infty u^{-\sigma-k} du = a^{-\sigma-k} + a^{1-\sigma-k}/(\sigma + k - 1)$,
since the terms decrease, and $|\zeta(s + k, a)| \le \zeta(\sigma + k, a)$.
The ratio of consecutive terms is

$$
\frac{T(k+1)}{T(k)} = \frac{|x|}{a} \frac{k+\sigma-1}{k+\sigma} \frac{k+\sigma+a}{k+\sigma+a-1} \frac{|k+s|}{k+1}
\le \frac{|x|}{a} \left(1 + \frac{1}{K+\sigma+a-1}\right) \left(1 + \frac{|s-1|}{K+1}\right) = C
$$

for $k \ge K$, using $|k + s| \le k + 1 + |s - 1|$. The tail is then
dominated by a geometric series.
:::

With $a \ge A$ and $|x| \le 1/(2N)$, $C \approx 1/(2AN)$, so the error
is $O((2AN)^{-K})$; `acb_dirichlet_hurwitz_precomp_bound` evaluates this
bound and `acb_dirichlet_hurwitz_precomp_choose_param` minimises the
total cost $NK$ (precomputation) plus (number of evaluations) $\times$
$(K + A)$ subject to it.

## Polylogarithms {#polylog}

The polylogarithm $\operatorname{Li}_s(z) = \sum_{k=1}^\infty z^k k^{-s}$
for $|z| < 1$, continued analytically with a branch cut on
$z \in (1, \infty)$, is computed by `acb_polylog`, and its expansion in
$s$ by `acb_poly_polylog_series`.

**Small $|z|$.** For $|z| < 1/2$, the series is summed directly, also
for the derivatives with respect to $s$:
$\operatorname{Li}_{s+x}(z) = \sum_{d} x^d \frac{(-1)^d}{d!} \sum_{k \ge 1} T_d(k)$
with $T_d(k) = z^k \log^d(k) k^{-s}$. The tail after $N$ terms is bounded
by `mag_polylog_tail`: if $U(k)$ is nonincreasing and
$T_d(k+1)/T_d(k) \le U(k)$, then
$\sum_{k \ge N} T_d(k) \le T_d(N) / (1 - U(N))$ when $U(N) < 1$. One can
take $U(k) = |z| \, B(k, \max(0, -\sigma)) \, B(k\log k, d)$ with
$B(m, n) = (1 + 1/m)^n$, because $(k/(k+1))^{\sigma} \le (1 + 1/k)^{\max(0,-\sigma)}$
and $(\log(k+1)/\log k)^d \le (1 + 1/(k \log k))^d$.

**General $z$.** Otherwise, with $s = 1 - v$,

$$
\operatorname{Li}_s(z) = \frac{\Gamma(v)}{(2\pi)^v}
\left[ i^v \zeta\left(v, \frac{1}{2} + \frac{\log(-z)}{2\pi i}\right)
+ i^{-v} \zeta\left(v, \frac{1}{2} - \frac{\log(-z)}{2\pi i}\right) \right],
$$

with the principal branch of $\log(-z)$, which gives the conventional
branch cut. The Hurwitz zeta functions are evaluated by Euler–Maclaurin
summation with complex parameter. To expand in $s$, $v$ is replaced by
$v + x$, the right side becomes $\Gamma(v+x) [E_1 Z_1 + E_2 Z_2]$ with
$E_1 = (i/(2\pi))^{v+x}$ and $E_2 = (1/(2\pi i))^{v+x}$, and $x$ is
mapped to $-x$ at the end. At $v = 1$ the Hurwitz zeta functions have
poles; with the deflated $\tilde Z_j = Z_j - 1/x$,
$E_1 Z_1 + E_2 Z_2 = (E_1 + E_2)/x + E_1 \tilde Z_1 + E_2 \tilde Z_2$,
where $(E_1+E_2)/x$ is a power series since $E_1 + E_2$ vanishes at
$x = 0$. At $v = 0, -1, -2, \ldots$, where $\Gamma$ has a pole, the
product $x\Gamma(v+x) \cdot [E_1 Z_1 + E_2 Z_2]/x$ is formed instead.
These cancellations require $s$ to be an exact integer; a small ball
around $0, 1, 2, \ldots$ gives an indeterminate result. The method is
inefficient for $|z| \to 0$ (where the direct series is used) and for
$|z| \to \infty$, where an asymptotic expansion would be better. For
$s = 2$, `acb_hypgeom_dilog` is used instead
([](algorithms/special-functions/dilogarithm)).

## Dirichlet L-functions {#l-functions}

For a character $\chi$ modulo $q$,

$$
L(s, \chi) = \sum_{k=1}^{\infty} \frac{\chi(k)}{k^s} = q^{-s} \sum_{k=1}^{q} \chi(k) \, \zeta\!\left(s, \frac{k}{q}\right),
$$

and `acb_dirichlet_l` uses this decomposition (`acb_dirichlet_l_hurwitz`)
in general, computing the $\phi(q)$ Hurwitz zeta values with the Taylor
expansion of [](#hurwitz-taylor) when that is cheaper. At $s = 1$, the
deflated values $\zeta(s, a) - 1/(s-1)$ are used, since the poles cancel
for non-principal $\chi$. `acb_dirichlet_l_vec_hurwitz` computes
$L(s, \chi)$ for all $\chi$ modulo $q$ at once: the map from the vector
$(\zeta(s, k/q))_k$ to $(L(s, \chi))_\chi$ is a discrete Fourier
transform on the group of characters (`acb_dirichlet_dft_index`), so the
cost is $\phi(q)$ Hurwitz zeta values plus $O(\phi(q) \log \phi(q))$
operations.

For large $\operatorname{Re}(s)$, from about $8 + p/(2\log p)$, the Euler
product $1/L(s,\chi) = \prod_p (1 - \chi(p) p^{-s})$ is used; the
truncation at $p < P$ leaves an error at most
$\zeta(\operatorname{Re}(s), P)$ in $1/L$, because the truncated product
contains all terms $\mu(k)\chi(k)k^{-s}$ with $P$-smooth $k$.

For primitive $\chi$ of parity $\delta$, the completed function
$\Lambda(s, \chi) = (q/\pi)^{(s+\delta)/2} \Gamma((s+\delta)/2) L(s,\chi)$
satisfies $\Lambda(s, \chi) = \epsilon(\chi) \Lambda(1 - s, \overline\chi)$,
which moves evaluations from $\operatorname{Re}(s) < -1/2$ to the right
half-plane. The root number is
$\epsilon(\chi) = G(\chi)/(i^\delta \sqrt{q})$, where the Gauss sum
$G(\chi) = \sum_x \chi(x) e^{2\pi i x/q}$ factors over the prime powers
dividing $q$ by the Chinese remainder theorem. It can also be obtained
from the theta series
$\Theta(\chi, t) = \sum_{n \in \Z} \chi(n) n^\delta e^{-\pi n^2 t^2/q}$,
which satisfies $t^{1/2+\delta}\Theta(\chi, t) = \epsilon(\chi)\Theta(\overline{\chi}, 1/t)$,
so that $\epsilon(\chi) = \Theta(\chi, 1)/\overline{\Theta(\chi, 1)}$
whenever $\Theta(\chi, 1) \ne 0$; this needs $O(\sqrt{pq})$ terms.

The Mellin transform of the theta series gives the *approximate
functional equation*, a sum of incomplete gamma functions,

$$
\Lambda(s, \chi) = \sum_{n=1}^{\infty} \chi(n)
\left(\frac{q}{\pi}\right)^{\frac{s+\delta}{2}} \frac{\Gamma\bigl(\frac{s+\delta}{2}, \frac{\pi n^2}{q}\bigr)}{n^s}
+ \epsilon(\chi) \sum_{n=1}^{\infty} \overline{\chi}(n)
\left(\frac{q}{\pi}\right)^{\frac{1-s+\delta}{2}} \frac{\Gamma\bigl(\frac{1-s+\delta}{2}, \frac{\pi n^2}{q}\bigr)}{n^{1-s}},
$$

whose terms decrease like $e^{-\pi n^2/q}$, so $O(\sqrt{pq})$ terms
suffice at precision $p$ [@Coh2007b]. `acb_dirichlet_l_fmpq_afe` uses it
for rational $s$, for which the upper incomplete gamma function with
rational parameter can be computed efficiently (by its asymptotic series
when possible), with each term computed at the precision its magnitude
requires; it is faster than Euler–Maclaurin summation at very high
precision or for large $q$. `acb_dirichlet_l_fmpq_afe_jet` computes the
Taylor expansion in $s$ in the same way, with the incomplete gamma functions
expanded with respect to their parameter; this is used for jets of
$\zeta(s)$ and $L(s,\chi)$ at exact rational $s$ (`acb_dirichlet_zeta_jet`,
`acb_dirichlet_l_jet`), and thereby for Glaisher's constant and the
Stieltjes constants $\gamma_n$ at very high precision. The cost is
$\tilde O(p^{1.5})$ and the memory $O(p)$, compared with $\tilde O(p^2)$
and $O(p^2)$ for Euler–Maclaurin summation. Since
$\Gamma(s, z)$ is computed by binary splitting over rational
parameters, only $s$ of small height can be handled; from an `acb_t`,
these are the dyadic numbers.

## Other functions {#other}

**Lerch transcendent.** $\Phi(z, s, a) = \sum_{k \ge 0} z^k (k+a)^{-s}$
reduces to the Hurwitz zeta function for $z = 1$, to the Dirichlet eta
function, digamma or Hurwitz zeta functions for $z = -1$, and to
polylogarithms for $a = 1$. Otherwise `acb_dirichlet_lerch_phi` sums the
series directly when $|z| \le 0.875$, and else evaluates the Hankel
contour integral
$\Phi(z,s,a) = -\frac{\Gamma(1-s)}{2 \pi i} \int_C \frac{(-t)^{s-1} e^{-a t}}{1 - z e^{-t}} dt$
with `acb_calc_integrate`, deforming the contour around the poles
$t = \log z + 2\pi i k$.

**Stieltjes constants.** The generalised Stieltjes constants
$\gamma_n(a)$ are the coefficients of the Laurent expansion of
$\zeta(s, a)$ at $s = 1$. For moderate $n$, they are read off from the
Euler–Maclaurin expansion of $\zeta(s + x, a)$ at $s = 1$ with deflation
(for $a = 1$ and very high precision, from the approximate functional
equation instead); for large $n$, `acb_dirichlet_stieltjes` integrates a representation
along a contour that avoids the cancellation of the power series method
[@JB2018], which makes very large $n$, such as $10^{100}$, feasible.

## Use in FLINT {#use}

- `acb_dirichlet_zeta`, `acb_dirichlet_zeta_jet`, `acb_dirichlet_hurwitz`,
  `acb_dirichlet_eta`, `acb_dirichlet_xi` ([`acb_dirichlet`](acb_dirichlet));
  `acb_zeta` and `acb_hurwitz_zeta` are wrappers.
- `_acb_poly_zeta_cpx_series`, `_acb_poly_zeta_em_choose_param`,
  `_acb_poly_zeta_em_bound`, `_acb_poly_zeta_em_sum`,
  `_acb_poly_zeta_em_tail_bsplit`, `acb_poly_zeta_series`
  ([`acb_poly`](acb_poly)): Euler–Maclaurin summation.
- `acb_dirichlet_powsum_sieved`, `acb_dirichlet_powsum_smooth`,
  `_acb_poly_powsum_series_tree`: power sums.
- `arb_zeta_ui`, `arb_zeta_ui_euler_product`, `arb_zeta_ui_vec_borwein`,
  `arb_zeta_ui_borwein_bsplit`, `arb_zeta_ui_bernoulli` ([`arb`](arb)).
- `acb_dirichlet_hurwitz_precomp_init`, `acb_dirichlet_hurwitz_precomp_bound`,
  `acb_dirichlet_hurwitz_precomp_eval`.
- `acb_polylog`, `acb_poly_polylog_series`, `mag_polylog_tail`.
- `acb_dirichlet_l`, `acb_dirichlet_l_hurwitz`,
  `acb_dirichlet_l_euler_product`, `acb_dirichlet_l_vec_hurwitz`,
  `acb_dirichlet_l_fmpq_afe`, `acb_dirichlet_l_fmpq_afe_jet`, `acb_dirichlet_root_number`,
  `acb_dirichlet_gauss_sum`, `acb_dirichlet_theta_arb`.
- `acb_dirichlet_lerch_phi`, `acb_dirichlet_stieltjes`.

The tutorial [](tutorials/numerics/zeta) evaluates the zeta function,
its zeros and L-functions in complete programs.
