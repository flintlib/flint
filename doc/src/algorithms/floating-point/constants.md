---
title: Constants in mp_real
summary: The series and formulas behind the mp_real constants (π, log 2, e, Euler's constant, log 10, Catalan's constant, ζ(3), ζ(5), Γ(1/3), Γ(1/4)), the dual-number binary splitting for Euler's constant, and the cache of exact floors.
complexity: $O(M(n)\log^2 n)$ for each constant; $\Gamma(1/4)$ from $2^{16}$ limbs in $O(M(n)\log n)$ plus the cost of $\pi$.
order: 70
---

`mp_real` computes the constants used throughout FLINT's floating-point
code by binary splitting of hypergeometric series
([](algorithms/floating-point/hypergeometric-series)). The formulas were
taken from the collection of y-cruncher [@Yee2025], and for each
constant the formula measured fastest in this implementation was
chosen; the mathematical background of most of them is described in
[](algorithms/special-functions/constants).

## The series {#series}

| Constant | Formula |
|---|---|
| $\pi$ | Chudnovsky [@CC1988], $\pi = S/\sqrt{10005}$ (see below), 47.1 bits per term |
| $\log 2$ | Zuniga [@Zun2025], $\frac{1}{2160}\sum_{k \ge 0}(1497 + 1794k)\prod_{j=1}^k \frac{j(2j-1)}{216(6j+1)(6j+5)}$, 11.9 bits per term |
| $e$ | $\sum_{k \ge 0} 1/k!$ |
| $G$ | Pilehrood [@PP2010], $\frac{1}{64}\sum_{k \ge 1}\frac{256^k(580k^2 - 184k + 15)(2k)!^3(3k)!^2}{k^3(2k-1)(6k)!^2}$ |
| $\zeta(3)$ | a series of Zuniga (2023) with a numerator of degree 11 |
| $\zeta(5)$ | Sun's identity (2025), $\zeta(5) = (3S + 56\pi^2\zeta(3))/540$ with a ${}_3F_2$-type series $S$ |
| $\Gamma(1/3)$ | $(810^{1/4}\pi/X)^{1/3}$ with Guillera's series $X$ (2023) |
| $\Gamma(1/4)$ | $(\pi^6/(322S^4))^{1/8}$ with Ebisu's series $S$ (2016), or the AGM |

The Chudnovsky series in the format of `mp_real_hypgeom_series` is
$S = \bigl((c_Q + c_P\sum_{k \ge 1}\frac{P(k)}{Q(k)}\prod_{j<k}\frac{R(j)}{Q(j)})/c_D\bigr)^{-1}$
with $c_P = 1$, $c_Q = 13591409$, $c_D = 4270934400$,
$R(k) = (6k - 5)(2k - 1)(6k - 1)$, $P(k) = (13591409 + 545140134k)R(k)$
and $Q(k) = -(640320^3/24)k^3$; the example on [](mp_real/hypgeom) runs
it. Since `mp_real_hypgeom_series` gives $\bits(n - 1)$ correct bits,
the constants call it with one limb more than they need.

The parts of $\zeta(5)$ (the series, $\zeta(3)$ and $\pi$) and of the
$\Gamma$ values (the series and $\pi$) are independent and run on
separate threads. For $\Gamma(1/4)$ above $2^{16}$ limbs, the formula
$$
\Gamma(1/4) = \sqrt{\frac{(2\pi)^{3/2}}{\operatorname{agm}(1, \sqrt 2)}}
$$
is used instead, with the AGM ($O(M(n)\log n)$) computed in parallel
with $\pi$. $\log 10$ is $\log 2 + \log 5$ from the three-term
Machin-type formula of the first three primes
($\operatorname{atanh}(1/31)$, $\operatorname{atanh}(1/49)$ and
$\operatorname{atanh}(1/161)$, each by Zuniga's series where faster;
[](algorithms/floating-point/diophantine-reduction#machin)).
$2/\pi$ and $1/\log 2 - 1$ are derived from $\pi$ and $\log 2$.

## Euler's constant {#euler}

The Brent–McMillan algorithm [@BM1980] with the error bound of
[@BJ2013] gives, for a parameter $m$,
$$
\gamma = \frac{S_0}{I_0} - \frac{K_0}{I_0^2} - \log m + O(e^{-8m}),
\qquad
I_0 = \sum_{k \ge 0}\Bigl(\frac{m^k}{k!}\Bigr)^2, \quad
S_0 = \sum_{k \ge 0}\Bigl(\frac{m^k}{k!}\Bigr)^2 H_k,
$$
with $K_0 = I_0(2m)K_0(2m)$ evaluated by its asymptotic series and
$8m \ge p\log 2$ for $p$ bits. $S_0$ and $I_0$ come from one binary
splitting over about $4.97m$ terms in dual numbers.

::: proposition {#prop-dual}
Let $F(\varepsilon) = \sum_{k \ge 0} m^{2k}\big/\prod_{j=1}^k (j + \varepsilon)^2$
over $\mathbb{Q}[\varepsilon]/(\varepsilon^2)$. Then $I_0 = F(0)$ and
$S_0 = -\tfrac12 F'(0)$.
:::

::: proof
$\frac{d}{d\varepsilon}\prod_{j=1}^k (j + \varepsilon)^{-2}\big|_{\varepsilon=0}
= -2\sum_{j=1}^k \frac1j \prod_{j=1}^k j^{-2} = -2H_k/k!^2$.
:::

The splitting of $F$ carries a scalar $P$, a dual number
$D = D_0 + D_1\varepsilon = \prod(k + 1 + \varepsilon)$ (so that the
denominator is $D^2$) and $T = T_0 + T_1\varepsilon$, merged as
$$
T = P_1 T_2 + D_2^2 T_1, \qquad D = D_1D_2, \qquad P = P_1P_2,
$$
five full-size and five half-size products per merge. At the end, with
$X = T_0 + D_0^2$,
$$
\frac{S_0}{I_0} = \frac{2T_0D_1 - T_1D_0}{2D_0X}, \qquad
\frac{K_0}{I_0^2} = \frac{D_0^4\, T_K}{X^2 Q_K},
$$
where $T_K/Q_K$ is a second splitting over $2m$ terms at half the
precision. The main sum runs first and leaves only its quotient and two
half-precision values, so that the peak memory is that of the main sum.

The parameter $m$ is rounded up so that $\log m$ is a combination of a
few fast logarithms: over the primes $\{2\}$, $\{2, 3\}$, $\{2, 3, 5\}$
or $\{2, 3, 5, 7\}$, using $\log(9/8)$ and $\log(256/243)$; $\log(81/80)$,
$\log(32805/32768)$ and $\log(25/24)$; or $\log(2401/2400)$,
$\log(4375/4374)$, $\log(225/224)$ and $\log(64/63)$, each by Zuniga's
series. The splittings cost about $m(61 + \log_2 o)$ for the odd part $o$
of $m$ (powers of two in $m$ are free), so every $m \in [m_0, 1.5m_0]$
smooth over $\{2, 3, 5, 7\}$ is scored by this cost plus that of its
logarithms, and the best is taken (`_mp_real_const_euler_tune` forces
the set of primes). Below 3456 bits $\gamma$ comes from a static table.

## The cache of floors {#cache}

The constants are needed over and over by the argument reductions, at
precisions that grow during a computation. They are kept as exact
floors, which have the property that makes a single cache entry serve
every precision.

::: lemma {#lem-nest}
For $c \ge 0$ and $n \le N$,
$\lfloor cB^n \rfloor = \lfloor \lfloor cB^N \rfloor / B^{N-n} \rfloor$.
That is, the floor at $n$ limbs is the top part of the floor at $N$ limbs.
:::

::: proof
For any real $y$ and integer $q \ge 1$,
$\lfloor \lfloor y \rfloor / q \rfloor = \lfloor y/q \rfloor$; take
$y = cB^N$ and $q = B^{N-n}$.
:::

So a request at $n$ limbs reads the top limbs of the cached entry, and
its result does not depend on the precision at which the entry was
computed. A floor is obtained from a ball at a few guard limbs, retried
with more guard limbs if the ball contains a multiple of $B^{-n}$; for an
irrational constant this terminates. When a request exceeds the cached
precision $N$, the entry is recomputed at $\max(n + 5, 1.5N)$ limbs, so
that a sequence of growing requests costs at most a constant factor more
than its last one. The entry is a process-wide shared cache
([](thread_support)): one version at a time, replaced under a lock and
published with a release store, so readers never wait.

## Use in FLINT {#use}

- The ball and limb functions `mp_real_const_pi4`, ...,
  `_mp_real_const_inv_log2_frac` and `_mp_real_const_clear_cache`
  implement the formulas and the cache.
- `arb_const_pi`, `arb_const_log2`, `arb_const_e`, `arb_const_euler`,
  `arb_const_log10`, `arb_const_catalan`, `arb_const_apery`,
  `arb_const_zeta5` and `arb_hypgeom_gamma_fmpq` (at $1/3$ and $1/4$)
  call `_mp_real_const_arb`.
- The argument reductions of the elementary functions read $\pi/4$,
  $2/\pi$, $\log 2$ and $1/\log 2 - 1$ from the cache (static tables up
  to 64 limbs), and the AGM logarithm $\pi/4$ and $\log 2$.

::: functions mp_real_const_pi4 _mp_real_const_pi4 _mp_real_const_arb :::
