---
title: Mathematical constants
summary: The formulas and algorithms behind arb_const_pi, arb_const_e, arb_const_log2, arb_const_log10, arb_const_euler, arb_const_catalan, arb_const_apery, arb_const_zeta5, arb_const_khinchin, arb_const_glaisher and arb_const_reciprocal_fibonacci, with their error bounds, and how the values are cached.
complexity: $O(M(n) \log^2 n)$ for the hypergeometric constants by binary splitting and $O(M(n) \sqrt{n} \log n)$ for the reciprocal Fibonacci constant; Khinchin's constant and Glaisher's constant (below about $10^5$ to $3.5 \cdot 10^5$ bits, depending on the number of threads) use sums of many full-precision terms and are much slower; Glaisher's constant at higher precision costs $\tilde O(n^{1.5})$.
order: 20
---

Most of the mathematical constants provided by `arb` are sums of
*hypergeometric series*, series whose term
ratio is a rational function of the index, and are evaluated by binary
splitting. The series for $\pi$, $e$, $\log 2$, $\log 10$, $\gamma$,
Catalan's constant, $\zeta(3)$ and $\zeta(5)$ are summed by the `mp_real`
module, and the `arb_const_` functions call it through
`_mp_real_const_arb`; the formulas were chosen as the fastest among
those collected by y-cruncher [@Yee2025]. The implementation side
(splitting trees, threads, the shared cache of exact floors) is described
in [](algorithms/floating-point/constants); the sections below give the
mathematics and the remaining constants.

## Binary splitting {#binary-splitting}

Consider a sum $S = \sum_{k=0}^{N-1} a(k)/b(k) \cdot
\prod_{j<k} p(j)/q(j)$ with polynomials $a, b, p, q$ with integer
coefficients. For $0 \le u < v \le N$ define $P(u, v) = \prod_{j=u}^{v-1}
p(j)$, $Q(u, v) = \prod_{j=u}^{v-1} q(j)$, and $T(u, v) = Q(u, v)
\sum_{k=u}^{v-1} (\ldots)$, the partial sum scaled to an integer. Then for
$u < w < v$

$$
P(u, v) = P(u, w) P(w, v), \quad Q(u, v) = Q(u, w) Q(w, v), \quad
T(u, v) = T(u, w) Q(w, v) + P(u, w) T(w, v)
$$

(with the obvious extension for $a/b$), so the exact sum is obtained by a
balanced recursion and one final division $T/Q$.

::: proposition {#prop-bsplit}
If the series converges geometrically, so that $N = O(n)$ terms give $n$
bits, and the integers $P(u, v), Q(u, v), T(u, v)$ have $O((v - u)
\log v)$ bits, binary splitting computes $S$ to $n$ bits in $O(M(n \log
n) \log n)$ bit operations, which is $O(M(n) \log^2 n)$.
:::

::: proof
At each of the $O(\log N)$ levels of the recursion, the numbers multiplied
have a total size of $O(N \log N) = O(n \log n)$ bits, and multiplication
is superadditive, so each level costs $O(M(n \log n))$. The final
division costs $O(M(n))$.
:::

The series below add between 2 and 47 bits per term. The constant
factor is what distinguishes the formulas, and it is governed by the
number of bits per term and the size of the coefficients
[@HP1998] [@BZ2011].

## Caching {#cache}

Each `arb_const_` function keeps a thread-local cached ball. A request at
precision $p$ above the cached precision recomputes the constant at $p +
32$ bits; the result is then rounded to $p$ bits, and later requests at
precision at most $p$ only round the cached value. $\pi$, $\log 2$ and
$\gamma$ are moreover read from static tables (of $\pi/4$ and $\log 2$ to
4608 bits, of $\gamma$ to 3456 bits) below these sizes, correctly rounded
to nearest with a radius of half a unit in the last place. The caches are
freed by `flint_cleanup`.

## $\pi$ {#pi}

$\pi$ is computed by the Chudnovsky series [@CC1988]

$$
\frac{1}{\pi} = 12 \sum_{k=0}^{\infty} \frac{(-1)^k (6k)!\,(13591409 +
545140134 k)}{(3k)!\,(k!)^3\, 640320^{3k + 3/2}},
$$

which adds about 47 bits (14 digits) per term. In the hypergeometric form
used by `mp_real`, $\pi = S / \sqrt{10005}$ with
$R(k) = (6k - 5)(2k - 1)(6k - 1)$, $P(k) = (13591409 + 545140134 k)
R(k)$ and $Q(k) = -(640320^3/24) k^3$. Methods based on the
arithmetic-geometric mean [@Sal1976] [@Bre1976agm] have a better
asymptotic complexity, $O(M(n) \log n)$, but are about three times slower
in practice. `arb_const_sqrt_pi` and `arb_const_log_sqrt2pi` are computed
from `arb_const_pi` with `arb_sqrt` and `arb_log`.

## $e$ {#e}

$e = \sum_{k \ge 0} 1/k!$ is summed by binary splitting, with $P = 1$
and $Q(k) = k$; this is the simplest hypergeometric series and its cost
is dominated by the final division.

## Logarithms of integers {#log}

$\log 2$ is computed by Zuniga's series [@Zun2025]

$$
\log 2 = \frac{1}{2160} \sum_{k=0}^{\infty} (1497 + 1794 k)
\prod_{j=1}^{k} \frac{j (2j - 1)}{216 (6j + 1)(6j + 5)},
$$

which adds 11.9 bits per term. (Earlier versions used another series of
Zuniga [@Zun2023b], with terms decreasing like $3888^{-k}$.) $\log 10 =
\log 2 + \log 5$ uses the Machin-type formulas for the first three primes
in terms of
$$
\operatorname{atanh}(1/31), \quad \operatorname{atanh}(1/49), \quad
\operatorname{atanh}(1/161),
$$
for example $\log 10 = 46 \operatorname{atanh}(1/31) + 34
\operatorname{atanh}(1/49) + 20 \operatorname{atanh}(1/161)$, with the
three series evaluated in parallel. The logarithms of more primes, used
by the elementary functions ([](algorithms/special-functions/elementary)),
come from larger sets of such formulas (`arb_log_primes_vec_bsplit`).

## Euler's constant {#euler}

Euler's constant $\gamma$ is computed by the Brent-McMillan formula
[@BM1980] [@MPFR2012]

$$
\gamma = \frac{S_0(2n) - K_0(2n)}{I_0(2n)} - \log n,
$$

where $n$ is a free parameter and

$$
S_0(x) = \sum_{k=0}^{\infty} \frac{H_k}{(k!)^2} \left(\frac{x}{2}\right)^{2k},
\quad I_0(x) = \sum_{k=0}^{\infty} \frac{1}{(k!)^2}
\left(\frac{x}{2}\right)^{2k}, \quad
2x I_0(x) K_0(x) \sim \sum_{k=0}^{\infty} \frac{[(2k)!]^3}{(k!)^4 8^{2k}
x^{2k}}.
$$

The first two series are summed together up to $k = N - 1$ with $N \ge
\alpha n + 1$, where $\alpha \approx 4.9706257595442318644$ satisfies
$\alpha(\log \alpha - 1) = 3$, and the asymptotic series up to $k = 2n -
1$. With these parameters the error is bounded by $24 e^{-8n}$
[@BJ2013]. The current implementation sums $S_0$ and $I_0$ in one binary
splitting over dual numbers, and chooses $n$ smooth over $\{2, 3, 5, 7\}$
so that $\log n$ is cheap; see [](algorithms/floating-point/constants).

## Catalan's constant {#catalan}

Catalan's constant is computed by the series of Pilehrood and Pilehrood
[@PP2010]

$$
G = \frac{1}{64} \sum_{k=1}^{\infty} \frac{256^k (580 k^2 - 184 k + 15)
\,(2k)!^3\,(3k)!^2}{k^3 (2k - 1)\,(6k)!^2},
$$

which adds about 5.4 bits per term and was measured to be the fastest
available formula.

## $\zeta(3)$ and $\zeta(5)$ {#zeta}

Apéry's constant $\zeta(3)$ is computed by a hypergeometric series found
by Zuniga [@Zun2023],

$$
\zeta(3) = \frac{1}{48} \sum_{k=1}^{\infty} \frac{(-1)^{k-1}
P(k)}{k^5 (2k-1)^3 (3k-1)(3k-2)(4k-1)(4k-3)(6k-1)(6k-5)
\binom{5k}{k}\binom{5k}{2k}\binom{9k}{4k}\binom{10k}{5k}\binom{12k}{6k}},
$$

with $P(k)$ a polynomial of degree 11 given in the source. $\zeta(5)$ is
computed from Z.-W. Sun's identity

$$
\zeta(5) = \frac{3 S + 56 \pi^2 \zeta(3)}{540}, \quad
S = \sum_{n=1}^{\infty} \left(-\frac{1}{27}\right)^n \frac{(560 n^4 -
640 n^3 + 408 n^2 - 136 n + 17)\, n!^2}{(1/3)_n (2/3)_n\, n^5 (2n-1)^4},
$$

with the three parts $S$, $\zeta(3)$ and $\pi$ computed in parallel when
threads are available.

## Khinchin's constant {#khinchin}

Khinchin's constant $K_0$ is computed from the formula [@BBC1997]

$$
\log K_0 = \frac{1}{\log 2} \left[ -\sum_{k=2}^{N-1} \log
\left(\frac{k-1}{k}\right) \log \left(\frac{k+1}{k}\right) +
\sum_{n=1}^{\infty} \frac{\zeta(2n, N)}{n} \sum_{k=1}^{2n-1}
\frac{(-1)^{k+1}}{k} \right],
$$

where $N \ge 2$ is a free parameter (`arb_const_khinchin` uses $N \approx
p^{0.35}$). The finite sum has a minus sign: each product of logarithms
is negative.

::: proposition Truncation error {#prop-khinchin}
If the series over $n$ is truncated after $n = M$, the error inside the
brackets is at most $1/N^{2M}$.
:::

::: proof
The inner alternating sums lie in $(0, 1]$, so the remainder is bounded
in absolute value by $\sum_{n > M} \zeta(2n, N)/n \le \sum_{n > M}
\zeta(2n, N)$. By the integral bound,

$$
\sum_{n=M+1}^{\infty} \zeta(2n, N) = \sum_{n=M+1}^{\infty}
\sum_{k=0}^{\infty} (k + N)^{-2n} \le \sum_{n=M+1}^{\infty}
\left(N^{-2n} + \int_0^{\infty} (t + N)^{-2n}\, dt\right)
= \sum_{n=M+1}^{\infty} \frac{1}{N^{2n}} \left(1 + \frac{N}{2n -
1}\right),
$$

and since $N/(2n - 1) \le N$, this is at most $\sum_{n > M} (N + 1)
N^{-2n} = \frac{1}{N^{2M}(N - 1)} \le \frac{1}{N^{2M}}$.
:::

For an error of $2^{-p}$ it therefore suffices to take $M \ge p / (2
\log_2 N)$. The Hurwitz zeta values are computed as $\zeta(2n, N) =
\zeta(2n) - \sum_{k=1}^{N-1} k^{-2n}$, with `arb_zeta_ui` and a vector of
powers $k^{-2n}$ that is updated from one $n$ to the next.

## Glaisher's constant {#glaisher}

The Glaisher-Kinkelin constant $A = \exp(1/12 - \zeta'(-1))$ is computed
directly from this formula, with $\zeta'(-1)$ from `acb_dirichlet_zeta_jet`
(see [](algorithms/special-functions/zeta)). Up to a cutoff
(`_acb_dirichlet_l_jet_use_afe`: $3.5 \cdot 10^5$ bits with one thread,
about $8 \cdot 10^4$ bits with 8 threads), this uses Euler-Maclaurin
summation directly at $s = -1$; the
reflection formula is not used for $-4 < \operatorname{Re}(s) < 0$, since
the arithmetic of Euler-Maclaurin summation is cheaper at $s = -1$ than at
$s = 2$ (by a factor 1.3–1.7 for $\zeta'(-1)$).

Euler-Maclaurin summation needs $O(p)$ terms and $O(p)$ Bernoulli numbers
at precision $p$, which costs $\tilde O(p^2)$ time and $O(p^2)$ bits of
memory (several gigabytes at $p = 4 \cdot 10^5$). Above the cutoff,
`acb_dirichlet_zeta_jet` instead computes $\zeta(-1+x)$ to
$O(x^2)$ with the approximate functional equation
(`acb_dirichlet_l_fmpq_afe_jet`), a sum of $O(\sqrt{p})$ incomplete
gamma function jets evaluated by binary splitting, in time
$\tilde O(p^{1.5})$ and with $O(p)$ memory.

## The reciprocal Fibonacci constant {#fibonacci}

$\sum_{n \ge 1} 1/F_n$ is computed from Gosper's series [@Gos1974],
corrected in [@Arn2012],

$$
\sum_{n=1}^{\infty} \frac{1}{F_n} = \sum_{n=0}^{\infty} \frac{(-1)^{n(n-1)/2}
(F_{4n+3} + (-1)^n F_{2n+2})}{F_{2n+1} F_{2n+2} L_1 L_3 \cdots L_{2n+1}},
$$

where $L_n = 2F_{n-1} + F_n$ is a Lucas number. The truncation error
after $N \ge 1$ terms is bounded by $(1/\phi)^{N^2}$, so $N \approx 1.2
\sqrt{p}$ terms give $p$ bits. The series is not hypergeometric, but
factoring out the multiplicative recurrence for $L_1 L_3 \cdots$ turns it
into a product of $O(\sqrt{p})$ matrices with $O(\sqrt{p})$-bit entries,
which binary splitting evaluates in quasilinear time.

## Use in FLINT {#use}

::: functions arb_const_pi arb_const_e arb_const_log2 arb_const_log10 arb_const_euler arb_const_catalan arb_const_apery arb_const_zeta5 arb_const_khinchin arb_const_glaisher arb_const_reciprocal_fibonacci arb_const_sqrt_pi arb_const_log_sqrt2pi :::

The complex versions in `acb` (for example `acb_const_pi`) call these
functions. The constants $\Gamma(1/3)$ and $\Gamma(1/4)$, used by
`arb_hypgeom_gamma_fmpq`, are also computed by `mp_real`; see
[](algorithms/floating-point/constants) and
[the gamma function at rational arguments](algorithms/special-functions/gamma#rational).

The tutorial [](tutorials/numerics/constants) computes $\pi$ and other
constants to many digits in complete programs.
