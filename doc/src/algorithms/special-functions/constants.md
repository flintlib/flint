---
title: Mathematical constants
summary: The formulas and algorithms behind the constants of FLINT (π, e, log 2, log 10, Euler's constant, Catalan's constant, ζ(3), ζ(5), Γ(1/3), Γ(1/4), Khinchin's, Glaisher's and the reciprocal Fibonacci constant) - binary splitting of hypergeometric series, the dual-number splitting for Euler's constant, sums of zeta values and, at very high precision, an approximate functional equation for Khinchin's constant - with their error bounds, and how the values are cached.
complexity: $O(M(n) \log^2 n)$ for the hypergeometric constants by binary splitting ($\Gamma(1/4)$ from $2^{16}$ limbs in $O(M(n)\log n)$ plus the cost of $\pi$) and $O(M(n) \sqrt{n} \log n)$ for the reciprocal Fibonacci constant; Khinchin's constant in $O(n^3/\log^3 n)$ word operations by sums of zeta values and in $O(\sqrt n\, M(n) \log^2 n)$ by the approximate functional equation; Glaisher's constant (below about $10^5$ to $3.5 \cdot 10^5$ bits, depending on the number of threads) uses sums of many full-precision terms and is much slower, and costs $\tilde O(n^{1.5})$ at higher precision.
order: 20
---

Most of the mathematical constants of FLINT are sums of *hypergeometric
series*, series whose term ratio is a rational function of the index,
and are evaluated by binary splitting. The series for $\pi$, $e$,
$\log 2$, $\log 10$, $\gamma$, Catalan's constant, $\zeta(3)$, $\zeta(5)$,
$\Gamma(1/3)$ and $\Gamma(1/4)$, and Khinchin's constant, are summed by
the `mp_real` module, which `arb_const_pi` and the other `arb_const_`
functions call through `_mp_real_const_arb`; the same values, kept as
exact floors, serve the argument reductions of the elementary functions
([](algorithms/special-functions/elementary)). The formulas were taken
from the collection of y-cruncher [@Yee2025], and for each constant the
formula measured fastest in this implementation was chosen.

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

In `mp_real` the splitting is `mp_real_hypgeom_series`
([](algorithms/floating-point/hypergeometric-series)), which sums a
series given by integer polynomials $P$, $Q$, $R$ in the format of
y-cruncher's formula files, with rigorous tail bounds and threads. Since
it gives $\bits(n - 1)$ correct bits, the constants call it with one limb
more than they need.

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

## $\pi$ {#pi}

$\pi$ is computed by the Chudnovsky series [@CC1988]

$$
\frac{1}{\pi} = 12 \sum_{k=0}^{\infty} \frac{(-1)^k (6k)!\,(13591409 +
545140134 k)}{(3k)!\,(k!)^3\, 640320^{3k + 3/2}},
$$

which adds about 47 bits (14 digits) per term. Methods based on the
arithmetic-geometric mean [@Sal1976] [@Bre1976agm] have a better
asymptotic complexity, $O(M(n) \log n)$, but are about three times slower
in practice. `arb_const_sqrt_pi` and `arb_const_log_sqrt2pi` are computed
from `arb_const_pi` with `arb_sqrt` and `arb_log`.

The Chudnovsky series in the format of `mp_real_hypgeom_series` is
$S = \bigl((c_Q + c_P\sum_{k \ge 1}\frac{P(k)}{Q(k)}\prod_{j<k}\frac{R(j)}{Q(j)})/c_D\bigr)^{-1}$
with $c_P = 1$, $c_Q = 13591409$, $c_D = 4270934400$,
$R(k) = (6k - 5)(2k - 1)(6k - 1)$, $P(k) = (13591409 + 545140134k)R(k)$
and $Q(k) = -(640320^3/24)k^3$, and $\pi = S/\sqrt{10005}$; the example
on [](mp_real/hypgeom) runs it.

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
three series evaluated in parallel (each by Zuniga's series where faster;
[](algorithms/special-functions/diophantine-reduction#machin)). The logarithms of more primes, used
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
[@BJ2013].

Below 3456 bits $\gamma$ comes from a static table. Otherwise, writing
$m$ for the parameter $n$, $S_0$ and $I_0$ come from one binary
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
the set of primes).

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

## $\Gamma(1/3)$ and $\Gamma(1/4)$ {#gamma}

$\Gamma(1/3) = (810^{1/4}\pi/X)^{1/3}$ with Guillera's series $X$ (2023),
and $\Gamma(1/4) = (\pi^6/(322S^4))^{1/8}$ with Ebisu's series $S$
(2016); the series and $\pi$ run on separate threads. For
$\Gamma(1/4)$ above $2^{16}$ limbs, the formula
$$
\Gamma(1/4) = \sqrt{\frac{(2\pi)^{3/2}}{\operatorname{agm}(1, \sqrt 2)}}
$$
is used instead, with the AGM ($O(M(n)\log n)$) computed in parallel
with $\pi$.

## Khinchin's constant {#khinchin}

### Zeta expansion

Up to about $1.5 \cdot 10^6$ digits,
Khinchin's constant $K$ is computed from the formula [@BBC1997]

$$
\log K = \frac{1}{\log 2} \left[ -\sum_{k=2}^{K_0-1} \log
\left(\frac{k-1}{k}\right) \log \left(\frac{k+1}{k}\right) +
\sum_{n=1}^{\infty} \frac{\zeta(2n, K_0)}{n} \sum_{k=1}^{2n-1}
\frac{(-1)^{k+1}}{k} \right],
$$

where $K_0 \ge 2$ is a free parameter. The finite sum has a minus sign:
each product of logarithms is negative. `arb_const_khinchin` calls
`mp_real_const_khinchin`, which takes $K_0$ between 6 and 12, evaluates the
finite sum as a quadratic form in the logarithms of the primes up to
$K_0$, and computes the $\zeta(2n, K_0)$ by direct summation with recycled
powers for large $n$ and from exact Bernoulli numbers for small $n$
([below](#khinchin-algorithm)).

::: proposition Truncation error {#prop-khinchin}
If the series over $n$ is truncated after $n = M$, the error inside the
brackets is at most $1/K_0^{2M}$.
:::

::: proof
The inner alternating sums lie in $(0, 1]$, so the remainder is bounded
in absolute value by $\sum_{n > M} \zeta(2n, K_0)/n \le \sum_{n > M}
\zeta(2n, K_0)$. By the integral bound,

$$
\sum_{n=M+1}^{\infty} \zeta(2n, K_0) = \sum_{n=M+1}^{\infty}
\sum_{k=0}^{\infty} (k + K_0)^{-2n} \le \sum_{n=M+1}^{\infty}
\left(K_0^{-2n} + \int_0^{\infty} (t + K_0)^{-2n}\, dt\right)
= \sum_{n=M+1}^{\infty} \frac{1}{K_0^{2n}} \left(1 + \frac{K_0}{2n -
1}\right),
$$

and since $K_0/(2n - 1) \le K_0$, this is at most $\sum_{n > M} (K_0 + 1)
K_0^{-2n} = \frac{1}{K_0^{2M}(K_0 - 1)} \le \frac{1}{K_0^{2M}}$.
:::

For an error of $2^{-p}$ it therefore suffices to take $M \ge p / (2
\log_2 K_0)$.

### The algorithm {#khinchin-algorithm}

With
$f(x) = -\log(1-x)\log(1+x) = \sum_{n \ge 1} c_n x^{2n}$,
$c_n = A_n/n$ and $A_n = 1 - \tfrac12 + \tfrac13 - \cdots + \tfrac{1}{2n-1}$,
$$
S = \log K \log 2 = \sum_{k \ge 2} f(1/k)
  = \sum_{k=2}^{K_0-1} f(1/k) + \sum_{n \ge 1} c_n\,\zeta(2n, K_0),
$$
the classical sum over $\zeta(2n) - 1$ [@BBC1997] being the case
$K_0 = 2$. Each cell $(k, n)$ of the double sum $\sum_k\sum_n c_n k^{-2n}$
is covered in one of three ways.

*Rows* $k < K_0$ use $f(1/k) = (\ell_k - \ell_{k-1})(\ell_{k+1} - \ell_k)$
with $\ell_m = \log m$, an integer combination of the logarithms
$\lambda_j$ of the primes up to $K_0$, which come from one set of
Machin-type series ([](algorithms/special-functions/diophantine-reduction#machin)).
The sum of the rows is a quadratic form $\sum C_{ij}\lambda_i\lambda_j$
with small integer $C$, so it costs one full product per prime rather
than one per row.

*Columns* $n \ge n_1$ sum $\zeta(2n, K_0) = \sum_{k \ge K_0} k^{-2n}$ directly,
in fixed point with $W$ fraction limbs, over the $k$ with
$k^{-2n} \ge B^{-W}$. The powers are recycled from column to column,
downward: $k^{-2n}B^{W+e}$ is kept as an integer with $e$ extended limbs
and stepped by one `mpn_mul_1` by $k^2$, several times faster than the
division by $k^2$ an upward recurrence needs. A downward step multiplies
the absolute error by $k^2$ but keeps the relative error, so a power is
created (binary powering and one division) with the relative accuracy
needed at the bottom of its range, its bottom limb is dropped whenever
its error count passes $2^{20}B$ (so that the truncations, one unit of
the next limb each, cost a negligible relative amount over thousands of
drops), and the range is cut into up to three segments, each created
afresh, so that the products shorten with the accuracy actually
needed. Only $k \in [K_0, 2K_0)$ and the odd $k \ge 2K_0$ are summed:
the even $k \ge 2K_0$ contribute $4^{-n}\zeta(2n, K_0)$, so the partial
sum is divided by $1 - 4^{-n}$ by a few shifted additions.

Streaming every power, and the column sum, once per column is bound by
the memory bandwidth once the powers outgrow the L2 cache, and the more
so the more threads share it. The columns are therefore taken in batches of 16 with one
partial sum per column, and the powers eight at a time by blocks of 256
limbs of the column sums: each block of a power is stepped through the
16 columns while in cache, and the blocks of the 16 sums serve the
eight powers. The drops of bottom limbs depend only on the error counts,
so they are scheduled before the products, which are those of the
column-by-column recurrence; the carries between blocks are kept per
power and column. This roughly halves the time per limb and step on one
thread for powers of $2 \cdot 10^5$ limbs, at a cost of 16 column sums
of $W$ limbs per task. The subtracted powers $k < K_0$ of the Bernoulli
columns (below) are summed the same way.

*Columns* $n < n_1$ use $\zeta(2n, K_0) = \zeta(2n) - 1 - \sum_{k<K_0} k^{-2n}$
with $\zeta(2n) = |B_{2n}|(2\pi)^{2n}/(2(2n)!)$ and the exact Bernoulli
numbers of the reverse iterator `bernoulli_rev` [@Blo2009], which recycles odd
powers downward at the precision of $B_{2n}$, below the working
precision for these $n$. The factor $t_n = (2\pi)^{2n}/(2n)!$ would
cost one full product per column; instead the columns are taken in
blocks $[lo, hi]$ of $b \approx \min(2\sqrt{n_1}, 128)$ and computed in
units of $t_{hi}$,
$$
\frac{\zeta(2m)}{t_{hi}} = \frac{|B_{2m}|\,(2m+1)\cdots(2\,hi)}{2}\,(2\pi)^{-2(hi-m)},
$$
with the powers $(2\pi)^{-2j}$, $j \le b$, computed once: a column then
costs one product of the integer by a full-precision power, and a block
four full products.

The weights $c_n$ are applied by summation by parts: with $a_1 = 1$,
$a_i = -1/((2i-2)(2i-1))$, so that $A_n = a_1 + \cdots + a_n$, and
$R_i = \sum_{lo \le n \le hi,\, n \ge i} X_n/n$,
$$
\sum_{n=lo}^{hi} c_n X_n = A_{lo-1}R_{lo} + \sum_{i=lo}^{hi} a_i R_i,
$$
two single-limb divisions and two additions per column instead of a
product by the real $A_n$. Since $\zeta(2n, K_0) \approx K_0^{-2n}$, the
suffix sums of the direct columns only occupy the low limbs, and their
cost tapers off with the column. All fixed-point errors are counted in
ulps and bounded rigorously, including the omitted terms
$\sum_{k \ge k_1} k^{-2n} \le k_1^{-2n}(1 + k_1/(2n-1))$ and the columns
beyond the last one.

A direct column $n$ costs about $(K(n) - K_0)/2$ products of single limbs
by numbers of up to $W$ limbs, with $K(n) = B^{W/2n}$, and a Bernoulli
column about two products of the size $s_n$ of $B_{2n}$ (in the
iterator) and one of size $s_n$ by $W$ limbs; $n_1$ is chosen where
$K(n_1)$ is between $2^3$ and $2^{13}$ (measured), growing with the
precision up to $3 \cdot 10^6$ bits. Beyond, $n_1$ is chosen so that
$B_{2n_1}$ has about as many bits as the working precision, so that the
direct column $n_1$ takes about as many powers as the iterator does
there ($K(n_1) \approx 2^{16.8}$ at $10^7$ digits); otherwise the
Bernoulli columns, whose cost grows like the cube of the precision with
$n_1$ fixed, would take nearly all the time. For
$n$ limbs, the full products cost $O(n\,M(n)/\log n)$ in all, and the
power sums of the iterator and the direct columns $O(n^3/\log^3 n)$ word
operations, which dominate at very high precision.

The sums being linear, the direct columns are split into groups of $k$
(two per thread), the subtracted powers $\sum_{k<K_0} k^{-2n}$ into ranges
of columns (with $A_{lo-1} = H_{2lo-3} - H_{lo-2}$ from exact harmonic
numbers) and the Bernoulli columns into ranges of $n$, each with its own
summation by parts, which run as independent tasks together with the
rows. Each range of Bernoulli columns has its own iterator, whose power
table takes about $(2\,hi)^2/24$ bits ([](algorithms/number-theory/bernoulli#reverse)),
which dominates the memory at high precision: with one range per thread,
some five times the table of the top range at once (about 5 GB at
$2 \cdot 10^6$ digits). The ranges are therefore run one per thread only
while their tables total below 1 GiB; beyond, as many ranges run at once
as stay below it (at least one), each with the threads left over, used
inside the iterator and for the columns: a group of 16 columns takes its
16 numbers from the iterator, then their products by the powers of
$(2\pi)^{-2}$ in parallel, the prefix sums $R_i$, and the divisions by
$(2i-2)(2i-1)$ in parallel, the $\sum a_i$ of a block being summed
exactly. The direct columns, the subtracted powers and the rows then
follow. `_mp_real_const_khinchin_tune` forces $K_0$, $n_1$ and the
number of Bernoulli ranges run at once.

### The approximate functional equation {#khinchin-afe}

From `MP_REAL_KHINCHIN_AFE_CUTOFF` limbs (about $1.5 \cdot 10^6$
digits), `mp_real_const_khinchin` uses `_mp_real_const_khinchin_afe`,
which needs $O(\sqrt n)$ evaluations of cost $O(M(n)\log^2 n)$ each
instead of the cubic power sums above. It starts from the dilogarithm
series of [@BBC1997, §3],
$$
\log K \log 2 = \log^2 2 + \operatorname{Li}_2(-\tfrac12)
  + \tfrac12 \sum_{k \ge 2} (-1)^k \operatorname{Li}_2(4/k^2),
$$
and the Laplace transform $\operatorname{Li}_2(4/K) = \int_0^\infty
e^{-Ku} G(u)\,du$ for $K > 4$, with $G(u) = E(4u)/u$ and $E(x) =
\sum_{n \ge 1} x^n/(n \cdot n!)$. The sum over $k \ge 3$ is then an
integral of $G$ against a theta function. Splitting the integral at
$u_0 = 4^{-a}$ and applying the transformation of $\theta_4$ on
$(0, u_0]$ gives, for every $u_0 > 0$,
$$
\log K \log 2 = \log^2 2 + \operatorname{Li}_2(-\tfrac12) + \frac{\pi^2}{12}
  + \frac12\Bigl(\sum_{k \ge 3} (-1)^k I_k + \sum_{m \ge 0} J_m + M\Bigr)
$$
with
$$
I_k = \int_{u_0}^\infty e^{-k^2 u} G(u)\,du, \quad
J_m = \sqrt\pi \int_0^{u_0} u^{-1/2} e^{-\pi^2(m+1/2)^2/u} G(u)\,du, \quad
M = \int_0^{u_0} \bigl(e^{-u} - e^{-4u} - \tfrac12\bigr) G(u)\,du.
$$
For $p$ bits, about $\sqrt{p \log 2/u_0}$ terms $I_k$ and
$\sqrt{p u_0 \log 2}/\pi$ terms $J_m$ are needed; the default
$u_0 = 1/16$ balances their costs. $\operatorname{Li}_2(-\tfrac12)$
comes from $\operatorname{Li}_2(1/9)$ and $\log 3$ by Landen's identity,
and $M$ is a rational series.

This is a generalization of the approximate functional equation of the
Dirichlet $L$-functions
([](algorithms/special-functions/zeta#l-functions)), which is the same
construction for the summand $k^{-s}$: there the function $G$ is a power
of $u$, so that every term of the two series is a single incomplete
gamma function. Here $G$ is an entire function, and the terms $I_k$ and
$J_m$ are themselves series, which can still be evaluated by binary
splitting since their coefficients are rational.

*The $k$-terms.* With $K = k^2$ and $x = Ku_0$,
$I_k = e^{-x}\sum_{n \ge 1} (4/K)^n e_n(x)/n^2$, $e_n$ the truncated
exponential series. The terms satisfy a triangular $2 \times 2$
recurrence with rational coefficients, summed by matrix binary
splitting; the factors $e^{-k^2 u_0}$ of consecutive $k$ are propagated
by multiplications by $e^{-(2k+1)u_0}$. The tail after $N$ terms is
bounded through the Poisson weights in $e^{-x}e_n(x)$.

*The dual terms.* $J_m = \sqrt\pi\,\nu_{3/2}(c)$ with
$c = \pi^2(m+\tfrac12)^2$ and the moments
$\nu_t(c) = \int_0^{u_0} u^{-t} e^{-c/u} E(4u)\,du$.
Integration by parts relates $\nu_{t+1}$ to $\nu_t$, to the similar
moments $Y_t$ of $e^{4u}$ and $C_t$ of 1, and to the boundary values
at $u_0$, which carry the factor $e^{-c/u_0}$; and $-d/dc$ shifts $t$ by
one. So $\nu_{3/2}$ is a holonomic function of $c$ of order four whose
only finite singularity is $c = 0$, and it is evaluated by the bit-burst
algorithm: starting values at a 64-bit dyadic $c_0$ near $c$, from a
shift series in terms of
$\Gamma(\tfrac12, c_0/u_0)$ (by `mp_real_gamma_upper_half_frac`) and
$e^{-c_0/u_0}$, then Taylor steps to dyadic points with $64 \cdot 4^j$
fractional bits converging to $c$. The parts of the shift series cancel
by up to $e^{c_0/u_0}$ individually but their weighted sums only by
about $e^{4\sqrt{c_0}}$, so $O(\sqrt{c_0})$ guard bits suffice. All
Taylor remainders, the tails of the starting values and the final
distance to $c$ are bounded rigorously through upper incomplete gamma
functions; the bounds carry the factor $e^{-c_0/u_0}$, so that the late
stages are computed only to the significant bits of $J_m$, which is of
order $e^{-c/u_0}$.

*Matrix binary splitting.* The $k$-terms, the dual terms (starting
values and Taylor steps), $M$ and the continued fraction for
$\Gamma(\tfrac12, y)$ all use one internal engine
(`_mp_real_linrec_eval`) for recurrences
$X_{l+1} = (A_l X_l + f_l B_l v)/q_l$ with sums
$S \mathrel{+}= (S^A_l X_l + f_l S^B_l v)/q_l$ and a scalar geometric
forcing $f_{l+1} = \varphi_l f_l/q_l$, given by integer leaf matrices.
The range is cut into chunks whose exact products (integer matrices over
a common denominator) have about half the working precision in bits.
Within a chunk, binary splitting combines the two halves in one
augmented product that shares the transforms of each entry among the
products (`fmpz_mat_mul_fft_small`, since `fmpz_mat_mul` does not choose
the transformed product for such small matrices) when the entries have
at least about $1.2 \cdot 10^4$ bits and the inner dimension $d+g+s$ is
at least 6; for the narrow products of the $k$-terms, separate products
were faster at all sizes measured. Each chunk's matrix is then applied
to the state vector $(X, f v)$: above the chunks there are only
matrix-vector products. The state is kept scaled by the product of the
denominators so far, which is divided out once at the end. For large
entries the application is one truncated transformed product
(`fmpz_mat_mul_fft_small_trunc`, which transforms each state entry once
for all the rows and accumulates the rows in the transformed domain),
each state entry converted to fixed point with its own exponent, so
that it keeps its relative accuracy as in ball arithmetic; for small
entries it is done in ball arithmetic.

*Parallelism and memory.* The terms are independent tasks. The dual
terms are the most expensive and take about four times the memory of a
$k$-term, so they start first and at most `maxdual` of them run at once
(by default half the threads, until the last chunk of $k$-terms has
started; then any thread may take one), while the other threads take
chunks of consecutive $k$-terms, starting each chunk with $e^{-k^2 u_0}$ and
$e^{-(2k+1) u_0}$ as powers of a shared $e^{-u_0}$ (about $3 \log_2 k$
multiplications instead of two exponentials). The exact products of
large chunks also fork onto idle threads. Each finished term is added to one fixed-point integer
with an integer error count, so no term is kept and the result does not
depend on the order in which the terms finish. The memory is then that
of the running terms, linear in the precision times the number of
threads, and dominated by the transforms of the shared-transform
products, which fft_small limits to
`flint_fft_small_max_transformed_ring_size` per product (blocking the
products beyond). At $10^6$ digits, a dual term takes about 370 numbers of the
working precision with all transforms resident and a $k$-term about 90.
The user can reduce memory in exchange for a slight increase
in running time by decreasing `flint_fft_small_max_transformed_ring_size`.

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

## Caching {#cache}

Each `arb_const_` function keeps a thread-local cached ball. A request at
precision $p$ above the cached precision recomputes the constant at $p +
32$ bits; later requests at precision at most $p$ only round the cached
value, to nearest with a radius of half a unit in the last place plus
the cached error. $\pi$ and $\log 2$ up to about 4000 bits are read from
the static floors of $\pi/4$ and $\log 2$ (64 limbs on 64-bit machines),
and $\gamma$ up to 3456 bits from a static table, correctly rounded to
nearest with a radius of half a unit in the last place. The caches are
freed by `flint_cleanup`.

Inside `mp_real` the constants $\pi/4$, $\log 2$, $2/\pi$ and
$1/\log 2 - 1$ (the last two derived from the first two) are needed over
and over by the argument reductions, at precisions that grow during a
computation. They are kept as exact
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

::: functions arb_const_pi arb_const_e arb_const_log2 arb_const_log10 arb_const_euler arb_const_catalan arb_const_apery arb_const_zeta5 arb_const_khinchin arb_const_glaisher arb_const_reciprocal_fibonacci arb_const_sqrt_pi arb_const_log_sqrt2pi mp_real_const_pi4 _mp_real_const_pi4 _mp_real_const_arb :::

The complex versions in `acb` (for example `acb_const_pi`) call these
functions, and `arb_hypgeom_gamma_fmpq` uses $\Gamma(1/3)$ and
$\Gamma(1/4)$ ([the gamma function at rational
arguments](algorithms/special-functions/gamma#rational)). The argument
reductions of the elementary functions read $\pi/4$, $2/\pi$, $\log 2$
and $1/\log 2 - 1$ from the cache of floors (static tables up to 64
limbs), and the AGM logarithm $\pi/4$ and $\log 2$.
`_mp_real_const_khinchin_afe` implements the approximate functional
equation for Khinchin's constant, with `mp_real_gamma_upper_half_frac`
for its starting values.

The tutorial [](tutorials/numerics/constants) computes $\pi$ and other
constants to many digits in complete programs.
