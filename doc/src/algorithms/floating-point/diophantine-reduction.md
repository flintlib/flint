---
title: Multi-prime argument reduction
summary: Reduction of exp and sin/cos arguments by integer combinations of logarithms of primes (or arguments of Gaussian primes) found from precomputed integer relations, with the logarithms computed by Machin-type formulas whose terms are summed by Zuniga's fast series.
complexity: Precomputation of $k$ constants in $O(k\,M(n)\log^2 n)$; each evaluation $O(M(n)\log^2 n)$ with a smaller constant than without reduction.
order: 50
---

The bitwise reduction of [](algorithms/floating-point/elementary-functions#bitwise)
needs a table of $r$ logarithms $\log(1 + 2^{-i})$ at the full
precision, each a series of its own; at 600 limbs and $r = 768$ this
table already takes noticeable time to build. Van der Hoeven and the
author [@HJ2024] observed that a much smaller table suffices if one
reduces by integer combinations of a few constants: the logarithms of
the first primes for the exponential, and the arguments of the first
Gaussian primes for the trigonometric functions. `mp_real` uses this
*diophantine* reduction from 600 to 65536 limbs.

## The exponential {#exp}

Let $p_0 = 2, p_1 = 3, \ldots, p_{k-1}$ be the first $k$ primes. Given
$x \in [0, 1)$, the reduction chooses integers $c_j$ so that
$$
t = x - \sum_{j<k} c_j \log p_j
$$
is tiny, and then
$$
\exp(x) = 2^{c_0}\,\frac{P}{Q}\,\exp(t), \qquad
P = \prod_{j \ge 1,\, c_j > 0} p_j^{c_j}, \quad
Q = \prod_{j \ge 1,\, c_j < 0} p_j^{-c_j}.
$$
Beyond $\exp(t)$, which the reduced series evaluates
([](algorithms/floating-point/elementary-functions#bit-burst)), the work
is computing the integers $P$ and $Q$, one multiplication by $P$ and one
division by $Q$, all of moderate size.

**Relations.** The coefficients come from a precomputed table of
integer relations among the logarithms: rows $d_i \in \Z^k$ with small
entries and
$$
\sum_j d_{ij}\log p_j = \varepsilon_i, \qquad
|\varepsilon_1| > |\varepsilon_2| > \cdots,
$$
found by lattice reduction and independent of the precision
(`_mp_real_rel_table`; tables for common $k$ are stored in the library).

**Descent.** After the free step by $\log 2$ (whose power is a shift),
each row in turn subtracts the nearest integer multiple of
$\varepsilon_i$ from the residual and adds the same multiple of $d_i$
to the coefficient vector. The descent stops when the weight
$$
w(c) = \sum_{j \ge 1} |c_j| \log_2 p_j,
$$
which bounds $\log_2(PQ)$, would exceed a budget. Larger budgets give
deeper reductions (a smaller $t$, hence a cheaper series) and larger
$P$ and $Q$, which have about half the weight in bits each. The
residual is tracked in double precision and resynchronised every eight
rows from a short fixed-point evaluation of $x - \sum c_j\log p_j$.
Rounding to nearest leaves $t$ of either sign, and a negative $t$ is
lifted by one more relation with positive $\varepsilon$, which leaves
$0 < t \le \varepsilon$.

**Prime powers.** $P$ and $Q$ are computed by a vector exponentiation
on the bits of the exponent vector,
$$
\Pi(c) = \Pi(\lfloor c/2 \rfloor)^2 \prod_{j:\, c_j \text{ odd}} p_j,
$$
one squaring per exponent bit, with the product over the odd
exponents a single limb. The total costs about one squaring at the final
size. Products beyond the working precision are carried truncated.

::: proposition Error of the reduction {#prop-error}
Let the cached logarithms be floors $\ell_j$ with $\log p_j - B^{-n} <
\ell_j \le \log p_j$, and let $\hat t = x - \sum_j c_j \ell_j$ be
computed exactly. Then $|\hat t - t| < (\sum_j |c_j|)\,B^{-n}$, and the
relative error of $\exp(\hat t)$ as an approximation of $\exp(t)$ is
below $(\sum_j |c_j| + 1) B^{-n}$ for $n$ large enough.
:::

::: proof
$\hat t - t = \sum_j c_j(\log p_j - \ell_j)$, each term below $|c_j|
B^{-n}$ in absolute value; and
$|\exp(\hat t - t) - 1| \le |\hat t - t|(1 + |\hat t - t|)$.
:::

The multiplication by $P$ is exact and the division by $Q$ a floor, and
the factor $2^{c_0}P/Q = \exp(x - t) < e$ amplifies the errors by less
than 3, so one guard limb makes the total error at most 3 ulps
(`_mp_real_exp_diophantine`). The default uses 13 primes and a weight
budget equal to the precision in bits.

## The trigonometric functions {#trig}

For the sine and cosine, the constants are $\alpha_0 = \pi/2$ and
$\alpha_j = 2\arg\pi_j$ for the first nonreal Gaussian primes
$\pi_1 = 1 + i$, $\pi_2 = 1 + 2i$, $\pi_3 = 2 + 3i$, ... With
$x = c_0\pi/2 + \sum_j c_j\, 2\arg\pi_j + t$,
$$
e^{ix} = i^{c_0}\, e^{it}\,\frac{A^2}{N}, \qquad
A = \prod_j \pi_j^{c_j}, \quad N = |A|^2,
$$
where $\pi_j^{c_j}$ means $\bar\pi_j^{|c_j|}$ for $c_j < 0$, since
$\arg(\pi_j^2) = 2\arg\pi_j$ and $|\pi_j^2| = N(\pi_j)$. The norm $N$ is a
rational integer, so the normalisation is one reciprocal of an integer
and two multiplications; the tangent is the ratio of the two parts of
$e^{it}A^2$ and needs no normalisation at all. The reduced sine and
cosine dominate the cost and need a deeper reduction than the
exponential to be efficient, so the default uses 32 Gaussian primes
and a budget of four times the precision.

## Computing the logarithms {#machin}

The constants $\log p_j$ and $\arg\pi_j$ are computed together, to the
working precision, by Machin-type formulas
$$
\log p_i = \frac{1}{\mathit{den}}\sum_j c_{ij}\operatorname{atanh}\frac{1}{x_j},
\qquad
\arg\pi_i = \frac{1}{\mathit{den}}\sum_j c_{ij}\operatorname{atan}\frac{1}{x_j},
$$
with square coefficient matrices and large arguments $x_j$ (the $x_j$
for which $(x_j \pm 1)$, resp. $x_j^2 + 1$, are smooth over the primes
of the set). Sets exist for 3 to 32, 40 and 48 values; they were found
with the program at https://github.com/fredrik-johansson/machin. A set
of $k$ terms costs about as much as one series of $\sum_j n/\log_2
x_j$ terms. The terms are independent, so they run on separate threads,
costliest first, and the cached values are exact floors read off the
resulting balls.

### Zuniga's series {#zuniga}

Each term of a Machin-type formula for logarithms is the logarithm of a
ratio of smooth numbers: $2\operatorname{atanh}(1/x) = \log(u/v)$ with
$u/v = (x + 1)/(x - 1)$ in lowest terms. Zuniga's Ramanujan-type series
[@Zun2025] gives, for rationals $u/v > 1$,
$$
\log\frac{u}{v} = \sum_{n=1}^{\infty} \rho^n
\frac{\alpha n + \beta}{\gamma\, n(2n - 1)}
\frac{(1)_n(1/2)_n}{(1/6)_n(5/6)_n},
\qquad \rho = \frac{(u - v)^6}{108\, u^2v^2(u + v)^2},
$$
$$
\alpha = -2(u + v)(u^2 - 14uv + v^2)(u^2 + 4uv + v^2), \quad
\beta = (u + v)^3(u^2 - 8uv + v^2), \quad
\gamma = 2(u - v)^5.
$$
For $u/v = 2$ this is the series with $\rho = 1/3888$ that the module
uses for $\log 2$. In terms of $x$,
$\rho = 4/(27x^2(x^2 - 1)^2)$, so one term gains
$6\log_2 x + \log_2(27/4)$ bits, as much as three terms of the Taylor
series of $\operatorname{atanh}(1/x)$ and 2.75 bits more. The same
series at $u, v = t \pm id$ gives $\log(u/v) = 2i\operatorname{atan}(d/t)$
[@Zun2025b], an alternating series in the real quantities $s = 2t$ and
$p = t^2 + d^2$ with $\rho = -4d^6/(27t^2p^2)$.

In binary splitting, the cost per bit is governed by how fast the
integers grow compared with the precision gained per term: about
$\log_2(2n)$ bits against $2\log_2 x$ for the Taylor series, and about
$2\log_2 n + 1$ against $6\log_2 x + 2.75$ for Zuniga's series, which
is two thirds as much relatively, and much less for small $x$. The terms
are three times longer, which costs more at low precision and for large
$x$; the module uses Zuniga's series for arguments of $b$ bits from a
precision of about $b^2/2$ limbs and the Taylor series
(`_mp_real_atanh_frac_bsplit`, `_mp_real_atan_frac_bsplit`) below.
Measured on the logarithms of the first 2 to 8 primes, the speedup is
1.2 to 1.5 at $10^6$ to $10^8$ bits; for 32 to 48 primes (arguments of
62 to 87 bits) it is a few percent. For the angles of the first 4
Gaussian primes the speedup is 1.7 to 1.9.

## Precomputation against evaluation {#tradeoff}

Compared with the bitwise reduction, the precomputation is much cheaper
(one Machin-type evaluation for all $k$ constants instead of $r$ series
at the full precision), at the price of a shallower reduction and hence
a slower evaluation. The bitwise kernels are therefore used up to 600
limbs and the diophantine ones above; beyond 65536 limbs even the
Machin-type constants do not pay for a single evaluation, and the
table-free kernels take over. The program `profile/p-diophantine` prints
the trade-off.

## Use in FLINT {#use}

- `_mp_real_exp_diophantine`, `_mp_real_sin_cos_diophantine` and
  `_mp_real_tan_diophantine` implement the reductions, with tunable
  variants taking the number of primes and the budget; they are the
  kernels of `mp_real_exp_bits`, `mp_real_sin_cos_bits` and
  `mp_real_tan_bits` from 600 to 65536 limbs, and the forward functions
  of the Newton–Taylor logarithm and arctangent.
- `_mp_real_rel_table` holds the relations, `_mp_real_machin_table` the
  Machin-type sets, and `_mp_real_log_primes_vec` and
  `_mp_real_atan_gauss_vec` compute the constants;
  `_mp_real_log_ratio_zuniga` and `_mp_real_atan_ratio_zuniga` sum
  Zuniga's series through `mp_real_hypgeom_series`.
- `arb` shares the relation tables and Machin-type sets
  (`arb_exp_arf_log_reduction`, `arb_sin_cos_arf_atan_reduction`,
  `arb_log_primes_vec_bsplit`), and `mp_real_const_log10` and Euler's
  constant use the logarithms of 2, 3, 5 and 7 from the same sets.

::: functions _mp_real_exp_diophantine _mp_real_sin_cos_diophantine _mp_real_rel_table _mp_real_machin_table :::
