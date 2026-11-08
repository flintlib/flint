---
title: Series evaluation and the bit-burst algorithm
summary: How the elementary-function kernels evaluate the Taylor series of a reduced argument - rectangular splitting with tapered windows and integer coefficient chains, the chunked tangent series - and, at high precision, Brent's bit-burst algorithm and the table-free kernels that halve the argument.
complexity: About $2\sqrt{N/3}$ full products for $N$ terms by tapered rectangular splitting; $O(M(n)\log^2 n)$ by the bit-burst algorithm.
order: 12
---

After an argument reduction
([](algorithms/special-functions/elementary#reduction),
[](algorithms/special-functions/elementary-tables)), the kernels of the
elementary functions are left with the value of a function at a small
fixed-point argument $t < 2^{-r}$. This page describes how that value is
computed: by the Taylor series with rectangular splitting at moderate
precision, and by Brent's bit-burst algorithm at high precision, which
also serves the table-free kernels for single evaluations at huge
precision. Throughout, $B = 2^{\bits}$ and $n$ is a number of limbs.

## Rectangular splitting {#series}

A Taylor series $\sum_{k<N} a_k v^k$ with rational coefficients is
evaluated at a fixed-point $v$ by rectangular splitting [@Smi1989]
[@PS1973]: with blocks of $m$ terms,

$$
S_b = \sum_{i<m} a_{mb+i}\, v^i + v^m S_{b+1},
$$

the powers $v, \ldots, v^m$ are computed once ($m$ full products), every
term is one scalar multiply-add of a precomputed power by a one-limb
coefficient, and each block boundary costs one full product. The
coefficients $1/k!$ are generated as ratios of consecutive integers,
chained downwards so that one `mpn_divrem_1` closes each run of factors
that fills a limb; for $1/(2k + 1)$ the integers $D/(2k+1)$ for a group
of odd numbers with product $D$ are used, with one rescaling per group.

**Tapering.** Block $b$ is multiplied by $v^{mb} < 2^{-rmb}$ when it
enters the sum, so only its top $\lceil(\bits n + g - rmb)/\bits\rceil$
limbs matter. Computing each block in its shortened window reduces the
total cost to about $m$ products for the powers and $N/(3m)$ for the
block boundaries; with $m \approx \sqrt{N/1.5}$ this is about
$2\sqrt{N/3}$ full products, against $N/3$ for the tapered Horner scheme
(`_mp_real_series_tapered`), which is faster only below about a dozen
limbs. Errors committed in block $b \ge 1$ are damped by
$2^{-\min(g, rm)}$, and the error of the whole sum, a few units at its
own precision, by the final multiplication by $t$ or $t^2$ (for the
odd functions and $1 - \cos$), which is below $2^{-r}$. What remains is
the error of that last product, taken directly over the significant
limbs of its operands: one unit for the high product and up to one for
each operand truncated below it. The results are within 3.1 ulps for
$\sin$, $\sinh$, $\operatorname{atan}$ and $\operatorname{atanh}$, 3.6
ulps for $1 - \cos$ and $\cosh - 1$ (with the error of $t^2$ times a
sum near $1/2$), and 5 ulps for exp, whose sum is the result.

**The tangent** has coefficients $T_{k+2}/(2k+3)!$ involving the
tangent numbers, which are not ratios of small integers. Their leading
denominators are small, though: in `_mp_real_series_rs_tan`, the terms
are grouped into chunks whose common denominator $Q_j$ fits $j + 1$
limbs, and the sum is carried as $Q_j$ times the partial sum, so that a
term costs $j + 1$ scalar passes and moving to the next chunk one
division by $Q_j/Q_{j-1}$.

For the sizes on the hot paths (arguments below $2^{-32}$ or $2^{-64}$
after a reduction, up to a few dozen limbs), the series are generated
straight-line code; a generic engine making the same calls in the same
order serves the other sizes with identical results.

## Bit-burst evaluation {#bit-burst}

At high precision the Taylor series of the reduced argument has many
terms, and evaluating it with full-precision arithmetic costs
$O(\sqrt{N})$ full products at best. Brent's bit-burst algorithm
[@Bre1976b] instead splits the argument into slices of increasing
length,
$$
t = t_0 + t_1 + \cdots + t_K, \qquad t_k = \text{bits } 2^k\ldots 2^{k+1}
\text{ of } t \text{ (roughly)},
$$
and computes $\exp(t) = \prod_k \exp(t_k)$, where each $\exp(t_k)$ is
the sum of a series with a short rational argument, evaluated exactly
by binary splitting ([](algorithms/floating-point/hypergeometric-series)).
A slice of $2^k$ bits starting at bit $2^k$ needs about $p/2^k$ terms
whose integers grow by $2^k$ bits per term, so each slice costs
$O(M(p)\log p)$ and the whole $O(M(p)\log^2 p)$.

For the exponential the splitting computes
$T/(QB^{QE}) = \sum_k x^k/(k!\,B^{Dk})$ for a slice $x B^{-D}$. For the
sine and cosine, both are obtained from one tree over the series in
$y = x^2$:
$$
\sin x' = x'\,\frac{Q - B'}{Q}, \qquad 1 - \cos x' = \frac{A}{Q}, \qquad
-\frac{B'}{Q} = \sum_{k \ge 1}\frac{(-y)^k}{(2k+1)!}, \quad
-\frac{A}{Q} = \sum_{k \ge 1}\frac{(2k+1)(-y)^k}{(2k+1)!},
$$
using $1/(2k)! = (2k+1)/(2k+1)!$ so that both share the denominator
$Q = \prod (2i)(2i+1)$; the second track costs two extra products per
merge, half the work of a splitting of $\exp(ix)$ over all $2N$ terms.
Past a term threshold the cosine track is dropped and recovered from
the slice's denominator by one square root. The slice factors
$\exp(ix_k)$ are multiplied by `mp_real_mul_complex`.

These splittings have fixed midpoints, so every merge at a given level
multiplies by the same power $x^s$ of the argument; the powers are
tabulated, and their `fft_small` transforms are computed once per tree
for the entries used by at least two merges
([](algorithms/floating-point/hypergeometric-series#powtab)).

`_mp_real_exp_reduced` and `_mp_real_sin_cos_reduced` choose between the
series, one bit-burst step (the leading slice by binary splitting and
the rest by the series at the reduced size) and the full cascade, by
tuned thresholds that depend on the number of threads (a slice
parallelises where the series does not). The slices of a cascade run as
tasks on the available threads.

## Table-free evaluation {#notab}

A single evaluation at huge precision cannot amortise tables. The
table-free kernels halve the argument $h$ times, so that the reduced
argument is below $2^{-r}$ for a tuned depth $r$, evaluate the reduced
function by bit-burst, and undo the halving: the exponential by $h$
squarings, the sine and cosine by $h$ angle doublings on
$g = 1 - \cos$ alone,
$$
g(2\theta) = 2g(2 - g), \qquad \sin x = \sqrt{2g - g^2},
$$
with every intermediate a pure fraction and one square root at the end.

::: proposition {#prop-notab}
Squaring a value with relative error $\varepsilon$ gives relative error
at most $2\varepsilon + \varepsilon^2$; a doubling of $g$ with absolute
error $\delta$, for $g < 1/2$, gives an absolute error at most $4\delta +
2\delta^2$; and $\sqrt{2g - g^2}$ divides an absolute error of $g$ by
$2\sin x \cdot$ (to first order) $1/(1 - g)$.
:::

::: proof
$(y(1 + \varepsilon))^2 = y^2(1 + 2\varepsilon + \varepsilon^2)$. For
the doubling, $\frac{d}{dg} 2g(2 - g) = 4 - 4g \in (2, 4]$ and the
second derivative is $-4$. For the square root,
$\frac{d}{dg}\sqrt{2g - g^2} = (1 - g)/\sqrt{2g - g^2} = \cos x/\sin x$,
and $\cos x \le 1$.
:::

Hence $h$ squarings amplify the relative error by $2^h$, and $h$
doublings followed by the square root amplify the absolute error by at
most $2^{2h}/\sin x$; since any path that doubles at all has
$x \ge 2^{-z-1}$ for an argument with $z$ leading zero bits, $h + \log_2
n + 8$ resp. $2h + z + \log_2 n + 8$ guard bits keep the amplified
error below one output ulp, and the bound written is 128 ulps.

`mp_real_exp_notab_squaring` applies the same idea to balls, with $s$
ball squarings, and `mp_real_exp_agm` improves a starting value computed
this way by one Newton–Taylor step on the [AGM logarithm](algorithms/special-functions/elementary#agm).

## Use in FLINT {#use}

- `_mp_real_series_rs` and `_mp_real_series_rs_sin_cos` evaluate the
  reduced series by tapered rectangular splitting,
  `_mp_real_series_tapered` by the tapered Horner scheme and
  `_mp_real_series_rs_tan` the chunked tangent series.
- `_mp_real_exp_reduced` and `_mp_real_sin_cos_reduced` choose between
  the series and the bit-burst algorithm; `_mp_real_exp_notab` and
  `_mp_real_sin_cos_notab` are the table-free kernels, and
  `mp_real_exp_notab_squaring` and `mp_real_exp_agm` the table-free ball
  exponentials.

::: functions _mp_real_series_rs _mp_real_series_tapered _mp_real_series_rs_tan _mp_real_exp_reduced _mp_real_sin_cos_reduced _mp_real_exp_notab _mp_real_sin_cos_notab mp_real_exp_notab_squaring mp_real_exp_agm :::
