---
title: Double-word and multi-double arithmetic
summary: Error-free transformations, floating-point expansions of two to four doubles, their renormalisation, multiplication, long division and square roots, and how dfloat turns the dropped terms into rigorous ball radii.
complexity: $O(N^2)$ double operations for a product of $N$-term expansions, $O(N^2)$ for a long division with tapered remainders; about $N^2$ fused multiply-adds dominate.
order: 85
---

A *floating-point expansion* represents a real number as an unevaluated
sum $x = x_0 + x_1 + \cdots + x_{N-1}$ of machine doubles, each much
smaller than the previous one. With $N = 2$ (double-word or
double-double arithmetic) this gives about 106 bits of precision, with
$N = 4$ (quad-double) about 212 bits, using only hardware floating-point
operations. The idea goes back to Møller [@Mol1965] and Dekker
[@Dek1971]; general expansions are due to Priest [@Pri1991] and
Shewchuk [@She1997], the classical quad-double library is that of Hida,
Li and Bailey [@HLB2001], and tight error bounds for double-word
operations were proved by Joldes, Muller and Popescu [@JMP2017]. A
comprehensive treatment is in [@MBD2018, chapters 4 and 14].

FLINT implements expansions with $N = 1, \ldots, 4$ in the `dfloat`
module. Its kernels are designed for ball arithmetic: every operation
computes the exact terms of the result with error-free transformations,
keeps $N$ of them, and returns the magnitude of everything it dropped,
which becomes the radius of a ball.

Throughout, doubles are IEEE 754 binary64 numbers with round-to-nearest
($\operatorname{RN}$), $u = 2^{-53}$ is the unit roundoff, and
$\operatorname{ulp}(x)$ is the weight of the last bit of the double
$x$. For any real $t$ such that $\operatorname{RN}(t)$ is a normal
double, $\operatorname{RN}(t) = t(1 + \delta)$ with $|\delta| \le u$.

## Error-free transformations {#eft}

The sum and the product of two doubles are not doubles in general, but
their rounding errors are, and can be computed exactly.

::: lemma Sterbenz {#lem-sterbenz}
If $a$ and $b$ are doubles with $b/2 \le a \le 2b$, then $a - b$ is a
double, so it is computed exactly [@Ste1974].
:::

::: algorithm TwoSum and FastTwoSum {#alg-twosum}
**TwoSum**$(a, b)$ [@Knu1997; @Mol1965]: $s \gets \operatorname{RN}(a + b)$,
$b' \gets \operatorname{RN}(s - a)$,
$e \gets \operatorname{RN}(\operatorname{RN}(a - \operatorname{RN}(s - b')) + \operatorname{RN}(b - b'))$.

**FastTwoSum**$(a, b)$ [@Dek1971]: $s \gets \operatorname{RN}(a + b)$,
$e \gets \operatorname{RN}(b - \operatorname{RN}(s - a))$.
:::

::: theorem {#thm-twosum}
If $a + b$ does not overflow, TwoSum returns $s = \operatorname{RN}(a + b)$
and $e$ with $s + e = a + b$ exactly, using six operations and no
branch. FastTwoSum gives the same result in three operations if the
exponent of $a$ is at least that of $b$ (in particular if
$|a| \ge |b|$) or $a = 0$. In both cases $|e| \le \operatorname{ulp}(s)/2$.
:::

::: proof
For FastTwoSum with $|a| \ge |b|$, we first show that $z = s - a$ is
computed exactly. If $a$ and $b$ have the same sign, then
$|a| \le |s| \le 2|a|$ with $s$ of the sign of $a$, so this follows from
[](#lem-sterbenz). If they have opposite signs and $|b| \ge |a|/2$, then
$a + b$ is exact by the same lemma, so $s = a + b$ and $z = b$. If they
have opposite signs and $|b| < |a|/2$, then $|a|/2 \le |s| \le |a|$
(rounding is monotone and $|a|/2$ is a double), and again
[](#lem-sterbenz) applies. Now $b - z = (a + b) - s$ is the rounding
error of a nearest-rounded sum of two doubles, which is itself a double
(its magnitude is at most $\operatorname{ulp}(s)/2$ and it is a multiple
of the ulp of the operand with the smaller exponent), so it is also
computed exactly.
TwoSum reconstructs the same quantities without knowing which operand is
larger; see [@MBD2018, §4.3] for the complete proof. Sums are exact in
the subnormal range, so no underflow condition is needed.
:::

::: algorithm TwoProd {#alg-twoprod}
**TwoProd**$(a, b)$: $p \gets \operatorname{RN}(ab)$,
$e \gets \operatorname{RN}(ab - p)$ computed with one fused multiply-add.
:::

::: theorem {#thm-twoprod}
If $ab$ does not overflow and $|ab| \ge 2^{-969}$ or $ab = 0$, TwoProd
returns $p + e = ab$ exactly.
:::

::: proof
The exact product $ab$ has at most 106 significant bits and is a
multiple of $\operatorname{ulp}(a)\operatorname{ulp}(b)$. Its rounding
error $ab - p$ is at most $\operatorname{ulp}(p)/2$ in magnitude and is
a multiple of the same quantity, so it fits in 53 bits; the condition
on $|ab|$ ensures that this multiple is not below the subnormal range.
An FMA rounds the exact $ab - p$ once, which is exact for a double.
:::

Without an FMA, TwoProd can be computed with Dekker's product, which
splits each operand into two 26-bit halves (Veltkamp splitting) and
costs 17 operations [@Dek1971]; `_d_vec_dot_thrice` uses this form.
In C, all of these transformations break if the compiler contracts a
separate multiplication and addition into an FMA, or evaluates
intermediate results in extended precision (x87). The `dfloat` module
is therefore compiled with `-ffp-contract=off` and with SSE2 arithmetic,
and `dfloat_is_supported` reports whether these requirements hold.

## Expansions and renormalisation {#renormalisation}

An expansion $(x_0, \ldots, x_{N-1})$ is *nonoverlapping* if
$|x_{i+1}| \le \operatorname{ulp}(x_i)/2$ for all $i$ with $x_i \ne 0$,
and *canonical* if in addition zeros occur only at the end. A canonical
expansion is obtained from a real number $x$ by rounding term by term,
$x_0 = \operatorname{RN}(x)$, $x_1 = \operatorname{RN}(x - x_0)$, and so
on; it is unique up to ties.

Arithmetic produces a list of terms, ordered by decreasing magnitude
but overlapping, which must be converted to $N$ nonoverlapping
components.

::: algorithm VecSum (bottom-up renormalisation) {#alg-vecsum}
**Input:** doubles $t_0, \ldots, t_{m-1}$, roughly decreasing in
magnitude.

**Output:** $t_0', \ldots, t_{m-1}'$ with the same exact sum.

1. For $i = m - 2, m - 3, \ldots, 0$: $(t_i, t_{i+1}) \gets \text{TwoSum}(t_i, t_{i+1})$.
:::

Since each step is error-free, the sum is preserved exactly, and after
the pass each $t_{i+1}$ is the rounding error of a sum whose rounded
value is $t_i$, so $|t_{i+1}| \le \operatorname{ulp}(t_i)/2$. The result
is not necessarily canonical: a cancellation in the sum that produced
$t_{i+1}$ can leave $t_{i+2}$ larger than $t_{i+1}$. FLINT calls this
*weak* normalisation. It is all that the error analysis needs, and it is
what the `dfloat` kernels produce; a second, top-down pass of TwoSums
(VecSumErr) followed by a check makes the expansion canonical
(`dN_canonicalise`, or every result in a context with `DFLOAT_STRONG`).
For $N = 2$ weak normalisation is already canonical. The general
renormalisation `_dfloat_renorm` of an arbitrary number of terms adds a
branchy pass (VecSumErrBranch of [@JMPT2016]) that removes zero
components before nonzero ones.

## Addition and multiplication {#mul}

The kernels organise the exact terms of a result by *order*: a term of
order $k$ has magnitude about $u^k$ times that of the head.

::: algorithm Addition of $N$-term expansions {#alg-add}
1. For $i = 0, \ldots, N-1$: $(s_i, e_i) \gets \text{TwoSum}(x_i, y_i)$.
2. $S_0 \gets s_0$. For $k = 1, \ldots, N-1$: accumulate
   $S_k \gets s_k + e_{k-1} + (\text{residuals of order } k-1)$ with
   TwoSums, keeping their residuals as terms of order $k$.
3. Renormalise $S_0, \ldots, S_{N-1}$ with [](#alg-vecsum).
4. The dropped part is $e_{N-1}$ plus the residuals of order $N-1$.
:::

::: algorithm Multiplication of $N$-term expansions {#alg-mul}
1. For $i + j < N$: $(h_{ij}, l_{ij}) \gets \text{TwoProd}(x_i, y_j)$.
2. $S_0 \gets h_{00}$. For $k = 1, \ldots, N-1$: accumulate with TwoSums
   $S_k \gets \sum_{i+j=k} h_{ij} + \sum_{i+j=k-1} l_{ij} + (\text{residuals of order } k-1)$.
3. Renormalise $S_0, \ldots, S_{N-1}$ with [](#alg-vecsum).
4. The dropped part consists of the residuals of order $N - 1$, the low
   parts $l_{ij}$ with $i + j = N - 1$, and the products $x_i y_j$ with
   $i + j \ge N$, which are not computed exactly and are bounded by
   $|\operatorname{RN}(x_i y_j)|$.
:::

The products need $N(N+1)/2$ TwoProds, each a multiplication and an
FMA. Squaring uses the symmetry $h_{ij} = h_{ji}$ (doubling is exact).
The product with a single double needs $N$ TwoProds.

::: theorem {#thm-tracked}
Let $\mathit{res}$ be the renormalised result and $E$ the sum of the
absolute values of the dropped terms of [](#alg-add) or [](#alg-mul).
If no product of nonzero components has magnitude below $2^{-969}$ and
nothing overflows, then $|x \circ y - \mathit{res}| \le E$, where
$\circ$ is the exact operation.
:::

::: proof
Every step is an error-free transformation by [](#thm-twosum) and
[](#thm-twoprod), so $x \circ y$ equals the sum of the kept and the
dropped terms exactly, except for the products of order $\ge N$, whose
exact values satisfy $|x_i y_j| \le |\operatorname{RN}(x_i y_j)|(1 + u)$;
the factor $1 + u$ and the rounding errors in forming $E$ itself are
absorbed by the slack factor of [](#lem-slack).
:::

For the plain formats, $E$ is not computed, and the error is of the
order of $u^N$ times the magnitude of the result, with small constants
[@HLB2001; @JMP2017]. In the ball formats, $E$ is added to the radius.
When the inputs are not exact, tracking the last-order residuals exactly
costs as much as the product itself for $N = 4$, so the ball products
then accumulate the last order with plain additions and bound their
rounding errors by $u$ times the accumulated magnitudes; this loses
only a constant factor in a radius that is nonzero anyway, and keeps
exact products of exact balls exact.

## Radius arithmetic {#radius}

The radius of a ball is a double computed in round-to-nearest
arithmetic from nonnegative quantities, and is made an upper bound by a
single multiplication at the end.

::: lemma {#lem-slack}
Let $r$ be a nonnegative real number obtained from nonnegative doubles
by at most $2^7$ additions and multiplications, and $\tilde r$ the
result of performing the same operations in round-to-nearest arithmetic,
where no operation has an inexact subnormal result. Then
$r \le \operatorname{RN}(\tilde r (1 + 2^{-45}))$.
:::

::: proof
Each operation on nonnegative numbers satisfies
$t \le \operatorname{RN}(t)/(1 - u)$, so by induction
$r \le \tilde r (1 - u)^{-2^7} \le \tilde r \exp(2^7 u (1 + u)) \le \tilde r (1 + 2^{-46} + 2^{-90})$.
On the other hand
$\operatorname{RN}(\tilde r (1 + 2^{-45})) \ge \tilde r (1 + 2^{-45})(1 - u) \ge \tilde r (1 + 2^{-45} - 2^{-52})$,
which is larger. Sums of doubles in the subnormal range are exact, and
products and quotients that could underflow are bumped up by $2^{-1070}$
(`_dfloat_rad_mul`), which preserves the inequality.
:::

The ball operations then use the standard propagation bounds. For the
product, $|x'y' - xy| \le |x| r_y + |y| r_x + r_x r_y$ for
$|x' - x| \le r_x$ and $|y' - y| \le r_y$, with $|x|$ bounded by the sum
of the absolute values of the components. For the quotient:

::: proposition {#prop-div}
Let $q$ be any number and $|y| > r_y$. For all $x'$ and $y'$ with
$|x' - x| \le r_x$ and $|y' - y| \le r_y$,
$$\left|\frac{x'}{y'} - q\right| \le \frac{|x - qy| + r_x + |q| r_y}{|y| - r_y}.$$
:::

::: proof
$x'/y' - q = (x' - qy')/y'$, with $|x' - qy'| \le |x - qy| + |x' - x| + |q||y' - y|$
and $|y'| \ge |y| - r_y$.
:::

The residual $x - qy$ for the computed quotient $q$ is enclosed with the
tracked kernels, so that an exact quotient of exact balls gets radius
zero. For the reciprocal square root, with $e = 1 - x r^2$ for the
computed $r$, the series
$(1 - e)^{-1/2} = \sum_k \binom{2k}{k} (e/4)^k$ has positive
coefficients and $\sum_{k \ge 2} \binom{2k}{k} 4^{-k} |e|^{k-2} \le 1/2$
for $|e| \le 1/4$, so that $|x^{-1/2} - r| \le |r|(|e|/2)(1 + |e|)$.

## Division and square root {#division}

::: algorithm Long division with tapered remainders {#alg-div}
**Input:** expansions $x$, $y$ with $y_0 \ne 0$.

**Output:** $q = q_0 + \cdots + q_{N-1} \approx x/y$.

1. $r \gets x$, $v \gets \operatorname{RN}(1/y_0)$.
2. For $k = 0, \ldots, N-1$:
   1. $q_k \gets \operatorname{RN}(r_0 v)$ (for the last digit, $\operatorname{RN}(r_0/y_0)$).
   2. If $k < N - 1$: $(p, e) \gets \text{TwoProd}(q_k, y_0)$,
      $d_0 \gets r_0 - p$ (exact), and
      $r \gets (d_0 - e) + (r_1 + \cdots) - q_k (y_1 + \cdots)$, kept with
      $N - k - 1$ components.
3. Renormalise $q_0, \ldots, q_{N-1}$.
:::

The digit $q_k$ is within a relative $2u + u^2$ of $r_0/y_0$, so that
$p = \operatorname{RN}(q_k y_0)$ is within a factor 2 of $r_0$ and
$d_0 = r_0 - p$ is exact by [](#lem-sterbenz): the cancellation of the
heads loses no information, and the next remainder absorbs the error of
the digit. After step $k$ the remainder is about $u^{k+1}|x|$ and is only
needed to absolute precision $u^N |x|$, so it is kept with fewer
components (*tapered* precision); the remainder update is a three-term
sum of expansions of the size of the new remainder. A
Newton–Karp–Markstein division [@KM1997] (a reciprocal at half
precision, one product and a correction) turned out to be slower at
every $N$, because the steps of the long division are cheap products of
a double with the divisor and the divisions are off the critical path.

The square root is the Newton iteration $s \gets s + (x - s^2)/(2s)$
from $s = \operatorname{RN}(\sqrt{x_0})$, doubling the number of correct
components per step. The residual $x - s^2$ is computed as
$(x_0 - s_0^2) - (s - s_0)(s + s_0) + (x_1 + \cdots)$ with the head
exact (TwoProd and [](#lem-sterbenz)), and the correction is a long
division producing only the digits it needs. The reciprocal square root
iterates $r \gets r + r(1 - x r^2)/2$ in the same way.

## Use in FLINT {#use}

- `dfloat`: the plain formats `dN` use [](#alg-add), [](#alg-mul),
  [](#alg-div) and the Newton iterations; the ball formats `dNb` add the
  dropped terms of [](#thm-tracked) to the radius with
  [](#lem-slack), and use [](#prop-div) for `dNb_div` and the bound
  above for `dNb_rsqrt`. The kernels are written once
  (`src/dfloat/kernels.inc`) and instantiated for doubles and for
  four-lane SIMD vectors, which the vector functions `_dN_vec_mul`,
  `_dN_vec_dot`, etc. use. The error-free transformations are exposed
  as the macros `DFLOAT_TWO_SUM`, `DFLOAT_FAST_TWO_SUM` and
  `DFLOAT_TWO_PROD`.
- The elementary functions of `dfloat` are built on these kernels; see
  [](algorithms/floating-point/dfloat-elementary).
- `_dfloat_platt_smk` computes the moments of Platt's multi-evaluation
  in plain double-double arithmetic with the algorithms DWTimesDW1 and
  AccurateDWPlusDW of [@JMP2017], using their proven error bounds a
  priori.
- `_d_vec_dot_thrice` in `d_vec` is the compensated dot product of
  Ogita, Rump and Oishi [@ORO2005], built on TwoSum and Dekker's product.

::: functions d2_mul d2b_mul d2_div _d2_mul_tracked _dfloat_renorm d2_canonicalise :::
