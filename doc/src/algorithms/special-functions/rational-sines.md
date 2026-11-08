---
title: Sines and cosines of rational multiples of π
summary: Evaluating cos(πR/Q), sin(πR/Q) and tan(πR/Q) as algebraic numbers - closed forms and half-angle steps by square roots, a Chebyshev step of arbitrary order with hypergeometric corrections that needs no division or square root, and the choice between these and the general kernels in mp_real, arb and the functions of πx at dyadic x.
complexity: $O(M(n)\log Q)$ for the Chebyshev iteration, about $k + 1$ square roots for the half-angle steps at $Q = 2^k Q_0$.
order: 15
---

The values $\cos(\pi R/Q)$ are algebraic of degree at most $Q$, and at
high precision it is cheaper to compute one as an algebraic number than
to evaluate the cosine with the general kernels of
[](algorithms/special-functions/elementary).
`mp_real_sin_cos_pi_ui_div_ui` (and `mp_real_tan_pi_ui_div_ui`) chooses
by the denominator and the precision between closed forms with square
roots, half-angle steps, a Chebyshev step whose order grows with the
precision, and the general kernel on $\pi R/Q$; `arb_sin_cos_pi_fmpq`,
`arb_tan_pi_fmpq` and the functions of $\pi x$ at a dyadic $x$ with a
small denominator use it where it is faster.

## The Chebyshev step {#chebyshev}

::: theorem Chebyshev step {#thm-chebyshev}
Let $\varphi^* = \pi R/Q \in (0, \pi/2)$ in lowest terms,
$x^* = \cos\varphi^*$, $\varepsilon = (-1)^R$, and let $x_0 = \cos\varphi_0$
with $\varphi_0 = \varphi^* + \delta$, $|Q\delta| < \pi/2$. Put
$C = T_Q(x_0)$, $A = x_0 T_Q(x_0) - T_{Q+1}(x_0)$, $a = \varepsilon A$,
$y = 1 - C^2$ and $b = 1/Q$. Then
$$
x^* = x_0 + x_0\, G(y) + \frac{a}{Q} F(y),
$$
$$
G(y) = {}_2F_1\bigl(-\tfrac b2, \tfrac b2; \tfrac12; y\bigr) - 1, \qquad
F(y) = {}_2F_1\bigl(\tfrac{1-b}{2}, \tfrac{1+b}{2}; \tfrac32; y\bigr).
$$
:::

::: proof
Since $T_k(\cos\varphi) = \cos(k\varphi)$, $C = \cos(Q\varphi^* +
Q\delta) = \varepsilon\cos(Q\delta)$, and
$x_0 T_Q(x_0) - T_{Q+1}(x_0) = \cos\varphi_0\cos(Q\varphi_0) -
\cos((Q+1)\varphi_0) = \sin\varphi_0 \sin(Q\varphi_0)$, so
$a = s_0\sigma$ with $s_0 = \sin\varphi_0$ and $\sigma = \sin(Q\delta)$,
and $y = \sigma^2$. Then $x^* = \cos(\varphi_0 - \delta) = x_0\cos\delta +
s_0\sin\delta$, and with $\delta = b\operatorname{asin}\sigma$ (valid on
the principal branch, $|Q\delta| < \pi/2$) the classical expansions
$\cos(b\operatorname{asin}\sigma) = {}_2F_1(-\frac b2, \frac b2; \frac12;
\sigma^2)$ and $\sin(b\operatorname{asin}\sigma) = b\sigma\,{}_2F_1(\frac{1-b}2,
\frac{1+b}2; \frac32; \sigma^2)$ give the claim, $s_0\sin\delta$
becoming $s_0\sigma\, b F(y) = (a/Q)F(y)$.
:::

No division or square root enters: $s_0$ appears only through $a$ and
$y$. The pair $(T_k(x_0), T_{k+1}(x_0))$ is computed by a Montgomery
ladder on $V_k = 2T_k$, $V_{2k} = V_k^2 - 2$, $V_{2k+1} = V_kV_{k+1} -
V_1$, one squaring and one product per bit of $Q$. The coefficients of
$G$ and $F$ decrease in magnitude from $1/(2Q^2)$ resp. 1, so for
$|y| \le 1/2$ the tails after $y^J$ are below $|y|^J/Q^2$ resp.
$2|y|^J$; in $y/(4Q^2)$ both have integer coefficients over common
denominators, so the series cost a few products and one division by a
word each. With $x_0$ correct to $e$ bits, $\sigma \approx Q2^{-e}$, and
summing to the target precision makes a step of any order $r$ from $x_0$
at about $(p + (r-1)\log_2 Q)/r$ bits; the order used is 12. The
recursion runs in ball arithmetic and checks the branch condition
$|Q\delta| < \pi/2$ from the rigorous radius of the previous level. The
cost is $O(M(n)\log Q)$, against the $O(M(n)\log^2 n)$ of the general
kernel; the iteration is used for denominators up to a bit length that
grows with the precision (none below 48 limbs, 12 bits at 512 limbs, 48
at 65536 limbs, and all denominators beyond). Of $\cos(\pi t/q)$ and
$\sin(\pi t/q) = \cos(\pi(q - 2t)/(2q))$ the one with the smaller
denominator is iterated and the other is $\sqrt{1 - x^{*2}}$.

## Square roots and half angles {#half-angle}

The denominators up to 6 have closed forms ($\cos(\pi/3) = 1/2$,
$\cos(\pi/4) = \sqrt 2/2$, $\cos(\pi/5) = (\sqrt 5 + 1)/4$,
$\cos(\pi/6) = \sqrt 3/2$), and $Q = 10, 12$ quadratic ones. For
$Q = 2^k Q_0$ with $Q_0 = 1, 3, 5$, the half-angle formulas
$$
\cos\frac\theta2 = \sqrt{\frac{1 + \cos\theta}{2}} \quad
\Bigl(\theta \le \frac\pi2\Bigr), \qquad
\cos\frac\theta2 = \sqrt{\frac{1 - \cos(\pi - \theta)}{2}} \quad
\Bigl(\theta > \frac\pi2\Bigr)
$$
lead from a closed form at the denominator $2^j Q_0 \le 6$ to
$\cos(\pi R/Q)$ in $k$ or fewer square roots. The second form cancels:
$1 - \cos(\pi s/(Q/2)) \ge 1 - \cos(2\pi/Q)$, and the absolute error of
the step grows by at most a factor $Q/4$, so each step costs at most
$\operatorname{bits}(Q)$ guard bits. Of $\sin(\pi R/Q)$ and
$\cos(\pi R/Q)$, one value alone is its own chain ($\sin(\pi R/Q) =
\cos(\pi(Q - 2R)/(2Q))$, one step more); for both, the larger is the
chain and the other $\sqrt{1 - x^2}$, which costs $2\operatorname{bits}(Q)$
more guard bits for the cancellation. With $L = k + 1$, the chain beats
the general kernel from about $L^2$ limbs ($Q = 8$ from 16 limbs,
$Q = 16$ from 24, $Q = 1024$ from 128), and for $L \le 6$ it stays
10-20% faster than the Chebyshev step at any size; for larger $L$ it is
used up to 1024 limbs, and beyond for the cosine alone (5-25% faster
than the Chebyshev step there, which wins for the sine and for both).

## Choosing the method {#choice}

`_mp_real_sin_cos_pi_ui_div_ui` takes, in this order, the exact values
($Q \le 2$), the closed forms ($Q = 3, 4, 6$; $Q = 5, 10, 12$ from 12
limbs), the half-angle steps for $Q = 2^k Q_0$ where they win, the
Chebyshev step for denominators up to the tuned bit length, and the
kernel on $\pi R/Q$ (or its complement) otherwise; the tangent uses the
closed forms $\sqrt 3$, $\sqrt 2 \mp 1$, $2 \mp \sqrt 3$, and the
sine and cosine of the half-angle steps or the Chebyshev step divided,
where they win, and the tangent kernel otherwise.
`_mp_real_sin_cos_pi_ui_div_ui_algebraic` tells whether the choice is
other than the kernel.

The callers use this predicate:

- `_mp_real_sin_cos_pi_fp` and `_mp_real_tan_pi_fp`, the fixed-point
  kernels behind `arb_sin_pi`, `arf_sin_pi`, `nfloat_sin_pi` and their
  relatives, reduce $x$ exactly to $t \in [0, 1/4]$; when $t = m/2^j$
  with $j \le \bits - 8$ and the predicate holds, $\sin(\pi t)$ and
  $\cos(\pi t)$ (or $\tan(\pi t)$) come from the algebraic methods:
  $\sin(\pi/4)$ is one square root at any precision, $\sin(3\pi/16)$
  half-angle steps from about 24 limbs.
- `arb_sin_cos_pi_fmpq`, `arb_sin_pi_fmpq`, `arb_cos_pi_fmpq` and
  `arb_tan_pi_fmpq` reduce $x$ exactly to $v/w \in [0, 1/4]$ by
  symmetries. Below 12288 bits, square roots through `arb_sqrt` are
  faster than those of mp_real: the closed forms for $w \le 6$, and the
  half-angle steps for $w = 2^k Q_0$ from a precision growing with the
  number $r = k - 1, k, k + 1$ ($Q_0 = 1, 3, 5$) of square roots (768
  bits for $r \le 3$ to 3072 bits for $r = 8$ for both values, 256 to
  6144 bits up to $r = 10$ for one of them; not for the tangent at $w =
  8, 12$, whose closed forms in mp_real are faster). Otherwise a
  word-size denominator goes to `mp_real_sin_cos_pi_ui_div_ui` resp.
  `mp_real_tan_pi_ui_div_ui` where the predicate holds, and the
  remaining cases evaluate the functions of an enclosure of $\pi v/w$,
  with full relative accuracy since the reduced argument is small or
  away from the zeros; the tangent is the quotient of the sine and
  cosine (or its inverse).

## Use in FLINT {#use}

::: functions mp_real_sin_cos_pi_ui_div_ui _mp_real_sin_cos_pi_ui_div_ui arb_sin_cos_pi_fmpq arb_tan_pi_fmpq :::
