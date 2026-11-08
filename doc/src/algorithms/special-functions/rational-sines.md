---
title: Sines and cosines of rational multiples of π
summary: Evaluating cos(πR/Q) and sin(πR/Q) as algebraic numbers - a Chebyshev step of arbitrary order with hypergeometric corrections that needs no division or square root (mp_real), and Newton refinement of a root of the minimal polynomial (arb).
complexity: $O(M(n)\log Q)$ for the Chebyshev iteration; Newton refinement costs $O(M(n))$ times the degree of the minimal polynomial.
order: 15
---

The values $\cos(\pi R/Q)$ are algebraic of degree at most $Q$, and at
high precision it is cheaper to refine one as a root of a polynomial
than to evaluate the cosine with the general kernels of
[](algorithms/special-functions/elementary). FLINT does this in two
ways: `mp_real_sin_cos_pi_ui_div_ui` refines an approximation by a
Chebyshev step whose order grows with the precision, and
`arb_sin_cos_pi_fmpq` refines a root of the minimal polynomial by
Newton iteration for small denominators.

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

## Minimal polynomials {#minpoly}

`arb_sin_cos_pi_fmpq` reduces $x$ to $[0, 1/4]$ by symmetries and then
evaluates $\cos(\pi R/Q)$ either as $\cos$ of an enclosure of
$\pi R/Q$ or, when the denominator is small against the precision
(up to $Q = 1000$, from about $1500 + 150 Q$ bits), as a root of the
minimal polynomial of $2\cos(2\pi/Q)$ (`fmpz_poly_cos_minpoly`): a
starting value at about 100 bits plus the size of the coefficients is
refined by Newton iteration in ball arithmetic
(`_arb_poly_newton_refine_root`), whose convergence on an interval
around the start is checked first. Denominators up to 6 use closed
forms, and an even denominator is halved by the half-angle formula
$\cos(\theta/2) = \sqrt{(1 + \cos\theta)/2}$.

## Use in FLINT {#use}

::: functions mp_real_sin_cos_pi_ui_div_ui _mp_real_sin_cos_pi_ui_div_ui arb_sin_cos_pi_fmpq :::
