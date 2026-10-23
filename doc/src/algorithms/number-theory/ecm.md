---
title: The elliptic curve method
summary: Lenstra's elliptic curve method for finding a prime factor $p$ of $N$, with Montgomery curves, Suyama's parametrisation, the Montgomery ladder in stage 1 and the standard continuation in stage 2.
complexity: Heuristically $L_p[1/2, \sqrt 2]$ multiplications modulo $N$ to find a prime factor $p$.
order: 31
---

The elliptic curve method (ECM) of Lenstra [@Len1987] finds a prime
factor $p$ of $N$ in a time that depends mainly on the size of $p$. It is
the method of choice for factors of 20 to 40 digits in numbers too large
for the [quadratic sieve](algorithms/number-theory/quadratic-sieve), and
`fmpz_factor` runs it before the sieve to remove medium-sized factors.
The implementation follows Montgomery [@Mon1987] and Zimmermann and
Dodson [@ZD2006].

## The idea {#idea}

Pollard's $p - 1$ method computes $a^E \bmod N$ for a highly composite
exponent $E$, and finds $p$ if the order of $a$ modulo $p$, a divisor of
$p - 1$, divides $E$. It thus succeeds only if $p - 1$ is smooth. ECM
replaces the multiplicative group of $\mathbb F_p$ by the group of points
of an elliptic curve over $\mathbb F_p$. By Hasse's theorem
[@Sil2009, Theorem V.1.1] its order lies in
$[p + 1 - 2\sqrt p, p + 1 + 2\sqrt p]$, and as the curve varies the order
is distributed over this interval roughly like a random integer of that
size. Trying many curves gives many independent chances of a smooth group
order.

One computes with a curve and a point $P$ defined modulo $N$, as if
$\Z/N\Z$ were a field. In projective coordinates, a multiple $kP$ is a
point $(X : Z)$; if the order of $P$ modulo $p$ divides $k$, then $kP$ is
the point at infinity modulo $p$, so $p$ divides $Z$ and $\gcd(Z, N)$
reveals $p$, unless the same happens modulo all prime factors of $N$ at
once, which is unlikely.

## Montgomery curves {#montgomery}

FLINT uses curves in Montgomery form $B y^2 = x^3 + A x^2 + x$ and
represents a point by its $x$-coordinate in projective form $(X : Z)$,
which identifies $P$ and $-P$ but suffices for scalar multiplication.
With $a_{24} = (A + 2)/4$, the doubling and the *differential* addition
(which needs the difference of the two points) are

$$
\begin{aligned}
X_{2P} &= (X + Z)^2 (X - Z)^2, \\
Z_{2P} &= 4XZ \bigl((X - Z)^2 + a_{24} \cdot 4XZ\bigr), \\
X_{P+Q} &= Z_{P-Q} \bigl((X_P - Z_P)(X_Q + Z_Q) + (X_P + Z_P)(X_Q - Z_Q)\bigr)^2, \\
Z_{P+Q} &= X_{P-Q} \bigl((X_P - Z_P)(X_Q + Z_Q) - (X_P + Z_P)(X_Q - Z_Q)\bigr)^2,
\end{aligned}
$$

about ten multiplications modulo $N$ for a doubling and an addition
together, and no inversion. The *Montgomery ladder* computes $kP$ by
keeping the pair $(mP, (m+1)P)$, whose difference is always $P$, and
replacing it by $(2mP, (2m+1)P)$ or $((2m+1)P, (2m+2)P)$ for each bit of
$k$.

Curves are chosen by Suyama's parametrisation [@Suy1985; @ZD2006]: for a
random $\sigma$, put $u = \sigma^2 - 5$, $v = 4\sigma$, the starting point
$x_0 = u^3/v^3$ and $a_{24} = (v - u)^3(3u + v)/(16 u^3 v)$. The group
order of these curves is divisible by 12, which increases the probability
that it is smooth. The inversion needed for $a_{24}$ may itself reveal a
factor of $N$.

## Stage 1 {#stage1}

::: algorithm ECM stage 1 {#alg-stage1}
**Input:** $N$, a curve and a point $P$ modulo $N$, a bound $B_1$.

**Output:** a factor of $N$, or the point $Q = kP$ with
$k = \prod_{p \le B_1} p^{\lfloor \log_p B_1 \rfloor}$.

1. For each prime $p \le B_1$, replace $P$ by $p^{\lfloor \log_p B_1
   \rfloor} P$ with the Montgomery ladder, and if
   $1 < \gcd(Z_P, N) < N$, return the GCD.
2. Return $Q = P$.
:::

If the order of the starting point modulo $p$ is $B_1$-smooth (a product
of prime powers at most $B_1$), it divides $k$ and stage 1 finds $p$.
The cost is about $1.44 B_1$ ladder steps. FLINT takes the GCD after
every prime, which costs little since $B_1$ is small in its applications.

## Stage 2 {#stage2}

Often the order of $P$ is $B_1$-smooth except for one larger prime $\ell$
with $B_1 < \ell \le B_2$. Stage 2 tests all these primes at once with
the *standard continuation*. Choose a primorial $D$ (the product of the
first few primes) close to $\sqrt{B_2}$; every prime $\ell$ in
$(B_1, B_2]$ can be written $\ell = mD \pm j$ with $1 \le j < D/2$ and
$\gcd(j, D) = 1$.

::: lemma {#lem-stage2}
Let $Q$ be a point modulo a prime $p$ and $\ell = mD \pm j$. If
$\ell Q = O$, then $x(mDQ) = x(jQ)$, that is,
$X_{mDQ} Z_{jQ} - X_{jQ} Z_{mDQ} \equiv 0 \pmod p$.
:::

::: proof
$\ell Q = O$ means $mDQ = \mp jQ$, and $x(-R) = x(R)$ for every point
$R$, since the $x$-coordinate does not distinguish a point from its
negative.
:::

::: algorithm ECM stage 2 (standard continuation) {#alg-stage2}
**Input:** $N$, the point $Q$ from stage 1, bounds $B_1 < B_2$, the
primorial $D$.

**Output:** a factor of $N$ or failure.

1. Compute $jQ$ for all odd $j < D/2$ (by differential additions of
   $2Q$), and $DQ$.
2. $g \gets 1$. For $m$ from about $B_1/D$ to about $B_2/D$ (computing
   $mDQ$ by differential additions of $DQ$): for each $j$ with
   $\gcd(j, D) = 1$ such that $mD + j$ or $mD - j$ is prime, set
   $g \gets g \cdot (X_{mDQ} Z_{jQ} - X_{jQ} Z_{mDQ}) \bmod N$.
3. Return $\gcd(g, N)$ if it is a nontrivial factor.
:::

By [](#lem-stage2), if the order of $Q$ modulo $p$ is a prime in
$(B_1, B_2]$, then $p \mid g$. The cost is about $D/4 + B_2/D$ point
operations and three multiplications per prime $\ell$, and a single GCD.
Asymptotically faster versions of stage 2, using fast polynomial
arithmetic [@ZD2006], are not implemented.

## Parameters and complexity {#complexity}

The probability that a random integer near $p$ is $L_p[1/2, \alpha]$-smooth
is heuristically $L_p[1/2, -1/(2\alpha)]$. Taking
$B_1 = L_p[1/2, 1/\sqrt 2]$ balances the cost of a curve against the
number of curves needed, for a total of $L_p[1/2, \sqrt 2]$ operations
modulo $N$ [@Len1987]. In practice the parameters are tuned
experimentally: `fmpz_factor_smooth` uses $B_2 = 100 B_1$ and the values
of $B_1$ and of the number of curves shown on
[](algorithms/number-theory/integer-factoring), which find a factor of
the given size with probability about $1/3$ per size class. For example,
to find a factor of 40 bits it runs 18 curves with $B_1 = 610$, and for
60 bits 42 curves with $B_1 = 14000$.

## Use in FLINT {#use}

- `fmpz_factor_ecm` for multiprecision $N$, with the arithmetic modulo $N$
  done by `flint_mpn_mulmod_preinvn` on normalised limb arrays. The curves
  are distributed over the available threads; the factor returned is the
  one found by the curve of lowest index, so the result does not depend
  on the number of threads.
- `n_factor_ecm` for one-word and `n_ll_factor_ecm` for two-word $N$,
  with the same curves and stages in word arithmetic.
- `fmpz_factor_smooth` and `fmpz_factor`, which call `fmpz_factor_ecm`
  with the tuned parameters.

The stages are exposed as `fmpz_factor_ecm_select_curve`,
`fmpz_factor_ecm_stage_I` and `fmpz_factor_ecm_stage_II`.

::: functions fmpz_factor_ecm n_factor_ecm n_ll_factor_ecm fmpz_factor_smooth :::
