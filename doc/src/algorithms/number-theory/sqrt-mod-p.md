---
title: Square roots modulo a prime
summary: How `nmod_sqrt` finds a square root modulo a word-size prime $n$. It uses one exponentiation when $n \equiv 3 \pmod 4$ or $n \equiv 5 \pmod 8$, and otherwise Tonelli–Shanks in two forms or Müller's variant of Cipolla–Lehmer, depending on the power of 2 dividing $n - 1$. Every candidate root is checked by squaring.
complexity: $O(\log n)$ modular multiplications if $n \not\equiv 1 \pmod 8$. Otherwise, with $n - 1 = q 2^s$ and $q$ odd, $O(\log n + s^2)$ multiplications for Tonelli–Shanks, or $O(\log n)$ multiplications and a few Jacobi symbols for Müller's method.
order: 30
---

Given a prime $n$ and $a \in \Z/n\Z$, we want $r$ with $r^2 = a$. Half of
the nonzero residues are squares. For $a$ a square, the two roots are $r$
and $n - r$, and `nmod_sqrt` returns the one in $[0, n/2]$. The choice of
method depends on $n \bmod 8$ and on the 2-adic valuation $s$ of $n - 1$.
Standard references are [@Coh1996] and [@CraPom2005].

## Overview {#overview}

`nmod_sqrt(res, a, mod)` proceeds as follows.

1. If $a \le 1$, return $a$. If $n$ is even, return failure.
2. If $n < 128$ (`NMOD_SQRT_SMALL_CUTOFF`), search $t = 1, 2, \ldots,
   \lfloor n/2 \rfloor$ with $t^2$ updated incrementally
   (`_nmod_sqrt_small`). This is exact for any odd $n$.
3. Otherwise choose the arithmetic: half-word lazy
   [Montgomery arithmetic](algorithms/arithmetic/montgomery) if
   $n < 2^{\bits/2 - 2}$, full-word lazy Montgomery arithmetic if
   $n < 2^{\bits - 2}$, and plain Montgomery arithmetic otherwise. All
   remaining steps work in Montgomery form, and the three versions are
   compiled separately so that each is fully inlined.
4. Compute a candidate root $r$:
   - $n \equiv 3 \pmod 4$: [](#alg-3mod4);
   - $n \equiv 5 \pmod 8$: [](#alg-atkin);
   - $n \equiv 1 \pmod 8$: [](#alg-mueller) if $6s \ge b + 120$, where $b$
     is the bit length of $n$, and if it finds a suitable parameter;
     otherwise [](#alg-ts) for $s < 12$ and [](#alg-ts-bits) for $s \ge 12$.
5. Return $r$ if $r^2 = a$, and failure otherwise.

No Jacobi symbol of $a$ is computed in advance. For prime $n$, every
method below returns a correct root when $a$ is a square, so the final
check is a complete test of whether $a$ is a square. For composite $n$, the
check guarantees that a returned root is genuine, but a failure proves
nothing.

The branches in the inner loops are replaced by masks and selections,
since their conditions depend on random data and would be mispredicted.
For the same reason, exponentiation multiplies by 1 instead of skipping a
multiplication when an exponent bit is 0.

## The case $n \equiv 3 \pmod 4$ {#mod4}

::: algorithm {#alg-3mod4}
**Input:** prime $n \equiv 3 \pmod 4$ and $a$.

**Output:** a candidate root $r$.

1. Return $r = a^{(n+1)/4}$.
:::

If $a$ is a square, Euler's criterion gives $a^{(n-1)/2} = 1$, so
$r^2 = a^{(n+1)/2} = a \cdot a^{(n-1)/2} = a$.

## The case $n \equiv 5 \pmod 8$ {#mod8}

Here $-1$ is a square and $2$ is not. Atkin's method [@Atk1992] uses a
square root of $-1$ built from $2a$.

::: algorithm Atkin {#alg-atkin}
**Input:** prime $n \equiv 5 \pmod 8$ and $a \ne 0$.

**Output:** a candidate root $r$.

1. $t \gets 2a$ and $b \gets t^{(n-5)/8}$.
2. $i \gets t b^2$.
3. Return $r = ab(i - 1)$.
:::

::: theorem {#thm-atkin}
If $a$ is a nonzero square modulo the prime $n \equiv 5 \pmod 8$, then
[](#alg-atkin) returns a square root of $a$.
:::

::: proof
We have $i = t^{(n-1)/4}$, so $i^2 = (2a)^{(n-1)/2} = -1$, because 2 is a
nonsquare and $a$ is a square. Hence $(i - 1)^2 = i^2 - 2i + 1 = -2i$ and

$$
r^2 = a^2 b^2 (-2i) = -a i \cdot (2a b^2) = -a i^2 = a.
$$
:::

## Tonelli–Shanks {#tonelli-shanks}

For $n \equiv 1 \pmod 8$ write $n - 1 = q 2^s$ with $q$ odd, so $s \ge 3$.
The method of Tonelli [@Ton1891] and Shanks [@Sha1973] works in the
subgroup of order $2^s$ of $(\Z/n\Z)^\times$, which is cyclic.

**Nonresidue.** It needs some $k$ that is not a square. FLINT takes the
least odd prime $k \le 61$ with Jacobi symbol $(k \mid n) = -1$, and
otherwise the least odd $k \ge 67$ with $(k \mid n) = -1$. Only odd primes
need to be tried among the small values, because the least nonresidue is
prime and 2 is a square for $n \equiv 1 \pmod 8$. By quadratic
reciprocity $(k \mid n) = (n \mid k)$, and $(n \mid k)$ is read from a
bitmask indexed by $n \bmod k$, computed with two multiplications instead
of a division: with $m = n k^{-1} \bmod B$, the high word $h$ of $mk$
satisfies $mk = n + hB$, so $n \equiv -hB \pmod k$, and since $B$ is a
square, $(n \mid k) = (-h \mid k)$. The search gives up, which makes
`nmod_sqrt` fail, if $n$ is found to be a perfect square or no nonresidue
is found below $2^{20}$; neither happens for prime $n$.

**Setup.** With $u = a^{(q-1)/2}$ and $c = k^{(q-1)/2}$, both computed by one
exponentiation loop that interleaves the two independent squaring chains,
set

$$
r = ua = a^{(q+1)/2}, \qquad t = ur = a^q, \qquad g = c^2 k = k^q.
$$

Then $r^2 = a t$, and $g$ has order exactly $2^s$ since
$g^{2^{s-1}} = k^{(n-1)/2} = -1$. If $a$ is a square, then
$t^{2^{s-1}} = a^{(n-1)/2} = 1$. It remains to find $x$ with $t g^x = 1$;
then $x$ is even when $a$ is a square, and $r g^{x/2}$ is a root.

::: algorithm Tonelli–Shanks descent {#alg-ts}
**Input:** $r, t, g$ as above, with $r^2 = at$ and $g$ of order $2^s$.

**Output:** a root $r$ of $a$, or failure.

1. $m \gets s$.
2. While $t \ne 1$:
   1. Find the least $i$ with $0 < i < m$ and $t^{2^i} = 1$, by repeated
      squaring. If there is none, return failure.
   2. $b \gets g^{2^{m - i - 1}}$, $g \gets b^2$, $t \gets t g$,
      $r \gets r b$, $m \gets i$.
3. Return $r$.
:::

::: theorem {#thm-ts}
If $a$ is a nonzero square modulo the prime $n$, then [](#alg-ts)
terminates and returns a square root of $a$. If $a$ is not a square, it
returns failure.
:::

::: proof
The loop maintains $r^2 = at$, $g$ of order exactly $2^m$, and
$t^{2^{m-1}} = 1$. This holds initially. In step 2.1, $t \ne 1$ and
$t^{2^{m-1}} = 1$, so $i$ exists and $t$ has order exactly $2^i$. Then $b$
has order $2^{i+1}$, so $b^2$ has order $2^i$. In a cyclic group, the
elements of order exactly $2^i$ all satisfy $x^{2^{i-1}} = -1$, so
$(t b^2)^{2^{i-1}} = 1$. The updates give
$(rb)^2 = at b^2 = a t_{\text{new}}$, and the invariant holds with $m = i$.
Since $m$ decreases strictly, the loop terminates, and on exit $t = 1$ and
$r^2 = a$. If $a$ is not a square, then $t^{2^{s-1}} = -1$, so no $i < s$
exists in the first round.
:::

Each round costs up to $m$ dependent squarings, so the total is $O(s^2)$
multiplications in one long dependency chain. This is cheap for small $s$,
which is the common case: asymptotically, half of the primes
$n \equiv 1 \pmod 8$ have $s = 3$. FLINT uses it for $s < 12$.

### Bitwise variant {#ts-bits}

For larger $s$ the latency of the dependent squarings dominates. The
second form finds the bits $x_0, x_1, \ldots$ of $x$ from the bottom, as in
Pohlig–Hellman (see [](algorithms/number-theory/discrete-log)), using
products that do not depend on each other.

::: algorithm Tonelli–Shanks by bits {#alg-ts-bits}
**Input:** $r, t, g$ as above.

**Output:** a root $r$ of $a$, or failure.

1. For $0 \le i < s$: $D_i \gets t^{2^i}$ and $G_i \gets g^{2^i}$.
2. If $D_{s-1} \ne 1$, return failure.
3. For $j = 1, \ldots, s - 1$: if $D_{s-1-j} \ne 1$, set
   $D_i \gets D_i G_{i+j}$ for $0 \le i < s - 1 - j$, and $r \gets r G_{j-1}$.
4. Return $r$.
:::

::: theorem {#thm-ts-bits}
If $a$ is a nonzero square modulo the prime $n$, then [](#alg-ts-bits)
returns a square root of $a$. If $a$ is not a square, it returns failure.
:::

::: proof
Let $x \in [0, 2^s)$ with $t g^x = 1$, which exists because $t$ lies in the
cyclic group generated by $g$. Let $X_j$ be the part of $x$ below bit $j$.
We claim that before iteration $j$ of step 3, $D_i = (t g^{X_j})^{2^i}$ for
$i \le s - 1 - j$, and $r = a^{(q+1)/2} g^{X_j/2}$.

Since $t g^{X_j} = g^{-(x - X_j)}$ and $2^j$ divides $x - X_j$, we get
$D_{s-1-j} = g^{-2^{s-1} (x - X_j)/2^j} = (-1)^{x_j}$, so the test in
step 3 reads off bit $x_j$. Setting the bit multiplies $t g^{X_j}$ by
$g^{2^j}$, hence each $D_i$ by $G_{i+j}$, and $r$ by $g^{2^{j-1}} = G_{j-1}$.
Only $D_i$ with $i < s - 1 - j$ are used later. Step 2 is the case $j = 0$:
$D_{s-1} = a^{(n-1)/2}$, which is 1 if and only if $a$ is a square, and
then $x_0 = 0$. At the end $X = x$ and
$r^2 = a^{q+1} g^{x} = a \cdot t g^x = a$.
:::

Step 1 is two independent chains of $s - 1$ squarings, and step 3 performs
at most $s(s-1)/2$ multiplications. The work is still $O(s^2)$, but the
multiplications of each iteration are independent, so their cost is bounded
by multiplication throughput instead of latency.

## Müller's method {#mueller}

The cost of Tonelli–Shanks grows with $s$. The Cipolla–Lehmer method
works in $\mathbb{F}_{n^2}$ instead, with a cost that does not depend on
$s$. Müller [@Mue2004] gave a version for $n \equiv 1 \pmod 4$ that needs
a single Lucas sequence. For $P \in \Z/n\Z$, let $V_0 = 2$, $V_1 = P$ and
$V_{k+1} = P V_k - V_{k-1}$, so that $V_k = \alpha^k + \alpha^{-k}$ where
$\alpha, \alpha^{-1}$ are the roots of $X^2 - PX + 1$.

::: theorem Müller {#thm-mueller}
Let $n \equiv 1 \pmod 4$ be prime, let $a$ be a nonzero square, and let
$t \ne 0$ with $(a t^2 - 4 \mid n) = -1$. With $P = a t^2 - 2$,

$$
\left(\frac{V_{(n-1)/4}(P)}{t}\right)^2 = a.
$$
:::

::: proof
Let $u = t \sqrt a \in \Z/n\Z$, and let $w$ be a square root of
$u^2 - 4 = a t^2 - 4$ in $\mathbb{F}_{n^2}$; it is not in $\Z/n\Z$. Put
$\beta = (u + w)/2$ and $\beta' = (u - w)/2$. Then $\beta \beta' = 1$,
$\beta + \beta' = u$ and $\beta^2 + \beta'^2 = u^2 - 2 = P$, so
$\alpha = \beta^2$ is a root of $X^2 - PX + 1$. The Frobenius map fixes $u$
and sends $w$ to $-w$, so $\beta^n = \beta' = \beta^{-1}$ and
$\beta^{n+1} = 1$. Hence $\varepsilon = \beta^{(n+1)/2}$ satisfies
$\varepsilon^2 = 1$ and $\varepsilon = \pm 1$. Now

$$
V_{(n-1)/4} = \beta^{(n-1)/2} + \beta^{-(n-1)/2}
= \varepsilon \beta^{-1} + \varepsilon \beta = \varepsilon u
= \pm t \sqrt a.
$$
:::

::: algorithm Müller {#alg-mueller}
**Input:** prime $n \equiv 1 \pmod 4$ and $a$.

**Output:** a candidate root $r$, or no answer.

1. For $j = 0, 1, \ldots, 32$, let $t = 2^j$ and test whether
   $(a t^2 - 4 \mid n) = -1$; stop at the first success. If none succeeds,
   give no answer.
2. $P \gets a t^2 - 2$.
3. Compute $V = V_{(n-1)/4}(P)$ with the ladder
   $(V_m, V_{m+1}) \mapsto (V_{2m}, V_{2m+1})$ or $(V_{2m+1}, V_{2m+2})$
   over the bits of $(n-1)/4$, using
   $V_{2m} = V_m^2 - 2$ and $V_{2m+1} = V_m V_{m+1} - P$.
4. Return $r = V / t$.
:::

Choosing $t = 2^j$ makes the division by $t$ in step 4 a sequence of
halvings, each a shift and a conditional addition of $(n+1)/2$. The Jacobi
symbols in step 1 are computed for $a - 4$ ($j = 0$) and for
$c_j = 4^{j-1} a - 1$ ($j \ge 1$), because $a 4^j - 4 = 4 c_j$ and 4 is a
square; they satisfy $c_1 = a - 1$ and $c_{j+1} = 4 c_j + 3$. Each
candidate fails with probability about $1/2$, so giving no answer after 33
candidates is very unlikely; `nmod_sqrt` then falls back to Tonelli–Shanks.

Each step of the ladder costs one squaring and one independent
multiplication, so the latency is about that of one exponentiation. The
fixed cost of the Jacobi symbols makes the method pay off only for large
$s$. FLINT uses it when $6s \ge b + 120$, where $b$ is the bit length of
$n$. Since $s < b$, this happens only for $n \ge 2^{25}$, and for 64-bit
$n$ it requires $s \ge 31$.

## Use in FLINT {#use}

- `nmod_sqrt` implements all the methods above. `n_sqrtmod` is a wrapper
  that skips `nmod_init` for $n < 128$.
- Other callers include `nmod_poly_sqrt` (for the constant term), the
  square root method of `gr_ctx_init_nmod` contexts, and
  `flint_mpn_sqrtmod_preinv` for one-word moduli. For multi-word moduli,
  `flint_mpn_sqrtmod_preinv` has its own implementation of Müller's
  method.

::: functions nmod_sqrt _nmod_sqrt_small n_sqrtmod :::
