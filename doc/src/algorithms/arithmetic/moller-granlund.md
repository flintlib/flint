---
title: Division by an invariant word
summary: Divides a two-word number by a fixed word using one multiplication by a precomputed reciprocal, one low multiplication and two cheap corrections, instead of a hardware division.
complexity: One $1 \times 1 \to 2$ word product, one low product, a two-word addition and two comparisons. The reciprocal costs one division.
order: 10
---

To divide a two-word number by a fixed single word without a hardware
division, FLINT uses the algorithm of Möller and Granlund [@GraMol2010]
(Algorithm 4 in their paper), which multiplies by a precomputed
reciprocal of the divisor. Every `_preinv` function in `ulong_extras`
that reduces or divides by a single word uses it.

## The reciprocal {#reciprocal}

Let $B = 2^{\bits}$ and let $d$ be *normalised*, meaning $B/2 \le d < B$.
Its reciprocal is

$$
v = \left\lfloor \frac{B^2 - 1}{d} \right\rfloor - B.
$$

Since $B < (B^2 - 1)/d < 2B$ for normalised $d$, we have $0 < v < B$ and $v$
fits in a word. Equivalently, $v$ is the quotient of the two-word number
$(B - 1 - d) B + (B - 1)$ by $d$, which is how `n_preinvert_limb_prenorm`
computes it with a single hardware division. On platforms without a fast
two-by-one division instruction, a table lookup followed by Newton steps is
used instead, as in GMP.

An arbitrary divisor $d \ne 0$ is handled by shifting: with
$s$ the number of leading zero bits of $d$, divide $u 2^s$ by $d 2^s$ and
shift the remainder right by $s$ bits. The quotient is unchanged. This is
why `n_preinvert_limb` accepts any nonzero $n$: it returns the reciprocal of
$n 2^s$.

## The algorithm {#algorithm}

::: algorithm Two-by-one division with a precomputed reciprocal {#alg-div}
**Input:** normalised $d$, its reciprocal $v$, and $u = u_1 B + u_0$ with
$0 \le u_1 < d$ and $0 \le u_0 < B$.

**Output:** $q, r$ with $u = qd + r$ and $0 \le r < d$.

1. $(q_1, q_0) \gets v u_1 + (u_1 B + u_0)$, a two-word result.
2. $q_1 \gets (q_1 + 1) \bmod B$.
3. $r \gets (u_0 - q_1 d) \bmod B$.
4. If $r > q_0$: set $q_1 \gets (q_1 - 1) \bmod B$ and $r \gets (r + d) \bmod B$.
5. If $r \ge d$: set $q_1 \gets q_1 + 1$ and $r \gets r - d$.
6. Return $q = q_1$ and $r$.
:::

Step 4 can be implemented without branches using a mask, and step 5 is
rarely taken. Compared with the older method of Granlund and Montgomery
[@GraMon1994], which needs a high product and a full product in the
correction, this algorithm replaces one full product with a cheaper low
product.

## Correctness {#correctness}

Let $m = B + v = \lfloor (B^2 - 1)/d \rfloor$ and define $k$ by
$md = B^2 - k$. Then $1 \le k \le d$, because
$B^2 - 1 - d < md \le B^2 - 1$.

::: lemma {#lem-fits}
Step 1 computes $q_1 B + q_0 = m u_1 + u_0$ exactly, with $q_1 < B$.
:::

::: proof
Indeed $v u_1 + u_1 B + u_0 = m u_1 + u_0$. Since $md \le B^2 - 1$ and
$u_1 \le d - 1$, we get $m u_1 + u_0 \le B^2 - 1 - m + B - 1 < B^2$ because
$m > B$.
:::

Let $q_1, q_0$ denote the values after step 1, and let
$R = u - (q_1 + 1) d$ be the exact remainder belonging to the candidate
quotient $q_1 + 1$ of step 2.

::: lemma Bounds on the candidate remainder {#lem-bounds}
With $M = \max(q_0, B - d)$ we have
$q_0 - B < R < M$ and $R \ge -d$.
:::

::: proof
Multiply the identity of [](#lem-fits) by $d$ and use $md = B^2 - k$:
$d q_1 B = (B^2 - k) u_1 + d u_0 - d q_0$. Substituting into
$B R = B(u_1 B + u_0) - B d q_1 - B d$ gives

$$
B R = k u_1 + (B - d) u_0 + d q_0 - B d.
$$

For the lower bounds, the first two terms are nonnegative, so
$R \ge d (q_0 - B)/B$. As $q_0 - B < 0$ and $d < B$, this is greater than
$q_0 - B$; as $q_0 \ge 0$, it is at least $-d$.

For the upper bound, use $k \le d$, $u_1 \le d - 1$ and $u_0 \le B - 1$:

$$
\begin{aligned}
B R &\le d(d-1) + (B - d)(B - 1) + d q_0 - B d \\
&= (B - d)^2 - B + d q_0 < (B - d) M + d M = B M.
\end{aligned}
$$
:::

::: theorem {#thm-correct}
[](#alg-div) returns the correct quotient and remainder.
:::

::: proof
Step 3 computes $r = R \bmod B$, since $u - (q_1 + 1) d \equiv u_0 - (q_1 + 1) d$
modulo $B$. Since $d$ is normalised, $B \le 2d$. There are three cases.

- $R < 0$. Then $r = R + B > q_0$ by [](#lem-bounds), so step 4 decrements
  the quotient and sets $r = R + d$, which lies in $[0, d)$ because
  $-d \le R < 0$. Step 5 does nothing.
- $R \ge 0$ and $r = R \le q_0$. Step 4 does nothing, and
  $0 \le R < M < B \le 2d$, so step 5 performs any necessary correction.
- $R \ge 0$ and $r = R > q_0$. This can only happen when $M = B - d$, so
  $R < B - d$. Step 4 wrongly decrements the quotient and sets
  $r = R + d < B$, but then $d \le r < 2d$, and step 5 undoes the change.

In all cases the final remainder is in $[0, d)$. Every update of $q_1$ is
an addition of $\pm 1$ modulo $B$, and the true quotient lies in $[0, B)$
because $u_1 < d$, so computing the quotient modulo $B$ gives the right
answer even if step 2 wraps around.
:::

## Use in FLINT {#use}

The reciprocal of [](#reciprocal) is the `ninv` field of every `nmod_t`, so
this algorithm is the default modular reduction for all of FLINT's
arithmetic modulo word-size integers.

- `nmod_mul` multiplies and reduces the two-word product with one step;
  the macros `NMOD_RED2` and `NMOD_RED3` reduce two and three words with
  one and two steps. `NMOD_RED`, for one word, uses instead a Barrett step
  whose constant is derived from the same reciprocal
  ([](algorithms/arithmetic/barrett#nmod)). They are used throughout
  `nmod_vec`, `nmod_poly`, `nmod_mat`, `nmod_mpoly` and the `gr` contexts
  for $\Z/n\Z$; in particular every dot product in `nmod_vec` ends with one
  such reduction ([](algorithms/linear-algebra/dot-products)).
- `n_mod2_preinv`, `n_divrem2_preinv` and `n_div2_preinv` shift a one-word
  dividend into two words and run [](#alg-div) once. The condition
  $u_1 < d$ holds automatically after shifting.
- `n_ll_mod_preinv` first reduces the high word modulo $n$ if necessary and
  then runs the algorithm on the two-word dividend.
- `n_lll_mod_preinv` reduces the top two words and then the bottom two.
- `n_mulmod2_preinv` multiplies and reduces the two-word product with
  `n_ll_mod_preinv`; `n_mulmod_preinv` skips the normalisation shift by
  requiring a normalised modulus.
