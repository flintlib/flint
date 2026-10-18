---
title: Multiplication by a fixed operand (Shoup)
summary: Computes $ab \bmod n$ for fixed $a$ and $n$ with two low products, one high product and one comparison, using the precomputed quotient $\lfloor aB/n \rfloor$.
complexity: One high product, two low products and one conditional subtraction per multiplication; one division per fixed operand.
order: 20
---

When the same $a$ and $n$ are used for many products $ab \bmod n$, as in
scaling a vector or applying a fixed twiddle factor in a number-theoretic
transform, a quotient can be precomputed for $a$ so that each product needs
no division at all. The method was introduced in NTL by Shoup
[@Sho2021] and is analysed by Harvey [@Har2014].

Throughout, $B = 2^{\bits}$ and the modulus satisfies $n < B/2$. This
restriction is what makes the final correction a single comparison.

## Multiplication {#mul}

The precomputed data for $0 \le a < n$ is $\hat a = \lfloor aB/n \rfloor$,
which fits in a word because $a < n$. We write $aB = \hat a n + \check a$
with $0 \le \check a < n$.

::: algorithm Shoup multiplication {#alg-mul}
**Input:** $a < n < B/2$, $\hat a = \lfloor aB/n \rfloor$ and any
$0 \le b < B$.

**Output:** $ab \bmod n$.

1. $p \gets \lfloor \hat a b / B \rfloor$ (the high word of $\hat a b$).
2. $c \gets (ab - pn) \bmod B$ (two low products).
3. If $c \ge n$, set $c \gets c - n$. Return $c$.
:::

::: lemma {#lem-quot}
$p = \lfloor ab/n \rfloor$ or $p = \lfloor ab/n \rfloor - 1$.
:::

::: proof
From $aB = \hat a n + \check a$ we get
$\hat a b / B = ab/n - \check a b/(nB)$. Since $0 \le \check a < n$ and
$0 \le b < B$, the subtracted term lies in $[0, 1)$, so
$ab/n - 1 < \hat a b / B \le ab/n$. Taking floors gives the claim.
:::

::: theorem {#thm-mul}
[](#alg-mul) returns $ab \bmod n$.
:::

::: proof
By [](#lem-quot), the exact value $ab - pn$ equals either $ab \bmod n$ or
$(ab \bmod n) + n$. Both are less than $2n < B$, so the computation modulo
$B$ in step 2 is exact, and the comparison in step 3 tells the two cases
apart. If $n$ were allowed to be as large as $B$, then $c$ and $c + n$
could not be distinguished after reduction modulo $B$.
:::

## Updating the precomputation {#update}

To multiply by $ab \bmod n$ repeatedly, for example when building a table
of powers $a, a^2, a^3, \ldots$, we also need its precomputed quotient.
`n_mulmod_and_precomp_shoup` computes it without a division. Now let also
$b < n$ and $bB = \hat b n + \check b$.

::: algorithm Product with its precomputed quotient {#alg-update}
**Input:** $a, b < n < B/2$, $\hat a$, $\check a$ and $\hat b$.

**Output:** $c = ab \bmod n$ and $\hat c = \lfloor cB/n \rfloor$.

1. Compute the two-word product $\hat a b = p B + e$.
2. $c \gets (ab - pn) \bmod B$; if $c \ge n$, set $c \gets c - n$.
3. $q \gets \lfloor \hat b \check a / B \rfloor$.
4. $h \gets (b \check a - qn) \bmod B$; if $h \ge n$, set $q \gets q + 1$.
5. Return $c$ and $\hat c = (e + q) \bmod B$.
:::

::: theorem {#thm-update}
[](#alg-update) is correct.
:::

::: proof
Steps 1 and 2 are [](#alg-mul), so $c = ab \bmod n$. For the quotient,

$$
\left\lfloor \frac{abB}{n} \right\rfloor
= \left\lfloor \frac{(\hat a n + \check a) b}{n} \right\rfloor
= \hat a b + \left\lfloor \frac{\check a b}{n} \right\rfloor.
$$

Since $c = ab - \lfloor ab/n \rfloor n$, we have
$\lfloor cB/n \rfloor = \lfloor abB/n \rfloor - \lfloor ab/n \rfloor B$,
and so $\lfloor cB/n \rfloor \equiv e + \lfloor \check a b/n \rfloor \pmod B$.
As $c < n$, the left side is less than $B$, so it equals the right side
reduced modulo $B$.

It remains to see that steps 3 and 4 compute $\lfloor \check a b / n \rfloor$.
This is [](#lem-quot) with the roles of $a$ and $b$ exchanged (valid
because $b < n$): the value $q$ of step 3 is the correct quotient or one
less, and $\check a b - qn < 2n$ detects which.
:::

## Use in FLINT {#use}

The macro `NMOD_CAN_USE_SHOUP(mod)` tests the condition $n < 2^{\bits - 1}$
for an `nmod_t`. Besides the functions in `ulong_extras`, the method is
used by:

- `_nmod_vec_scalar_mul_nmod_shoup` and `_nmod_vec_scalar_addmul_nmod_shoup`,
  which `_nmod_vec_scalar_mul_nmod` and `_nmod_vec_scalar_addmul_nmod`
  select above a tuned length (see [](nmod_vec/scalar));
- `_nmod_vec_invert_shoup`, which builds the prefix products of
  [](algorithms/arithmetic/batch-inversion) together with their
  precomputed quotients using [](#alg-update);
- scalar multiplication in `nmod_mat`, evaluation of polynomials at a
  point and at geometric progressions in `nmod_poly`, and modular
  multiplication by a fixed multiplier and Chinese remaindering in
  `mpn_extras`.

On SIMD hardware the same multiply-by-a-constant pattern is usually
handled in floating point instead
([](algorithms/simd/fp-modular-arithmetic)).

::: functions n_mulmod_shoup _nmod_vec_scalar_mul_nmod_shoup _nmod_vec_invert_shoup :::
