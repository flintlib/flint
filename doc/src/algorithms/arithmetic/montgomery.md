---
title: Montgomery multiplication (REDC)
summary: Represents $x$ modulo an odd $n$ by $xR \bmod n$, so that a modular product needs two word multiplications and a subtraction instead of a division. Lazy variants keep residues in $[0, 2n)$ and skip the final correction.
complexity: One full product, one low product and one high or full product per multiplication, plus one conditional correction (none in the lazy variants). Conversions cost one division step with a precomputed inverse.
order: 25
---

Montgomery's method [@Mon1985] replaces the division in $ab \bmod n$ by a
division by a power of two, which is free. The price is that residues are
kept in a scaled form, and that the modulus must be odd. FLINT implements
it in `nmod.h` for full words ($R = 2^{\bits}$) and half words
($R = 2^{\bits/2}$). A textbook treatment is in [@BZ2011].

Throughout, $B = 2^{\bits}$, $n$ is odd, and $R$ is a power of two with
$n < R \le B$.

## Montgomery form {#form}

The *Montgomery form* of $x \in \Z/n\Z$ is $\tilde x = xR \bmod n$. Since
$R$ is invertible modulo $n$, this is a bijection, and

$$
\widetilde{x + y} = \tilde x + \tilde y \bmod n, \qquad
\widetilde{xy} = \tilde x \tilde y R^{-1} \bmod n.
$$

Addition and subtraction are therefore unchanged, and multiplication needs
the map $T \mapsto T R^{-1} \bmod n$, called REDC. Conversion into
Montgomery form costs one ordinary reduction of $xR$, and conversion back
is one REDC of $\tilde x$.

## Reduction {#redc}

Let $n' = -n^{-1} \bmod R$, which exists because $n$ is odd. FLINT stores
it in the field `nred` of an `nmod_redc_ctx_t` and computes it with
`n_binvert`.

::: algorithm REDC {#alg-redc}
**Input:** $0 \le T < nR$, and $n' = -n^{-1} \bmod R$.

**Output:** $T R^{-1} \bmod n$.

1. $m \gets (T \bmod R)\, n' \bmod R$.
2. $t \gets (T + mn) / R$.
3. If $t \ge n$, set $t \gets t - n$. Return $t$.
:::

::: theorem {#thm-redc}
In [](#alg-redc) the division in step 2 is exact, $t \equiv T R^{-1}
\pmod n$, and $0 \le t < T/R + n$ before step 3. In particular, the result
is correct for $T < nR$.
:::

::: proof
By the choice of $m$, $mn \equiv -T n^{-1} n = -T \pmod R$, so $R$ divides
$T + mn$. Then $tR = T + mn \equiv T \pmod n$, and multiplying by
$R^{-1}$ gives $t \equiv T R^{-1}$. Since $0 \le m \le R - 1$, we have
$0 \le tR < T + Rn$, which is the bound. If $T < nR$, then $t < 2n$, so the
single subtraction in step 3 brings $t$ into $[0, n)$.
:::

With $R = B$ and $T = ab$ for reduced $a$ and $b$, step 1 is a low product,
$mn$ is a full product, and the division by $R$ is taking the high word.

## The subtractive variant {#subtractive}

`n_ll_redc`, used by `nmod_redc_mul`, works with $n^{-1} \bmod R = -n'$
instead of $n'$. This avoids the carry from the addition in step 2, which
matters when $n$ is close to $B$ and $T + mn$ does not fit in two words.

::: algorithm REDC with subtraction {#alg-redc-sub}
**Input:** $0 \le T < nB$ as two words $T = T_1 B + T_0$, and
$v = n^{-1} \bmod B$.

**Output:** $T B^{-1} \bmod n$.

1. $m \gets T_0 v \bmod B$.
2. $h \gets \lfloor mn / B \rfloor$ (a high product).
3. If $T_1 < h$, return $T_1 - h + n$; otherwise return $T_1 - h$.
:::

::: theorem {#thm-redc-sub}
[](#alg-redc-sub) is correct for every odd $n < B$.
:::

::: proof
We have $mn \equiv T_0 v n = T_0 \pmod B$, so the low word of $mn$ equals
$T_0$ and $T - mn = (T_1 - h) B$ exactly. Hence $T_1 - h \equiv T B^{-1}
\pmod n$. Both $T_1$ and $h$ are in $[0, n)$: the first because $T < nB$,
the second because $m < B$. So $T_1 - h \in (-n, n)$, and step 3 returns
the representative in $[0, n)$.
:::

Every intermediate value fits in a word, so this works for all odd moduli
up to $B - 1$. It is what the plain `nmod_redc_` functions use.

## Lazy reduction {#lazy}

Skipping step 3 of [](#alg-redc) saves a comparison and a conditional
subtraction per multiplication, at the cost of an output in $[0, 2n)$
instead of $[0, n)$. Such *lazy* residues can be fed back into the next
multiplication if the modulus is small enough. This observation is due to
Walter [@Wal1999].

::: theorem Lazy multiplication {#thm-lazy}
Let $4n \le R$ and let $0 \le a, b < 2n$. Then REDC of $T = ab$ without
step 3 returns $t \equiv ab R^{-1} \pmod n$ with $0 \le t < 2n$.
:::

::: proof
By [](#thm-redc), $t < T/R + n < 4n^2/R + n \le n + n = 2n$.
:::

The `fast` functions use this with $R = B$, so they require
$n < 2^{\bits - 2}$; `nmod_redc_can_use_fast` tests exactly this. The
condition is needed: for $2^{\bits - 2} < n < 2^{\bits - 1}$ and $a, b$ close
to $2n$, the output of `nmod_redc_fast_mul` can exceed $2n$, and the sum
$T + mn$ can overflow two words, so that the result is not even congruent
to $ab R^{-1}$. The bound also keeps all intermediate
values in range: $T + mn < 4n^2 + Bn < B^2$, and $2n < B$ for the lazy
addition and subtraction modulo $2n$.

A lazy residue is turned into an ordinary one by subtracting $n$ once if
needed (`nmod_redc_fast_normalise`). Comparing two lazy residues for
equality requires normalising both first.

## Half-word arithmetic {#half}

With $R = 2^{\bits/2}$ and $n < R$, the product $T = ab$ of two reduced
operands fits in one word, and so does $T + mn$ provided that
$T + (R - 1) n < B$. REDC then needs only single-word multiplications:
step 1 is a product reduced modulo $R$ (a mask), step 2 a product, an
addition and a shift.

- For $a, b < n$ (the `half` functions) we need $n^2 + (R - 1)n < B$,
  which holds for $n < R/2 = 2^{\bits/2 - 1}$.
- For lazy $a, b < 2n$ (the `half_fast` functions) we need
  $4n \le R$ by [](#thm-lazy), so $n < 2^{\bits/2 - 2}$. Then also
  $4n^2 + Rn < B$.

On 64-bit machines this limits the modulus to 31 or 30 bits. Whether the
half-word version is faster depends on the processor; it saves the
two-word arithmetic but not the number of multiplications.

## Costs {#costs}

| Operation | Plain (`n_mulmod2_preinv`) | `nmod_redc_mul` | `nmod_redc_fast_mul` |
|---|---|---|---|
| Products | 2 full, 1 low | 1 full, 1 high, 1 low | 2 full, 1 low |
| Corrections | up to 2, and normalisation shifts | 1 | none |
| Allowed $n$ | any | odd | odd, $n < 2^{\bits - 2}$ |

The number of products is the same as with
[Möller–Granlund](algorithms/arithmetic/moller-granlund) reduction, but
Montgomery reduction needs no normalisation shifts and fewer corrections.
This shortens the chain of dependent operations, which matters most in
exponentiation and other long chains of dependent products. The
conversions cost about one ordinary modular multiplication each, so
Montgomery arithmetic pays off only when many operations are done between
them.

## Use in FLINT {#use}

- `nmod_pow_ui` switches to Montgomery arithmetic for odd $n$ and
  exponents of at least $2^{11}$, and `nmod_2_pow_ui` for exponents of at
  least $2^{20}$. Both use the lazy variant when $n < 2^{\bits - 2}$.
- `n_is_prime` performs its strong probable prime tests for $n \ge 2^{32}$ in
  Montgomery form, and `n_factor_pollard_brent_single` iterates
  $x \mapsto x^2 + c$ in Montgomery form.
- The multimodular algorithms for Bell numbers
  (`arith_bell_number_multi_mod`) and Euler numbers (in
  `arb_euler_number_ui`) use lazy Montgomery arithmetic modulo each prime.
- `nmod_sqrt` does all its arithmetic in Montgomery form, using the
  `half_fast`, `fast` or plain functions depending on the size of $n$.
- `_nmod_vec_scalar_mul_nmod` uses `nmod_redc_mul` for full-word odd
  moduli and vectors of length at least 8, where [Shoup's method](algorithms/arithmetic/shoup) is not
  available.
- `nmod_poly` computes Euclidean gcds and resultants of medium length in a
  `gr_ctx_init_nmod_redc_fast` context when $2^{\bits/2 - 1} \le n < 2^{\bits - 2}$.

::: functions nmod_redc_mul nmod_redc_fast_mul nmod_redc_half_mul n_ll_redc n_redc :::
