---
title: Simultaneous inversion
summary: Inverts $k$ elements modulo $n$ with one modular inversion and $3(k - 1)$ multiplications, using prefix products (Montgomery's trick).
complexity: One inversion and $3(k-1)$ multiplications modulo $n$; $k$ words of temporary storage.
order: 40
---

A modular inversion by the extended Euclidean algorithm costs far more than
a multiplication. To invert many elements at once, Montgomery
[@Mon1987] observed that one inversion suffices: invert the
product of all elements and recover the individual inverses with
multiplications. See also [@BZ2011, §2.5.1].

## The algorithm {#algorithm}

::: algorithm Simultaneous inversion {#alg-batch}
**Input:** $v_0, \ldots, v_{k-1}$ invertible modulo $n$, with $k \ge 1$.

**Output:** $r_j = v_j^{-1} \bmod n$ for $0 \le j < k$.

1. $p_0 \gets v_0$; for $j = 1, \ldots, k - 1$: $p_j \gets p_{j-1} v_j \bmod n$.
2. $u \gets p_{k-1}^{-1} \bmod n$.
3. For $j = k - 1, \ldots, 1$: $r_j \gets u \, p_{j-1} \bmod n$, then
   $u \gets u \, v_j \bmod n$.
4. $r_0 \gets u$.
:::

::: theorem {#thm-batch}
[](#alg-batch) is correct. It performs one inversion and $3(k - 1)$
multiplications modulo $n$.
:::

::: proof
After step 1, $p_j = v_0 v_1 \cdots v_j$, which is invertible because
every factor is. We show that $u = p_j^{-1}$ at the start of iteration $j$
of step 3. This holds for $j = k - 1$ by step 2. If it holds for $j$, then
$r_j = p_j^{-1} p_{j-1} = v_j^{-1}$, and the update gives
$u = p_j^{-1} v_j = p_{j-1}^{-1}$, which is the claim for $j - 1$. After
the last iteration $u = p_0^{-1} = v_0^{-1}$, which step 4 stores. Step 1
uses $k - 1$ multiplications and step 3 uses $2(k - 1)$.
:::

If some $v_j$ is not invertible, then neither is $p_{k-1}$, and step 2
fails without revealing which element is at fault. Step 3 reads $v_j$
after writing $r_j$, so the output cannot overwrite the input.

## Use in FLINT {#use}

`_nmod_vec_invert` implements [](#alg-batch). `_nmod_vec_invert_generic`
stores the prefix products and multiplies with `nmod_mul`.
`_nmod_vec_invert_shoup`, used when $n < 2^{\bits - 1}$, stores with every
$p_j$ its precomputed [Shoup](algorithms/arithmetic/shoup) quotient
$\lfloor p_j 2^{\bits} / n \rfloor$, and with every $v_j$ its quotient
too. In step 1, the quotient of $v_j$ costs one division, and
`n_mulmod_and_precomp_shoup` computes $p_j$ together with its quotient
from those of $p_{j-1}$ and $v_j$, without a further division. All
multiplications in step 3 are then Shoup multiplications by $p_{j-1}$ and
$v_j$. This
costs $3k$ words of temporary storage. The inversion in step 2 is
`nmod_inv`, which raises a `FLINT_IMPINV` exception if $p_{k-1}$ is not
invertible.
