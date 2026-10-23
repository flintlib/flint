---
title: Multimodular reduction and Chinese remaindering
summary: Converting an integer to and from its residues modulo many word-size primes in quasi-linear time with a subproduct tree, using fractional cofactors so that Chinese remaindering needs no large inverses and no intermediate reductions.
complexity: $O(M(N) \log N)$ for $N$-bit products, for the precomputation and for each conversion.
order: 45
---

A multimodular algorithm computes a result $x$ (an integer, or the entries
of a matrix or the coefficients of a polynomial) by reducing its inputs
modulo many word-size primes $m_1, \ldots, m_k$, solving the problem
modulo each prime, and recovering $x$ from its residues by the Chinese
remainder theorem, which determines $x$ modulo $P = m_1 \cdots m_k$. The
two conversions, *multi-mod* (from $x$ to its residues) and *multi-CRT*
(back), each cost $O(k^2)$ limb operations if done naively; the subproduct tree makes them quasi-linear
[@BorMoe1974; @Ber2008].

The generic and fixed-size Chinese remaindering used inside the
[small-prime FFT](algorithms/fft/small-prime-fft) is described with
[`crt_helpers`](crt_helpers); multimodular algorithms for matrices and
polynomials are discussed in [](algorithms/linear-algebra/multimodular)
and [](algorithms/fft/multimodular).

## The subproduct tree {#tree}

Group the moduli into the leaves of a balanced binary tree, and label
each node with the product $v$ of the moduli below it, so that the root
carries $P$. With $N$ the bit size of $P$, each level of the tree has
total size about $N$ bits, there are $O(\log k)$ levels, and the tree is
computed bottom-up in time $O(M(N) \log k)$.

::: algorithm Remainder tree (multi-mod) {#alg-mod}
**Input:** $x \ge 0$ and the subproduct tree.

**Output:** $x \bmod m_i$ for all $i$.

1. At the root, $r \gets x \bmod P$.
2. Going down, for a node with residue $r$ and children with products
   $v_b$ and $v_c$, set $r_b \gets r \bmod v_b$ and $r_c \gets r \bmod v_c$.
3. At the leaves, return the residues.
:::

Since $v_b \mid v$, $(x \bmod v) \bmod v_b = x \bmod v_b$, so the
algorithm is correct; each level costs divisions of total size $O(N)$, so
the whole costs $O(M(N) \log k)$ with fast division.

## Chinese remaindering with fractional cofactors {#crt}

The classical solution $x \equiv \sum_i r_i \, c_i \, (P/m_i) \pmod P$ with
$c_i = (P/m_i)^{-1} \bmod m_i$ needs the cofactors $P/m_i$, of $N$ bits
each. The tree form combines the children $b$ and $c$ of a node $a$ with
$v_a = v_b v_c$ by $A = A_b v_c + A_c v_b$, which never forms anything
larger than $v_a$. For this to work the values must carry the right
scaling.

::: lemma {#lem-combine}
Let $x$ be an integer, and for a node $a$ let $R_a = (P/v_a) \bmod v_a$
(invertible modulo $v_a$ if the moduli are pairwise coprime). If
$A_b \equiv R_b^{-1} x \pmod{v_b}$ and $A_c \equiv R_c^{-1} x \pmod{v_c}$,
then $A = A_b v_c + A_c v_b$ satisfies $A \equiv R_a^{-1} x \pmod{v_a}$.
At the root, $R = 1$ and $A \equiv x \pmod P$.
:::

::: proof
Modulo $v_b$, $A \equiv A_b v_c \equiv R_b^{-1} v_c\, x$. Now
$P/v_b = (P/v_a) v_c$, so $R_b \equiv R_a v_c \pmod{v_b}$ and
$R_b^{-1} v_c \equiv R_a^{-1} \pmod{v_b}$. The same holds modulo $v_c$,
and since $v_b$ and $v_c$ are coprime, $A \equiv R_a^{-1} x$ modulo their
product $v_a$.
:::

The values $R$ themselves are computed from the root down by a remainder
tree: $R_b = \bigl((R_a \bmod v_b) \, v_c\bigr) \bmod v_b$, which follows
from $R_b \equiv R_a v_c \pmod{v_b}$ and $v_b \mid v_a$. At the bottom of
the tree only the inverses $R^{-1}$ modulo the leaves are needed, and a
leaf product fits in a word, so no multiprecision inverse is ever
computed. The precomputation therefore costs $O(M(N) \log k)$: a product
tree and a remainder tree.

::: algorithm Multi-CRT {#alg-crt}
**Input:** residues $x_i = x \bmod m_i$, the subproduct tree and the
leaf cofactor data.

**Output:** $x \bmod P$.

1. At each leaf group (a chunk of consecutive moduli with product $M$),
   form $y = \sum_j \bigl(x_j c_j \bmod m_j\bigr) (M/m_j)$ with
   precomputed $c_j = (M/m_j)^{-1} R^{-1} \bmod m_j$. Then
   $y \equiv R^{-1} x \pmod M$ and $y < k_M M$ for $k_M$ moduli in the
   chunk.
2. Combine upwards with [](#lem-combine), without reducing:
   $A = A_b v_c + A_c v_b$.
3. Reduce the root value modulo $P$; for a signed result, subtract $P$ if
   the value exceeds $P/2$.
:::

::: proposition {#prop-growth}
In [](#alg-crt), every value $A$ at a node with product $v$ satisfies
$A < \ell\, v$, where $\ell$ is the number of moduli below the node. In
particular the root value is less than $kP$, and the final reduction is a
division with a quotient of at most $\log_2 k$ bits.
:::

::: proof
By induction: at a chunk, each term of $y$ is less than $m_j (M/m_j) = M$.
If $A_b < \ell_b v_b$ and $A_c < \ell_c v_c$, then
$A_b v_c + A_c v_b < (\ell_b + \ell_c) v_b v_c$.
:::

So the values grow by about one bit per level, and the combination is a
pair of unbalanced products and one addition per node, again
$O(M(N) \log k)$ in total.

## Refinements in FLINT {#refinements}

- **Batched leaves.** Moduli of at most 32 bits are multiplied in pairs,
  and 20-bit moduli in triples, so that the tree works with full-word
  leaves; single-word code splits and combines the individual moduli at
  the bottom.
- **Basecases.** Below a chunk of about 2048 bits the CRT is the linear
  sum of step 1 (by `mpn_addmul_1`), and below a node of about 2048 bits
  the reduction computes $x \bmod m_j$ for all leaves at once as a dot
  product of the limbs of $x$ with precomputed powers $B^i \bmod m_j$,
  with a two-limb accumulator and a single reduction per leaf.
- **Fixed-length code.** When $P$ is a single small chunk, the CRT uses
  the straight-line kernels of `crt_helpers.h`.
- **Small values.** Results are often much smaller than $P$ (for example
  the small entries of a matrix whose bound is dominated by a few large
  ones). The reconstruction first reconstructs $x$ modulo the product $Q$
  of the first few moduli. If the result $x'$ has at least a limb of slack
  below $Q$ and $x' \bmod m_i = x_i$ for all the other moduli, then
  $x' \equiv x \pmod P$ and $0 \le x' < Q \le P$, so $x' = x$ and the tree
  is not traversed. A random residue passes the slack test with
  probability about $2^{-64}$, so this costs almost nothing when it fails.
- **Division.** Reductions down the tree use precomputed inverses of the
  node products (`flint_mpn_mod_preinvn`) above 16 limbs.
- **Vectors.** For many values, the tables are traversed once per block of
  values (`flint_mpn_multi_mod_vec`, `flint_mpn_multi_crt_vec`).
- **Memory.** The structure stores $O(N \log k)$ bits. The one-shot
  reconstruction `flint_mpn_multi_crt_once` instead traverses the tree
  depth first, recomputing subproducts as needed, with $O(N)$ memory,
  which makes reconstructions of billions of bits possible. Large nodes
  are processed in parallel.

## Use in FLINT {#use}

The functions are `flint_mpn_crt_init`, `flint_mpn_multi_mod`,
`flint_mpn_multi_crt` and their variants in `mpn_extras`
([](mpn_extras/multi_crt)). They are used through `fmpz_comb_t`
(`fmpz_multi_mod_ui`, `fmpz_multi_CRT_ui`) by the multimodular algorithms
for `fmpz_mat` and `mpn_mod` matrices (multiplication), `fmpz_poly` and
`fmpq_poly` (GCDs, interpolation, Taylor shifts), resultants of bivariate
polynomials, and by the multimodular computation of Bernoulli, Euler and
Bell numbers and Stirling numbers in `arith`.

::: functions flint_mpn_crt_init flint_mpn_multi_mod flint_mpn_multi_crt flint_mpn_multi_crt_once :::
