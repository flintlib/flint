---
title: Binary splitting and product trees
summary: Products of many integers by balanced product trees, sums of hypergeometric series by binary splitting, and simultaneous reduction and Chinese remaindering by product and remainder trees; how FLINT uses them for factorials, primorials, vector products, Chinese remaindering, radix conversion, constants and series.
complexity: A product of integers of total size $N$ bits in $O(M(N) \log k)$ for $k$ factors; $N$ bits of a hypergeometric series in $O(M(N) \log^2 N)$; reduction modulo $k$ moduli of total size $N$ and the inverse CRT in $O(M(N) \log k)$.
order: 50
---

Many computations with large integers produce the result as a product or
a sum of many small terms: a factorial is a product of $n$ word-size
integers, a partial sum of the series of $e$ or $\pi$ is a sum of
rational numbers with growing numerators and denominators, and Chinese
remaindering combines residues modulo many primes. Done term by term,
such computations multiply a large number by a small one at every step
and cost $\Theta(N^2)$ for an $N$-bit result. Organising the work as a
balanced binary tree, so that operands of similar size are multiplied
at each level, reduces the cost to $O(M(N) \log N)$ or $O(M(N) \log^2
N)$, where $M(N)$ is the cost of an $N$-bit multiplication; with FFT
multiplication this is quasi-linear.

## Product trees {#product-trees}

To compute $P = x_1 x_2 \cdots x_k$, split the list in halves, compute
the two partial products recursively and multiply them.

::: proposition {#prop-product-tree}
Let $x_1, \ldots, x_k$ be positive integers with $\sum_i \log_2 x_i \le
N$. The balanced product tree computes $\prod x_i$ in $O(M(N) \log k)$
bit operations, assuming that $M(n)/n$ is nondecreasing. The sequential
product costs $\Theta(Nk')$ word operations, where $k'$ is the number of
factors with at least one word, which is $\Theta(N^2)$ for word-size
factors.
:::

::: proof
At each level of the tree the operands are partial products of disjoint
sublists, so their total size is at most $N + O(k)$ bits. A level
multiplies pairs of operands of sizes $n_j$ and $n'_j$ with $\sum_j (n_j
+ n'_j) \le N$, at a cost of $\sum_j M(n_j + n'_j) \le M(N)$, by the
monotonicity of $M(n)/n$. There are $\lceil \log_2 k \rceil$ levels.
Sequentially, the $i$-th step multiplies a partial product of up to $N$
bits by a factor, which costs $\Theta(N/w)$ word operations for a factor
of one word $w$, and there are $k'$ such steps.
:::

In practice the leaves of the tree are blocks of factors multiplied
sequentially, since multiplying a number of a few words by one word is
cheapest done directly, and several small factors are packed into one
word before they enter the tree. If the factors have very different
sizes, a balanced split by count is not optimal; combining the two
smallest current operands, as in the construction of a Huffman code,
keeps the products balanced by size.

## Binary splitting for series {#series}

A sum of the form
$$
S(a, b) = \sum_{k=a}^{b-1} \frac{p(k)}{q(k)} \prod_{j=a}^{k-1} \frac{r(j)}{q(j)}
$$
with integer-valued $p$, $q$, $r$ covers the hypergeometric series used
for $e$, $\pi$, logarithms of rationals, Euler's constant and many
special functions at rational points [@Bre1976] [@HP1998]. Define, for a
range $[a, b)$,
$$
Q(a, b) = \prod_{k=a}^{b-1} q(k), \qquad
R(a, b) = \prod_{k=a}^{b-1} r(k), \qquad
T(a, b) = Q(a, b)\, S(a, b).
$$

::: theorem {#thm-bsplit}
$T(a, b)$ is an integer, and for $a < m < b$,
$$
Q(a, b) = Q(a, m) Q(m, b), \quad R(a, b) = R(a, m) R(m, b), \quad
T(a, b) = T(a, m) Q(m, b) + R(a, m) T(m, b).
$$
:::

::: proof
Splitting the sum at $m$,
$$
S(a, b) = S(a, m) + \frac{R(a, m)}{Q(a, m)} S(m, b),
$$
because the terms with $k \ge m$ carry the factor $\prod_{j=a}^{m-1}
r(j)/q(j)$ in addition to the factors of $S(m, b)$. Multiplying by $Q(a,
b) = Q(a, m) Q(m, b)$ gives the formula for $T$. Integrality follows by
induction from $T(k, k+1) = p(k)$.
:::

Applying the theorem recursively evaluates $S(0, n) = T(0, n)/Q(0, n)$
with a tree of $O(n)$ products, followed by one division. If $p$, $q$,
$r$ are polynomials, $T$, $Q$ and $R$ over a range of length $\ell$ have
$O(\ell \log n)$ bits, so the products at each level have total size
$O(n \log n)$ bits. For a series that converges geometrically, $n =
O(N)$ terms give $N$ bits, and the cost is $O(M(N \log N) \log N) =
O(M(N) \log^2 N)$ bit operations. The product $R$ is only needed for
left children and is skipped along the right spine of the tree; when $r$
divides $p$, the leaves can be evaluated with smaller numbers. Common
factors of $T$ and $Q$, such as powers of two, can be removed at each
node.

For example, $e = \sum_{k \ge 0} 1/k!$ has $p(k) = 1$, $q(k) = k$ for
$k \ge 1$ (with $q(0) = 1$) and $r(j) = 1$, so $Q(a, b) = \prod q(k)$ and
$T(a, b) = T(a, m) Q(m, b) + T(m, b)$; about $N / \log_2 N$ terms give
$N$ bits.

## Remainder trees and the Chinese remainder theorem {#remainder-trees}

Product trees also speed up the two directions of multimodular
arithmetic [@BorMoe1974] [@Ber2008]. Let $m_1, \ldots, m_k$ be pairwise
coprime with product $M$.

**Reduction.** To compute $x \bmod m_i$ for all $i$, build the product
tree of the moduli and reduce top-down: reduce $x$ modulo the product at
the root, then each remainder modulo the products of the two children,
and so on down to the leaves. Each level costs $O(M(N))$ for $N =
\log_2 M$, so the total is $O(M(N) \log k)$, instead of $O(Nk)$ word
operations for $k$ separate reductions of a number of $N$ bits by
word-size moduli. Bernstein's scaled remainder trees [@Ber2004] replace
the divisions by multiplications with precomputed approximations.

**Reconstruction.** With $M_i = M / m_i$ and $c_i = M_i^{-1} \bmod m_i$,
the integer
$$
x = \sum_{i=1}^{k} (r_i c_i \bmod m_i) M_i
$$
satisfies $x \equiv r_i \pmod{m_i}$ for all $i$. The sum is evaluated
bottom-up in the product tree: a node with children $L$ and $R$ combines
their values as $x_L \cdot \prod_{R} m + x_R \cdot \prod_{L} m$, which is
[](#thm-bsplit) with all ratios equal to 1. A single reduction modulo $M$
at the end puts $x$ in the desired range. The cost is again $O(M(N) \log
k)$.

## Use in FLINT {#use}

**Product trees.**

- `_fmpz_vec_prod` splits the vector in halves down to three entries,
  and `_fmpz_ui_vec_prod` down to 15 words, which are multiplied in
  sequence.
- `fmpz_primorial` multiplies the primes up to $n$ with a balanced tree
  over limb arrays (`mpn_prod_limbs_balanced`), with sequential products
  of fewer than 50 factors at the leaves.
- `_fmpz_rfac_ui`, behind `fmpz_rfac_ui` and `fmpz_rfac_uiui`, splits the
  range of factors in halves and packs as many consecutive factors into
  one word as fit at the leaves.
- `fmpz_fac_ui` (for $n \ge 21$) calls GMP's `mpz_fac_ui`, which uses
  Luschny's prime-swing algorithm and multiplies the factors in product
  trees [@Lus2008] [@GMP]; `fmpz_bin_uiui` calls GMP's `mpz_bin_uiui`.
- `_gr_vec_product` in the generic rings uses a parallel balanced product
  tree (`flint_parallel_binary_splitting`) for long vectors, and the
  polynomial analogue of product trees, the subproduct tree, is described
  in [](algorithms/polynomials/subproduct-tree).
- `fmpz_set_str` converts long decimal strings by splitting them in halves
  and combining the parts with powers of $10 = 2 \cdot 5$, in parallel;
  [](algorithms/integers/radix-conversion) describes radix conversion in
  both directions.

**Huffman trees.** `fmpz_multi_CRT_precompute` and
`fmpz_multi_mod_precompute` build their straight-line programs by
repeatedly combining the two smallest moduli, which keeps the products
balanced by size when the moduli have different sizes ([](fmpz/crt)).

**Remainder trees and CRT for word-size moduli.** `fmpz_comb_t` (behind
`fmpz_multi_mod_ui`, `fmpz_multi_CRT_ui` and `_fmpz_vec_multi_CRT_ui`)
wraps `flint_mpn_crt_t` of `mpn_extras`, which refines the scheme above
in three ways. Small primes are batched into leaves of one word (for
example pairs of 32-bit primes), and the leaf residues are split and
combined with word arithmetic. Consecutive leaves are grouped into
chunks of up to 2048 bits, which are the leaves of the product tree;
within a chunk, the reconstruction is a linear combination with
precomputed multipliers and the reduction is a dot product of the limbs
of the input with precomputed powers $2^{\bits j} \bmod l$. In the
reconstruction, the nodes are combined as above without intermediate
reductions, the values growing by at most one bit per level, and a single
reduction modulo the product of all primes is done at the end. Large
nodes are processed in parallel; [](algorithms/integers/multimodular-reduction)
describes this design in detail. The multimodular algorithms of
`fmpz_mat` (multiplication, characteristic polynomial), `fmpz_poly`
(GCD, resultant, interpolation, Taylor shift), `arith` and `bernoulli`
use this code for their Chinese remaindering.

**Binary splitting for series.**

- `mp_real_hypgeom_series` and `mp_real_hypgeom_series_int64` implement
  binary splitting for general hypergeometric series in the format of
  [](#thm-bsplit), with leaves of several terms evaluated in limb
  arithmetic and the final division done once. The constants $\pi$
  (Chudnovsky's series), $e$, $\log 2$ and Euler's constant, and the
  `arb` constants `arb_const_pi` and `arb_const_e` that call them, are
  computed this way, as are tables of $\log(1 + 2^{-i})$ and $\arctan
  2^{-i}$ for the elementary functions.
- The `arb` elementary functions use binary splitting of the Taylor
  series of $\exp$, $\sin$, $\cos$ and $\arctan$ at high precision
  (`_arb_exp_sum_bs_powtab`, `_arb_atan_sum_bs_powtab`), with a table of
  the powers of the argument needed at each level of the tree.
- `fmpq_harmonic_ui` computes the harmonic number $H_n = \sum_{k \le n}
  1/k$ exactly by binary splitting over the odd denominators only, each
  weighted by the contribution of its multiples by powers of two.
- `arb_hypgeom_sum_fmpq_arb` and `acb_hypgeom_pfq_sum_bs` sum
  hypergeometric series with rational parameters,
  `arb_hypgeom_rising_ui_bs` computes rising factorials, and
  `arith_stirling_number_1u` and related functions expand
  $(x + a)(x + a + 1) \cdots (x + b - 1)$ by binary splitting over
  polynomials.

::: functions _fmpz_vec_prod fmpz_primorial fmpz_rfac_ui fmpz_multi_CRT fmpz_multi_mod_ui :::
