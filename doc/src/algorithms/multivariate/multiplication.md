---
title: Sparse and dense multiplication
summary: Multiplication of sparse multivariate polynomials with Johnson's heap algorithm and the improvements of Monagan and Pearce (packed monomials, chaining, small heaps, three-word accumulation, parallel splitting of the output), dense array accumulation, Kronecker substitution for dense inputs, and powering by the FPS recurrence.
complexity: Heap method $O(\ell_B \ell_C \log \min(\ell_B, \ell_C))$ monomial comparisons and $\ell_B \ell_C$ coefficient multiplications with $O(\ell_B + \ell_C)$ extra memory besides the output; array method $O(\ell_B \ell_C + D)$ for an array of size $D$; Kronecker substitution $O(M(D))$ for the dense box of size $D$; FPS powering $O(k \ell_B \ell_P \log \ell_B)$ for $P = B^{k-1}$.
order: 20
---

Let $B$ and $C$ be polynomials with $\ell_B$ and $\ell_C$ terms in the
distributed representation of [](algorithms/multivariate/representation).
Their product has at most $\ell_B \ell_C$ terms, and in the sparse case
about that many; in the dense case far fewer, because many products
$b_i c_j$ share a monomial. The *classical* algorithm forms all products
and sorts or merges them. Done naively (adding $b_i C$ to an accumulator
one row at a time) it costs $O(\ell_B^2 \ell_C)$ monomial comparisons in
the worst case. FLINT's algorithms avoid this in different ways: a heap
(sparse inputs), a dense array (moderately dense inputs) or a reduction to
univariate multiplication (dense inputs).

## The heap method {#heap}

Johnson's algorithm [@Joh1974] merges the $\ell_B$ sorted sequences
$b_i C = (b_i c_0, b_i c_1, \ldots)$ with a priority queue holding the
next unused product of each sequence. Monagan and Pearce [@MonPea2007]
made it practical with three refinements: packed monomials (so that the
heap keys are words compared with one instruction), *chaining* (products
with equal monomials are linked in one heap node, so the heap only holds
distinct monomials) and a rule that keeps the heap small. We describe the
algorithm as implemented, with rows indexed by $B$ and columns by $C$.

::: algorithm Heap multiplication {#alg-heap}
Input: $B = \sum_{i < \ell_B} b_i x^{\beta_i}$ and $C = \sum_{j < \ell_C}
c_j x^{\gamma_j}$, sorted.
Output: $BC$, sorted.

1. Insert $(0, 0)$ with key $x^{\beta_0 + \gamma_0}$ into an empty heap.
   Mark row 0 as being in the heap at column 0.
2. While the heap is not empty:
   1. Let $m$ be the largest key. Pop all nodes $(i, j)$ with key $m$
      (each heap entry is a chain of such nodes), add $\sum b_i c_j$ to an
      accumulator and record the popped pairs.
   2. For each popped $(i, j)$: insert $(i + 1, j)$ if $i + 1 < \ell_B$
      and row $i + 1$ has used exactly the columns before $j$ and has no
      node in the heap; and insert $(i, j + 1)$ if $j + 1 < \ell_C$ and
      either $i = 0$ or row $i - 1$ has already used column $j + 1$.
   3. If the accumulator is nonzero, append the term with monomial $m$.
:::

::: lemma Heap invariant {#lem-heap}
In [](#alg-heap), the heap always contains the largest product not yet
used, each row has at most one node in the heap, and every product $b_i
c_j$ is used exactly once.
:::

::: proof
For $i \ge 1$, $x^{\beta_{i-1} + \gamma_j} > x^{\beta_i + \gamma_j}$
because the ordering is compatible with multiplication; similarly
$x^{\beta_i + \gamma_{j-1}} > x^{\beta_i + \gamma_j}$. So the products
form a matrix decreasing along rows and columns, and the unused product
$(i, j)$ with largest monomial must have both $(i-1, j)$ and $(i, j-1)$
used (when they exist). The insertion rule inserts $(i, j)$ exactly when
the second of these two neighbours is popped, so every such *frontier*
product is in the heap, and in particular the largest unused one. A row
$i$ gets a new node only when its previous node $(i, j - 1)$ is popped,
or when its first node $(i, 0)$ is inserted from row $i - 1$; hence each
row has at most one node and every product is inserted exactly once.
:::

::: theorem Cost of the heap method {#thm-heap}
[](#alg-heap) computes $BC$ with $\ell_B \ell_C$ coefficient
multiplications, at most $\ell_B \ell_C \lceil \log_2(\ell_B + 1) \rceil$
monomial comparisons in the heap, and $O(\ell_B)$ memory besides the
inputs and output. FLINT exchanges $B$ and $C$ if necessary so that
$\ell_B = \min(\ell_B, \ell_C)$.
:::

::: proof
By [](#lem-heap) the heap has at most $\ell_B$ entries, so each insertion
and extraction costs $O(\log \ell_B)$ comparisons, and there are $\ell_B
\ell_C$ of each (chained products cost less). The output is produced in
decreasing order because the popped key is always the largest unused
monomial.
:::

In the dense case, chaining makes the cost per product close to constant:
many products are popped together and the heap stays much smaller than
$\ell_B$. The packed monomial of each candidate product is a word sum
($N = 1$) or an $N$-word sum, and the comparisons use [the packed
comparison](algorithms/multivariate/representation#lem-cmp).

### Coefficient accumulation {#accumulation}

Over $\Z$, if all coefficients of $B$ and $C$ are small `fmpz`
(absolute value below $2^{62}$ on a 64-bit machine), each product $b_i
c_j$ is a signed two-word integer, and a sum of at most $\ell_B$ of them
has absolute value below $\ell_B 2^{124} < 2^{190}$; it is accumulated in
three words with `smul_ppmm` and three-word additions and converted to an
`fmpz` once per output term. Otherwise `fmpz_addmul` is used. Over
$\Z/n\Z$ (`nmod_mpoly`), products are accumulated unreduced in up to three
words and reduced once per output term, the delayed reduction of
[](algorithms/arithmetic/choosing). This accounts for much of the
speed of the specialised versions compared with the generic one.

### Parallel heap multiplication {#parallel}

Monagan and Pearce [@MonPea2009] parallelise by splitting the *output*
into ranges of monomials, each computed by a heap restricted to the
products in the range. FLINT (`fmpz_mpoly_mul_heap_threaded`) cuts the
monomials at $4t$ boundaries for $t$ threads: the $i$-th boundary is a
monomial $e$ such that about $((i+1)/4t)^2 \, \ell_B \ell_C$ of the
products $b_k c_l$ are at most $e$ (found by `mpoly_search_monomials`,
a binary search along the rows for each candidate). The quadratic spacing
makes the ranges for the largest monomials, which are processed first,
contain more products, so that no thread ends long before the others. For
each range and each row of $B$ the first and last columns of $C$ in the
range are known, and the restricted heap algorithm is run on these
segments. The pieces are concatenated in order at the end.

## Dense methods {#dense}

### Array accumulation {#array}

If the product is dense in a box of moderate size, it is faster to add the
products into a dense array indexed by the exponents than to sort them.
FLINT's array method (`fmpz_mpoly_mul_array`) slices the product:

- for `ORD_LEX`, by the exponent of $x_0$: with $d_i = \deg_{x_i} B +
  \deg_{x_i} C + 1$, the product of the slices of $B$ and $C$ of degrees
  $s$ and $t$ in $x_0$ contributes to slice $s + t$, and each slice is a
  dense array of $\prod_{i \ge 1} d_i \le 300000$ entries, with $d_0 \le
  300$;
- for the degree orderings, by total degree: within a slice the exponent
  of $x_{n-1}$ is determined by the others, and the array is indexed by
  the exponents of $x_0, \ldots, x_{n-2}$, each below a bound $d$ on the
  total degree of the product, so it has $d^{n-1}$ entries.

Each output slice is computed by adding all contributing products of
input slices into a zeroed array and then reading off the nonzero
entries in decreasing order.

::: lemma Accumulator size {#lem-array-bound}
If the products contributing to an output slice come from pairs of input
slices $(B_s, C_t)$, $s + t$ fixed, then every coefficient of the slice is
bounded in absolute value by
$$
\sum_{s + t \text{ fixed}} \min\bigl(\|B_s\|_1 \|C_t\|_\infty,
\|B_s\|_\infty \|C_t\|_1\bigr).
$$
:::

::: proof
A coefficient of $B_s C_t$ is $\sum b_i c_j$ over pairs with a fixed
monomial; each $b_i$ appears at most once in the sum, and so does each
$c_j$, so its absolute value is at most both $\|B_s\|_1 \|C_t\|_\infty$
and $\|B_s\|_\infty \|C_t\|_1$. Sum over the pairs.
:::

FLINT evaluates this bound in bits for each slice and accumulates in
one, two or three words if it fits (with one bit for the sign), and in
`fmpz` otherwise. Slices are distributed over threads in
`fmpz_mpoly_mul_array_threaded`. The method declines (returns 0) when
the array would be too large or the exponents do not fit one word.

### Kronecker substitution {#kronecker}

For dense inputs of moderate degree the fastest method is to substitute a
power of a single variable for each variable and multiply univariate
polynomials [@Kro1882].

::: lemma Kronecker substitution {#lem-kronecker}
Let $d_i > \deg_{x_i} B + \deg_{x_i} C$ for $i \ge 1$, and let $\kappa$ be
the substitution $x_{n-1} \mapsto y$, $x_{n-2} \mapsto y^{d_{n-1}}$,
$\ldots$, $x_0 \mapsto y^{d_1 \cdots d_{n-1}}$. Then $\kappa$ is injective
on the monomials of $BC$, so $BC$ is recovered from $\kappa(B)\kappa(C)$
by reading the exponents of $y$ in mixed radix $(d_1, \ldots, d_{n-1})$.
:::

::: proof
$\kappa(x^e) = y^{E}$ with $E = e_{n-1} + d_{n-1}(e_{n-2} + d_{n-2}(\cdots
+ d_1 e_0))$. For the monomials of $BC$, $0 \le e_i < d_i$ for $i \ge 1$,
so $E$ determines $e$ by division with remainder. Since $\kappa$ is a
ring homomorphism, $\kappa(BC) = \kappa(B)\kappa(C)$, and no two terms of
$BC$ are mapped to the same power of $y$.
:::

`fmpz_mpoly_mul_dense` converts $B$ and $C$ to dense arrays in this
order (no padding is needed for $x_0$), multiplies them with
`fmpz_poly_mul`, which uses the univariate algorithms of
[](algorithms/polynomials/overview) (Kronecker substitution into
integers, multimodular and Schönhage–Strassen methods), and converts back.
The cost is that of a univariate product of length $D = \prod_i d_i$, so
the method only pays when $D$ is not much larger than the number of
terms of the product. For `nmod_mpoly` the same is done with
`nmod_poly_mul`.

## Squaring and powering {#powering}

Squaring with the heap method can use symmetry: $B^2 = \sum_i b_i^2
x^{2\beta_i} + 2 \sum_{i < j} b_i b_j x^{\beta_i + \beta_j}$ needs only
about half the products (`gr_mpoly_sqr_commutative_heap`). For higher
powers, repeated multiplication computes $B^2, B^3, \ldots$ with a cost
dominated by the last product, and binary powering does not help much for
sparse polynomials because the last squaring of a dense intermediate
result is as expensive as all the rest. Monagan and Pearce
[@MonPea2012] instead compute $B^k$ directly by the recurrence for the
coefficients of a power of a power series (*FPS*, after the classical
method of Euler; see [@Knu1997, §4.7]).

::: lemma FPS recurrence {#lem-fps}
Let $D$ be a derivation of $R[x_0, \ldots, x_{n-1}]$ over a ring $R$ with
$D(x^e) = \phi(e) x^e$ for a group homomorphism $\phi: \Z^n \to \Z$
(for example $\phi(e) = \sum_i w_i e_i$). If $P = B^{k-1}$ and $A = B^k$,
then $B \, D(P) = (k-1) \, P \, D(B)$ and $A = BP$.
:::

::: proof
For a derivation, $D(B^{k-1}) = (k-1) B^{k-2} D(B)$. Multiply by $B$.
:::

Write $B = b_0 x^{\beta_0} + B'$ and compare the coefficients of a
monomial $x^{\mu}$ on both sides of $B D(P) - (k-1) P D(B) = 0$:
$$
\sum_{i,j:\ \beta_i + \pi_j = \mu} \bigl(\phi(\pi_j) - (k-1)\phi(\beta_i)\bigr)
\, b_i p_j = 0,
$$
where $P = \sum_j p_j x^{\pi_j}$. The term with $i = 0$ is $(\phi(\mu -
\beta_0) - (k-1)\phi(\beta_0)) b_0 p_{\mu - \beta_0} = (\phi(\mu) -
k\phi(\beta_0)) b_0 p_{\mu - \beta_0}$. So if $\phi(\mu) \ne k
\phi(\beta_0)$, the coefficient of $x^{\mu - \beta_0}$ in $P$ is
$$
p_{\mu - \beta_0} = \frac{1}{(k\phi(\beta_0) - \phi(\mu))\, b_0}
\sum_{i \ge 1} \bigl(\phi(\pi_j) - (k-1)\phi(\beta_i)\bigr) b_i p_j ,
$$
in terms of coefficients of $P$ at larger monomials. Generating the
products $b_i p_j$ with a heap in decreasing order of $\beta_i + \pi_j$,
the algorithm computes, for each monomial $\mu$, the coefficient
$a_\mu = \sum b_i p_j$ of $A = BP$ and the coefficient $p_{\mu - \beta_0}$
of $P$ at the same time; the latter is needed only when $x^{\beta_0}$
divides $x^\mu$, which is tested with
[the packed divisibility test](algorithms/multivariate/representation#lem-divides).

::: theorem FPS powering {#thm-fps}
With $\phi$ injective on the monomials that occur (for example $\phi(e) =
\Phi(e)$, the packed exponent read as an integer) and $\phi(\mu) \ne
k \phi(\beta_0)$ for all monomials $\mu \ne k\beta_0$ of $A$, the
algorithm computes $A = B^k$ and $P = B^{k-1}$ with $O(\ell_B
\ell_P)$ coefficient operations and $O(\ell_B \ell_P \log \ell_B)$
monomial comparisons, all divisions being exact over $\Z$.
:::

::: proof
The monomials of $A$ are at most $x^{k\beta_0}$, and those of $P$ at most
$x^{(k-1)\beta_0}$; the recurrence determines the coefficients of $P$ in
decreasing order, starting from $p_{(k-1)\beta_0} = b_0^{k-1}$. Each
product $b_i p_j$ enters the heap once, so the cost is that of a heap
multiplication of $B$ by $P$. The quotients are coefficients of
$P \in \Z[x]$, so the divisions are exact. If $\phi = \Phi$, injectivity
holds by construction, and for the orderings with zero `cmpmask` it is
also order preserving, so $\phi(\mu) < \phi(k\beta_0)$ for $\mu <
k\beta_0$. For `ORD_DEGREVLEX`, injectivity alone gives $\phi(\mu) \ne
\phi(k \beta_0)$.
:::

Since $\ell_P$ is the number of terms of $B^{k-1}$, this is about the cost
of the final multiplication in repeated multiplication, while avoiding the
intermediate products. `fmpz_mpoly_pow_ui` uses FPS for $k \ge 3$ (after
checking that the output length $k(\ell_B - 1) + 1$ cannot overflow); the
divisions by $(k\phi(\beta_0) - \phi(\mu)) b_0$ make it specific to
characteristic zero. Over $\Z/n\Z$, `nmod_mpoly_pow_ui` multiplies
repeatedly (`nmod_mpoly_pow_rmul`), and for prime $n = p$ uses $B^{p} =
\sum b_i x^{p\beta_i}$ (the Frobenius map, since $b^p = b$ in
$\F_p$) to handle the base-$p$ digits of $k$.

## Generic and specialised implementations {#use}

The generic implementation is `gr_mpoly_mul`, for any coefficient ring
given by a `gr` context: it uses `gr_mpoly_mul_monomial` for a monomial
factor, the heap method `gr_mpoly_mul_heap` (and
`gr_mpoly_mul_heap_threaded` when the product has more than $10^4$
term products and threads are available), and the symmetric squaring
`gr_mpoly_sqr_commutative_heap` for $B^2$ in commutative rings. It has no
array or Kronecker methods, and accumulates coefficients with the ring's
own operations. The specialised modules add:

| Module | Heap | Threaded heap | Array | Kronecker | Powering |
|---|---|---|---|---|---|
| `fmpz_mpoly` | `fmpz_mpoly_mul_johnson` | `fmpz_mpoly_mul_heap_threaded` | `fmpz_mpoly_mul_array` (and threaded) | `fmpz_mpoly_mul_dense` | FPS, `fmpz_mpoly_pow_fps` |
| `fmpq_mpoly` | via `fmpz_mpoly` | via `fmpz_mpoly` | via `fmpz_mpoly` | via `fmpz_mpoly` | via `fmpz_mpoly` |
| `nmod_mpoly` | `nmod_mpoly_mul_johnson` | `nmod_mpoly_mul_heap_threaded` | `nmod_mpoly_mul_array` (and threaded) | `nmod_mpoly_mul_dense` | repeated multiplication, Frobenius |
| `fmpz_mod_mpoly` | `fmpz_mod_mpoly_mul_johnson` | | | `fmpz_mod_mpoly_mul_dense` | repeated multiplication |
| `fq_nmod_mpoly` | `fq_nmod_mpoly_mul_johnson` | | | | repeated multiplication |

The choice of method in `fmpz_mpoly_mul`, with its thresholds, is described in
[](fmpz_mpoly/mul); `nmod_mpoly_mul` follows the same pattern.

- Heap: [](#alg-heap) in `_fmpz_mpoly_mul_johnson` (one-word exponents:
  `_fmpz_mpoly_mul_johnson1`), with the accumulation of
  [](#accumulation); parallel splitting in `fmpz_mpoly_mul_heap_threaded`
  with `mpoly_search_monomials`.
- Array: `fmpz_mpoly_mul_array`, `_fmpz_mpoly_mul_array_chunked_LEX` and
  `_fmpz_mpoly_mul_array_chunked_DEG`, the accumulators
  `_fmpz_mpoly_addmul_array1_slong1`, `_fmpz_mpoly_addmul_array1_slong2`,
  `_fmpz_mpoly_addmul_array1_slong` and `_fmpz_mpoly_addmul_array1_fmpz`
  ([](#lem-array-bound)).
- Kronecker: `fmpz_mpoly_mul_dense` ([](#lem-kronecker)).
- Powering: `fmpz_mpoly_pow_fps`, `fmpz_mpoly_pow_ui`,
  `fmpz_mpoly_pow_fmpz` ([](#thm-fps)).

::: functions fmpz_mpoly_mul fmpz_mpoly_mul_johnson fmpz_mpoly_mul_heap_threaded fmpz_mpoly_mul_array fmpz_mpoly_mul_dense fmpz_mpoly_pow_fps gr_mpoly_mul :::
