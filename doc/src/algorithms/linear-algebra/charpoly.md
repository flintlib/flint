---
title: Characteristic and minimal polynomials
summary: Berkowitz's division-free algorithm and Danilevsky's reduction to companion form for the characteristic polynomial over $\Z/n\Z$, and the Krylov method with least common multiples that FLINT uses for the minimal polynomial.
complexity: $O(r^4)$ ring operations for Berkowitz's algorithm, $O(r^3)$ field operations for Danilevsky's algorithm, and $O(r^3)$ field operations in typical cases (at most $O(r^4)$) for the minimal polynomial of an $r \times r$ matrix.
order: 50
---

Let $M$ be an $r \times r$ matrix over $\Z/N\Z$. Its characteristic
polynomial $\chi_M(x) = \det(xI - M)$ can be computed with
$O(r^\omega)$ operations over a field by the algorithms of Keller-Gehrig
[@KG1985] and Pernet and Storjohann [@PS2007], which reduce the problem
to matrix multiplication; see [@DPW2005] for a practical comparison. FLINT
uses two simpler algorithms instead, Danilevsky's cubic reduction to
companion form and Berkowitz's quartic division-free algorithm, chosen so
that the result is correct for every modulus $N$, prime or not. The
minimal polynomial is computed with a Krylov method ([](#minpoly)).

## Generic version {#generic}

Both characteristic polynomial algorithms are implemented once in `gr_mat`
(`_gr_mat_charpoly_berkowitz`, `_gr_mat_charpoly_danilevsky`), together
with the Faddeev-LeVerrier algorithm and its baby-step giant-step
variant (which divide by $1, \ldots, r$ and need characteristic zero or
larger than $r$) and reduction to Hessenberg form. Berkowitz's
algorithm, being division-free, is the default over any commutative
ring. Rings override `GR_METHOD_MAT_CHARPOLY` to use Danilevsky's
algorithm above a dimension: 8 for `nmod`, `nmod8` and `nmod32` and 15
for `fmpz_mod`, with Berkowitz's algorithm as the fallback when a pivot
is not invertible, and 12 for the `fq` family; `fmpz` and `fmpq` use
`fmpz_mat_charpoly`, which reduces modulo primes for large matrices
([](algorithms/linear-algebra/multimodular#charpoly)). `nmod_mat_charpoly`
uses its own implementations of the two algorithms, except in the case
described in [](#danilevsky). The minimal
polynomial is computed by the same Krylov method
over any field (`gr_mat_minpoly_field`), with the row reductions done
by the ring method `GR_METHOD_MAT_REDUCE_ROW`, which `fmpz_mod`,
`mpn_mod` and the `fq` and `fq_nmod` rings implement. Over exact real
and complex numbers, where Danilevsky's pivots require zero tests,
`ca_mat_charpoly` uses Berkowitz's algorithm by default
([linear algebra over exact numbers](algorithms/exact/linear-algebra#division-free)).

## Berkowitz's algorithm {#berkowitz}

Berkowitz's algorithm [@Ber1984] computes $\chi_M$ without any division,
so it works over every commutative ring. For $0 \le t \le r$, let $M_t$
be the leading principal $t \times t$ submatrix of $M$ and $\chi_t$ its
characteristic polynomial, with $\chi_0 = 1$. The algorithm computes
$\chi_1, \chi_2, \ldots, \chi_r$ in turn, each from the previous one and
a few entries of powers of $M_t$.

::: lemma {#lem-border}
Let $s_k = (M_{t+1}^k)_{t,t}$, the bottom right entry of the $k$-th power
of $M_{t+1}$, for $k \ge 0$. Then, as power series in $y$,

$$
y^t \chi_t(1/y) = y^{t+1} \chi_{t+1}(1/y) \cdot \sum_{k \ge 0} s_k y^k .
$$
:::

::: proof
By Cramer's rule, the bottom right entry of $(xI - M_{t+1})^{-1}$ is the
cofactor $\det(xI - M_t)$ divided by $\det(xI - M_{t+1})$, that is
$\chi_t(x)/\chi_{t+1}(x)$. On the other hand,
$(xI - M_{t+1})^{-1} = \sum_{k \ge 0} M_{t+1}^k x^{-k-1}$ as a power
series in $1/x$. Comparing bottom right entries and substituting
$x = 1/y$ gives $\chi_t(1/y) = \chi_{t+1}(1/y) \, y \sum_k s_k y^k$, and
multiplying by $y^t$ gives the claim.
:::

Since $s_0 = 1$, the series $\sum_k s_k y^k$ is invertible without
division, and the reversed polynomial
$\tilde\chi_{t+1}(y) = y^{t+1} \chi_{t+1}(1/y)$ is the product of the
reversed $\tilde\chi_t$ and the inverse series, truncated to degree
$t + 1$. Writing $\tilde\chi_t = \sum_i b_i y^i$ and
$\tilde\chi_{t+1} = \sum_i c_i y^i$, this is the triangular recurrence

$$
c_0 = 1, \qquad c_k = b_k - \sum_{j=1}^{k} s_j \, c_{k-j}
\quad (1 \le k \le t + 1),
$$

with $b_{t+1} = 0$. The numbers $s_1, \ldots, s_{t+1}$ come from the
vectors $u_k = M_{t+1}^k e_t$: $u_1$ is column $t$ of $M_{t+1}$, each
further vector is a matrix-vector product with $M_{t+1}$, and $s_k$ is
the last entry of $u_k$; the last one, $s_{t+1}$, only needs the last
row of $M_{t+1}$ times $u_t$.

Step $t$ costs about $t$ matrix-vector products of dimension $t + 1$, so
the algorithm needs about $r^4/4$ multiply-adds in all.
`_nmod_mat_charpoly_berkowitz` computes each entry of these products as a
dot product with a single reduction
([](algorithms/linear-algebra/dot-products)). Despite its quartic cost
it is the fastest method for $r \le 8$, and `nmod_mat_det` uses it for
$5 \le r \le 8$. Berkowitz designed the algorithm for parallel
computation: the matrix-vector products of all steps can be organised as
products of matrices, giving depth $O(\log^2 r)$.

## Danilevsky's algorithm {#danilevsky}

Danilevsky's method [@Dan1937] transforms $M$ by similarity transforms
into a block triangular matrix whose diagonal blocks are companion
matrices, whose characteristic polynomials can be read off; see also
[@Coh1996] for the related Hessenberg method.
`nmod_mat_charpoly_danilevsky` works on a copy of $M$ from the bottom row
upwards. At step $i$, the rows below $r - i$ already form the last rows of
a companion matrix, and the step makes row $r - i$ equal to the unit
vector $e_{r-i-1}$:

1. Let $h$ be the entry of row $r - i$ just left of the diagonal. If
   $h = 0$, look for a nonzero entry further left in the same row and
   move it into place by exchanging two rows and the corresponding two
   columns (a similarity transform by a permutation).
2. If the row has no nonzero entry left of the diagonal, the trailing
   block of size $i$ is a companion matrix and is decoupled from the rest:
   multiply the result by its characteristic polynomial, read off its
   bottom row, and continue with the leading block of size $r - i$.
3. Otherwise, apply the similarity transform by the identity matrix with
   row $r - i - 1$ replaced by row $r - i$ of $M$, whose inverse involves
   $h^{-1}$. This turns row $r - i$ into $e_{r-i-1}$ and costs one
   rank-one update of the matrix and one row of dot products, $O(r^2)$
   operations.

When the leading block is a full companion matrix, its characteristic
polynomial is multiplied into the result. There are $r - 1$ steps of
cost $O(r^2)$, so $O(r^3)$ operations in all. Every dot product is
computed with a single reduction.

Over a field the algorithm cannot fail. Over $\Z/N\Z$ with composite $N$,
step 3 needs the inverse of $h$; if $h$ is a nonzero zero divisor, the
function returns 0 and `nmod_mat_charpoly` falls back to Berkowitz's
algorithm. A zero divisor can appear even when the result is
well defined: the method depends on the arithmetic of a field, not just
on the ring.

`nmod_mat_charpoly` uses Danilevsky's algorithm for $r > 8$. For prime
$N \le 255$ and $r \ge 32$ it runs it through the `gr` context for
`nmod8`, on a copy of the matrix stored in bytes, which reduces the
memory traffic of the rank-one updates.

## Minimal polynomial {#minpoly}

The minimal polynomial $\mu_M$ of $M$ is the monic generator of the ideal
of polynomials $f$ with $f(M) = 0$. For a vector $v$, the minimal
polynomial $\mu_{M,v}$ of $v$ with respect to $M$ is the monic generator
of the ideal of $f$ with $f(M) v = 0$; its degree is the dimension of the
Krylov space $K(v) = \operatorname{span}\{v, Mv, M^2 v, \ldots\}$, and its
coefficients express the first vector $M^d v$ that is linearly dependent
on the previous ones.

::: theorem {#thm-lcm}
Over a field, if $K(v_1) + \cdots + K(v_s)$ is the whole space, then
$\mu_M = \operatorname{lcm}(\mu_{M, v_1}, \ldots, \mu_{M, v_s})$.
:::

::: proof
Each $\mu_{M, v_j}$ divides $\mu_M$, since $\mu_M(M) v_j = 0$, so their
least common multiple $L$ divides $\mu_M$. Conversely,
$L(M) M^i v_j = M^i L(M) v_j = 0$ for all $i$ and $j$, so $L(M)$
vanishes on a spanning set and $L(M) = 0$. Hence $\mu_M$ divides $L$.
:::

`nmod_mat_minpoly` applies this with unit vectors.

::: algorithm Minimal polynomial {#alg-minpoly}
**Input:** $M \in (\Z/p\Z)^{r \times r}$.

**Output:** $\mu_M$.

1. $\mu \gets 1$, $W \gets \{0\}$, $j \gets 0$.
2. While $W$ is not the whole space:
   1. compute $v = e_j, Mv, M^2 v, \ldots$, reducing each new vector
      against the previous ones (kept in echelon form, together with the
      combination that produced them), until one reduces to zero; the
      combination gives $\mu_{M, e_j}$;
   2. $\mu \gets \operatorname{lcm}(\mu, \mu_{M, e_j})$, computed as
      $\mu \mu_{M,e_j} / \gcd(\mu, \mu_{M,e_j})$;
   3. add the Krylov vectors to $W$, also kept in echelon form, and let
      $j$ be the first index such that $e_j$ is not in $W$.
3. Return $\mu$.
:::

The rows are reduced with `nmod_mat_reduce_row`, which accumulates the
updates of a row in up to three words and reduces each entry once. Each
Krylov vector costs a matrix-vector product and a row reduction,
$O(r^2)$ operations. If $e_0$ already generates the space, which is
often the case, there are $r + 1$ vectors and the cost is $O(r^3)$. A
matrix whose minimal polynomial has small degree, such as a scalar
matrix, needs up to $r$ generators, and the number of Krylov vectors is
bounded by $r (r + 1)$ in all cases. Algorithms with $O(r^\omega)$
operations, also based on Krylov spaces, are described in [@Sto2000].
