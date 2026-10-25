---
title: Optimal precision in p-adic linear algebra
summary: The precision of a result is bounded by its differential; how FLINT attains this bound for determinants and linear systems by full pivoting on valuations (blocked by valuation levels), and for characteristic polynomials by computing on an exact lift and assigning the precision of the differential, compared with division-free and partially pivoted algorithms.
complexity: Determinant and linear solving: the cost of a blocked LU decomposition with fast matrix multiplication, plus one rank-revealing LU decomposition over $\mathbf{F}_p$ per valuation level. Characteristic polynomial: $O(n^4)$ (Berkowitz) for matrices that are cyclic modulo $p$, otherwise an $O(n^3)$ Hessenberg reduction at doubled precision and $n$ matrix products modulo $p^N$.
order: 14
---

In the ball model ([](algorithms/p-adic/precision#balls)) every result
comes with a guaranteed precision, but different algorithms for the same
problem can guarantee very different precisions: the error terms
propagated operation by operation may be much larger than the actual
uncertainty of the result. Caruso, Roe and Vaccon [@CRV2014; @CRV2015;
@CRV2017] showed how to determine the actual uncertainty from the
differential of the function computed, and studied which algorithms
attain it; [@Car2017] surveys these results. This page describes the
methods `padic_radix` uses to attain the optimal precision for
determinants, linear systems and characteristic polynomials.

## Differential precision {#differential}

::: lemma First- and second-order terms {#lem-differential}
Let $f : \mathbf{Z}_p^m \to \mathbf{Z}_p$ be a polynomial with
coefficients in $\mathbf{Z}_p$, let $x \in \mathbf{Z}_p^m$, and let
$\varepsilon \in \mathbf{Q}_p^m$ with $\operatorname{ord}_p(\varepsilon_i)
\ge N_i \ge 0$. Then
$$
f(x + \varepsilon) - f(x) = \sum_i \frac{\partial f}{\partial x_i}(x)\,
\varepsilon_i + O(p^{2 N_{\min}}), \qquad N_{\min} = \min_i N_i.
$$
Hence $f(x + \varepsilon) \in f(x) + O(p^{M})$ with
$M = \min(2N_{\min}, \min_i (N_i + \operatorname{ord}_p \partial_i f(x)))$,
and if $N_i + \operatorname{ord}_p \partial_i f(x) < 2 N_i$ for the index
attaining the minimum, no better precision holds for all such
$\varepsilon$.
:::

::: proof
Expand $f(x + \varepsilon)$ by Taylor's formula: the coefficients of
the terms of degree $d \ge 2$ in $\varepsilon$ are values at $x$ of
divided derivatives $\partial^\alpha f / \alpha!$, which are polynomials
with coefficients in $\mathbf{Z}_p$, so these terms have valuation at
least $d N_{\min} \ge 2 N_{\min}$. For the last claim take
$\varepsilon = p^{N_i} e_i$: the change is $p^{N_i}\partial_i f(x)$ plus
terms of valuation at least $2N_i$.
:::

For a matrix function the partial derivatives are entries of a matrix:
$d \det(A) = \operatorname{tr}(\operatorname{adj}(A)\, dA)$, and
$d \chi_A = -\operatorname{tr}(\operatorname{adj}(xI - A)\, dA)$ for the
characteristic polynomial $\chi_A(x) = \det(xI - A)$. Non-integral inputs
are reduced to integral ones by scaling with a power of $p$. The
quantity $\min_i (N_i + \operatorname{ord}_p \partial_i f)$ is the
*optimal precision* below.

## Determinants and full pivoting {#det}

::: theorem Full pivoting attains the optimal precision {#thm-full-pivoting}
Let $A \in M_n(\mathbf{Z}_p)$ be known to $O(p^N)$ entrywise, with
elementary divisors $p^{a_1}, \ldots, p^{a_n}$, $a_1 \le \cdots \le a_n
< N$. Gaussian elimination in ball arithmetic in which each pivot has
minimal valuation in the remaining submatrix (full pivoting on
valuations) produces pivots of valuations $a_1, \ldots, a_n$ and
returns $\det A$ to precision $N + a_1 + \cdots + a_{n-1}$, which is the
optimal precision.
:::

::: proof
Let $d$ be the pivot of the current submatrix $S$, of minimal valuation
$a$, and let $y$ and $z$ be entries of its row and column. An entry $x$
becomes $x - yz/d$, and all four inputs are known to absolute precision
$N$. By [](algorithms/p-adic/precision#prop-ring) the error terms of the
update come from $x$ (precision $N$), from $y$ and $z$
($N + \operatorname{ord}_p(z) - a \ge N$ and likewise), and from $d$
($N + \operatorname{ord}_p(y) + \operatorname{ord}_p(z) - 2a \ge N$), so
every Schur complement is again known to precision $N$. Elimination is
a sequence of row and column operations in $\mathrm{GL}_n(\mathbf{Z}_p)$
followed by the removal of a pivot, so the elementary divisors of the
Schur complement are $p^{a_2}, \ldots, p^{a_n}$; since the minimal
valuation of the entries of a matrix is its smallest elementary divisor,
and $a_n < N$ makes the valuations visible at precision $N$, the pivots
have the valuations $a_k$. The determinant is the product of the pivots
up to sign; a product of balls of absolute precision $N$ and valuations
$a_k$ has precision $\min_k (N + \sum_{l \ne k} a_l) = N + a_1 + \cdots +
a_{n-1}$. Finally, with $A = P D Q$ for $P, Q \in
\mathrm{GL}_n(\mathbf{Z}_p)$ and $D$ diagonal,
$\operatorname{adj}(A) = \operatorname{adj}(Q) \operatorname{adj}(D)
\operatorname{adj}(P)$, whose minimal valuation is that of
$\operatorname{adj}(D)$, namely $a_1 + \cdots + a_{n-1}$. This is the
optimal precision of [](#lem-differential), attained since $a_n < N$.
:::

Other algorithms lose precision. A division-free algorithm (Berkowitz's
algorithm, cofactor expansion) computes only sums of products of
entries, each known to $O(p^N)$, so it can guarantee little more than
precision $N$, losing up to $\operatorname{ord}_p(\det A) - a_n$ digits.
Elimination that takes the first nonzero entry of a column as pivot
loses digits whenever that entry is not of minimal valuation, which for
small $p$ happens often even for random matrices. Choosing the entry of
minimal valuation in the pivot column (partial pivoting) is optimal for
matrices that are invertible modulo $p$, where all pivots are units, and
nearly optimal in general. Average digits lost against the optimal
precision for $p = 7$, entries known to $O(7^{20})$, and matrices
$U \operatorname{diag}(7^{a_i}) V$ with random $a_i \le 6$:

| $n$ | Berkowitz | first nonzero pivot | column pivoting | full pivoting |
|---|---|---|---|---|
| 12 | 28 | 2.7 | 0.6 | 0 |
| 40 | 115 | 10.3 | 1.2 | 0 |
| 60 | 173 | 16 | 0.7 | 0 |

### Blocked full pivoting {#levels}

Full pivoting needs the updated Schur complement before each pivot is
chosen, which rules out the blocked algorithms that make LU
decomposition fast. However, only valuations matter for the choice of
pivots, and these can be read off modulo $p$:

::: algorithm Elimination by valuation levels {#alg-levels}
**Input:** a square matrix $S$ over $\mathbf{Q}_p$.

**Output:** $\det S$.

1. If $S$ has no nonzero entry, return the determinant of the block of
   zeros (bounded by the precisions). Otherwise let $a$ be the minimal
   valuation of the entries of $S$.
2. Compute a rank-revealing LU decomposition of $T = (S/p^a) \bmod p$
   over $\mathbf{F}_p$; it selects $r$ rows and columns such that the
   corresponding submatrix of $T$ is invertible. Permute them to the
   front: $S = \begin{pmatrix} S_{11} & S_{12} \\ S_{21} & S_{22}
   \end{pmatrix}$ with $S_{11} = p^a U$, $U \in \mathrm{GL}_r(\mathbf{Z}_p)$.
3. Compute $\det S_{11}$ from an LU decomposition of $S_{11}$ (all of
   whose pivots have valuation $a$), and the Schur complement $S' =
   S_{22} - S_{21} (S_{11}^{-1} S_{12})$ by triangular solves and one
   matrix product.
4. Return $\pm \det(S_{11}) \det(S')$, with $\det S'$ computed
   recursively.
:::

The entries of $S'$ have valuation at least $a + 1$, because $T$ has
rank $r$ and $S'/p^a \equiv T_{22} - T_{21} T_{11}^{-1} T_{12} \equiv 0
\pmod p$. The $r$ pivots of a level have the same valuation $a$, which
is minimal in every intermediate Schur complement, so the computation is
a full pivoting elimination carried out in blocks, with the same
precision. A matrix $A$ such that $A/p^a$ is invertible modulo $p$ for
the minimal valuation $a$ of its entries (in particular a matrix
invertible modulo $p$) needs a single level, that is, one blocked LU
decomposition. The cost is that of the
blocked LU decomposition, which uses the fast matrix products of
[](algorithms/p-adic/block-products#matrices), plus $O(n^3)$ word
operations modulo $p$ per level.

## Linear systems {#solve}

For $AX = B$ the differential is $dX = A^{-1}(dB - dA\, X)$, so for
entries of $A$ and $B$ known to $O(p^{N_A})$ and $O(p^{N_B})$ the
optimal precision of $X_{ij}$ is
$$
\min_l \Big(\operatorname{ord}_p (A^{-1})_{il} + \min\big(N_B, \; N_A +
\min_m \operatorname{ord}_p X_{mj}\big)\Big).
$$
Gauss–Jordan elimination with full pivoting on valuations attains it
in experiments, while LU decomposition with column pivoting followed by
back substitution loses digits when $A$ is not invertible modulo $p$
(about 0.3 to 3.6 digits per entry of $X$ for $p = 2$ and $n = 5$ to
20). The difference lies in the order of operations: back substitution
computes $X_1 = S_{11}^{-1}(B_1 - S_{12} X_2)$ from the already rounded
$X_2$, whose errors are correlated with those of $S_{12}$ in a way the
ball arithmetic cannot see, whereas Gauss–Jordan elimination eliminates
the later pivots from the earlier pivot rows before any division by
them.

`padic_radix_mat_nonsingular_solve` therefore runs Gauss–Jordan
elimination by valuation levels on $[A \mid B]$: at each level the
pivot block of [](#alg-levels) is applied to its rows (LU decomposition
and triangular solves), and its columns are eliminated from all the
other rows, including the pivot rows of the previous levels, by one
matrix product. When $A/p^a$ is invertible modulo $p$ for the minimal
valuation $a$, all pivots have valuation $a$ and the LU solve with
column pivoting is already optimal, so it is used directly. The cost is
that of a blocked LU solve. `gr_mat_nonsingular_solve` and `gr_mat_inv`
call these functions for `padic_radix` matrices.

## Characteristic polynomials {#charpoly}

Write $\operatorname{adj}(xI - A) = \sum_{k=0}^{n-1} B_k x^k$, so that
$B_{n-1} = I$ and $B_{k-1} = A B_k + c_k I$ for $\chi_A = \sum c_k x^k$.
By [](#lem-differential), for an integral matrix with entries known to
$O(p^{N_{ij}})$ the coefficient $c_k$ is known to precision
$$
N_k = \min\big(2 N_{\min}, \; \min_{ij} (N_{ij} + \operatorname{ord}_p
(B_k)_{ji})\big),
$$
and to no better precision when the minimum is attained by the
first-order term. The higher-order terms are bounded by $2N_{\min}$
because their coefficients are minors of $xI - A$, which are integral.

::: proposition Cyclic matrices {#prop-cyclic}
If $A \in M_n(\mathbf{Z}_p)$ is cyclic modulo $p$ (its minimal
polynomial modulo $p$ has degree $n$) and all entries are known to
$O(p^N)$, then the optimal precision of every coefficient of $\chi_A$
is $N$, and any division-free algorithm, such as Berkowitz's, attains
it.
:::

::: proof
$B_k = q_k(A)$ for the monic polynomial $q_k = \sum_{i > k} c_i x^{i -
k - 1}$ of degree $n - k - 1 < n$. Since $I, A, \ldots, A^{n-1}$ are
linearly independent modulo $p$, $q_k(A) \not\equiv 0 \pmod p$, so some
entry of $B_k$ is a unit and $N_k = N$. A division-free algorithm
computes the coefficients from integral entries known to $O(p^N)$ by
additions and multiplications, which preserve absolute precision $N$
([](algorithms/p-adic/precision#prop-ring)).
:::

Matrices that are not cyclic modulo $p$ are where precision is lost:
for $U \operatorname{diag}(p^{a_i}) V$ with $a_i \le 5$ and $n = 16$,
Berkowitz's algorithm loses 35 digits on the constant coefficient and
11 digits on average, and the division-based algorithms (Hessenberg
reduction, Danilevsky's algorithm, see
[](algorithms/linear-algebra/charpoly)) lose a few digits on every
coefficient, even for generic matrices. `padic_radix_mat_charpoly`
therefore uses a division-free algorithm only for matrices that are
integral, of uniform precision $N$ and cyclic modulo $p$ (tested by the
rank of a random Krylov basis over $\mathbf{F}_p$). Even there it first
runs the Hessenberg reduction, which costs $O(n^3)$ operations instead of
$O(n^4)$ and loses nothing on matrices that are generic modulo $p$: when
all its coefficients come out with precision $N$, they are optimal by
[](#prop-cyclic), and only otherwise is Berkowitz's algorithm run (at
$10^4$ digits and $n = 30$ this is about 12 times faster for random
matrices). Other matrices get the optimal precision by computing on an
exact lift:

::: algorithm Characteristic polynomial with optimal precision {#alg-charpoly-lift}
**Input:** $A$ with entries known to $O(p^{N_{ij}})$.

**Output:** $\chi_A$, each coefficient to precision $N_k$.

1. Scale $A$ by a power of $p$ to make all entries and error terms
   integral, and let $M \in M_n(\Z)$ be the exact lift of the entries.
2. Compute $\chi_M \bmod p^{2N_{\min}}$: run the Hessenberg reduction with
   pivots of minimal valuation (similarity transformations in
   $\mathrm{GL}_n(\mathbf{Z}_p)$, $O(n^3)$ operations) in ball arithmetic
   on $M$ known to $O(p^{W})$ for some $W > 2N_{\min}$; if a coefficient
   comes out with precision below $2N_{\min}$, increase $W$, and as a
   last resort use Berkowitz's algorithm, which loses nothing.
3. Compute $B_{n-1}, \ldots, B_0$ by the recurrence modulo one word
   (the largest power of $p$ below $2^{64}$), which gives
   $\operatorname{ord}_p (B_k)_{ji}$ when it is smaller; repeat modulo
   $p^{N_{\min}}$ only if a valuation at that bound determines some
   $N_k$. These are $n$ matrix products over $\Z/p^K\Z$ (`nmod_mat` or
   `mpn_mod_mat`).
4. Return the coefficients of $\chi_M$ with the precisions $N_k$,
   scaled back.
:::

Where the bound $2N_{\min}$ for the higher-order terms limits the
result, Danilevsky's algorithm in ball arithmetic is also run, and the
more precise of the two valid results is kept for each coefficient. In
the experiments this gives the optimal precision up to a fraction of a
digit on average for $p = 2, 7$ and $1009$. An exact matrix gives the
exact characteristic polynomial (computed by `fmpz_mat_charpoly`),
rounded to the context precision. The cost is about that of Berkowitz's
algorithm for cyclic matrices and 1.5 to 5 times more for structured
ones.

The lift works with $K = N_{\min} + s$ digits, where $p^s$ is the
scaling, and with entries of $M$ of any size. Entries of high valuation
are harmless, since everything is computed modulo at most $p^{4K + 64}$
and $M$ can be reduced modulo it. A large $K$ is not: one entry of
valuation $-10^6$ among entries known to 20 digits gives $K \approx 10^6$,
and the lift then works at a million digits for results that carry a few
dozen. The lift gives up when $K$ exceeds $64(R + n) + 1024$ for the
largest relative precision $R$ of an inexact entry (and for exact
matrices whose entries are dominated by powers of $p$ or exceed what the
context keeps), and Berkowitz's algorithm in ball arithmetic, whose cost
follows the precisions of the entries, is used instead.

## Use in FLINT {#use}

- `padic_radix_mat_find_nonzero_pivot` chooses pivots of minimal
  valuation in a column for the generic algorithms of `gr_mat`.
- `padic_radix_mat_det_lu_full` and `padic_radix_mat_det_levels`
  implement [](#thm-full-pivoting) and [](#alg-levels);
  `padic_radix_mat_det` (the `gr_mat_det` method) uses the first for
  small matrices and the second above.
- `padic_radix_mat_nonsingular_solve` and `padic_radix_mat_inv`
  implement [](#solve) and are the `gr_mat_nonsingular_solve` and
  `gr_mat_inv` methods.
- `padic_radix_mat_charpoly_lift`, `padic_radix_mat_charpoly_hessenberg`
  and `_padic_radix_mat_charpoly_berkowitz` implement the algorithms of
  [](#charpoly), and `padic_radix_mat_charpoly` (the `gr_mat_charpoly`
  method) chooses between them.
- The context flag `PADIC_RADIX_FAST_LINALG` selects the fastest
  algorithms instead of the optimal ones (the Hessenberg reduction for
  characteristic polynomials, LU decompositions for linear systems), for
  computations at high working precision where a few lost digits are
  cheaper than the optimal algorithms.
- The tutorial [](tutorials/numerics/padic-interpolation) shows the
  optimal precision of a linear system at work.
