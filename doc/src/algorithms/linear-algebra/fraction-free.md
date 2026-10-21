---
title: Fraction-free elimination
summary: Bareiss' fraction-free Gaussian elimination over integral domains, Sylvester's identity, the fraction-free LU decomposition and how FLINT uses it for determinants, ranks, linear systems and echelon forms over $\Z$, $\Z[x]$ and $(\Z/p\Z)[x]$.
complexity: $O(mnr)$ ring operations for an $m \times n$ matrix of rank $r$; over $\Z$ with $b$-bit entries the operands have $O(n(b + \log n))$ bits, giving $O(n^3 M(n(b + \log n)))$ word operations for an $n \times n$ matrix.
order: 45
---

Let $R$ be an integral domain, such as $\Z$, $\Z[x]$ or $\mathbb{F}_p[x]$,
and $A$ an $m \times n$ matrix over $R$. Gaussian elimination over the
fraction field creates quotients whose numerators and denominators grow;
removing common factors needs gcd computations, which are expensive or
unavailable. Fraction-free elimination keeps all entries in $R$, uses
only exact divisions, and bounds every intermediate entry by a minor of
$A$ [@Bar1968; @NakTurWil1997]. It is the method of choice for small
matrices over $\Z$ and for matrices over polynomial rings, where modular
methods are less developed.

## Generic version {#generic}

`gr_mat_fflu` implements the algorithm below over any ring with exact
division, and `gr_mat_det_fflu`, `gr_mat_rank_fflu`,
`gr_mat_nonsingular_solve_fflu`, `gr_mat_nonsingular_solve_den_fflu` and
`gr_mat_rref_den_fflu` build on it. The specialised versions are separate
implementations with the same structure:

- `fmpz_mat_fflu` over $\Z$ chooses the first nonzero pivot, and while
  the entries fit in a word computes the updates in two-word arithmetic
  with a precomputed inverse of the previous pivot;
- `fmpz_poly_mat_fflu` and `nmod_poly_mat_fflu` over $\Z[x]$ and
  $(\Z/n\Z)[x]$ choose, among the nonzero entries of the column, one of
  least degree and, over $\Z[x]$, with the smallest coefficients
  (`fmpz_poly_mat_find_pivot_partial`), which keeps the entries small;
- `ca_mat_fflu` over exact real and complex numbers chooses pivots with
  zero tests that may fail
  ([linear algebra over exact numbers](algorithms/exact/linear-algebra#fflu)).

## Bareiss' algorithm {#bareiss}

::: algorithm Fraction-free elimination {#alg-bareiss}
**Input:** $A \in R^{m \times n}$.

**Output:** the rank $r$, a row permutation, and the matrix $B$ described
below; `den` $= p_{r-1}$.

1. $B \gets A$, $p_{-1} \gets 1$, $k \gets 0$ (pivot row), $c \gets 0$
   (column).
2. While $k < m$ and $c < n$:
   1. find a row $i \ge k$ with $B_{ic} \ne 0$; if there is none,
      $c \gets c + 1$ and continue; otherwise swap rows $i$ and $k$;
   2. $p_k \gets B_{kc}$, $c_k \gets c$;
   3. for all $i > k$ and $j > c$:
      $B_{ij} \gets (p_k B_{ij} - B_{ic} B_{kj}) / p_{k-1}$;
   4. $k \gets k + 1$, $c \gets c + 1$.
3. Return $r = k$.
:::

The divisions in step 2.3 are exact. This follows from Sylvester's
identity.

::: theorem {#thm-sylvester}
Assume that no row exchanges are needed and that the pivot columns are
$0, 1, \ldots, r - 1$. After the pivots $p_0, \ldots, p_{k-1}$ have been
processed, for $i, j \ge k$,
$$
B_{ij} = \det A[\{0, \ldots, k - 1, i\}, \{0, \ldots, k - 1, j\}],
$$
the minor of $A$ formed by the leading $k$ rows and columns bordered by
row $i$ and column $j$. In particular $p_k$ is the leading principal
$(k+1) \times (k+1)$ minor of $A$.
:::

::: proof
By induction on $k$; for $k = 0$ the claim is $B_{ij} = A_{ij}$. Suppose
it holds before pivot $k$ is processed, and let $T_0$ be the leading
$k \times k$ submatrix of $A$, with $\det T_0 = p_{k-1}$ (1 if $k = 0$),
and $T = A[\{0, \ldots, k, i\}, \{0, \ldots, k, j\}]$. The four
$(k+1) \times (k+1)$ minors of $T$ that contain $T_0$ are, by induction,
$p_k$, $B_{kj}$, $B_{ik}$ and $B_{ij}$, and Sylvester's determinant
identity states
$$
\det(T) \det(T_0) = p_k B_{ij} - B_{ik} B_{kj}.
$$
Hence the update of step 2.3 produces $\det T$, which is the claim for
$k + 1$.
:::

With row exchanges, the same holds for the row-permuted matrix, and with
skipped columns for the matrix restricted to the pivot columns and the
column $j$. Consequently every intermediate entry is a minor of $A$: over
$\Z$ it is bounded by Hadamard's inequality
([](algorithms/linear-algebra/multimodular#bounds)), with about
$n(b + \frac12 \log_2 n)$ bits for $b$-bit entries, and over $K[x]$ its
degree is at most $n$ times the largest degree. Gaussian elimination over
$\Q$ without gcds would instead produce entries whose sizes grow
exponentially with $k$.

For a nonsingular $n \times n$ matrix, the last pivot is
$p_{n-1} = \det(PA) = \operatorname{sgn}(P) \det A$, so
`fmpz_mat_det_bareiss` returns `den` with the sign of the permutation.
Each of the $O(mnr)$ updates costs two multiplications and one exact
division of minors, so over $\Z$ the cost is
$O(n^3 M(n(b + \log n)))$ word operations for an $n \times n$ matrix,
that is $O(n^5 b^2)$ with classical integer multiplication. This is
competitive with the modular methods of
[](algorithms/linear-algebra/multimodular) only for small $n$ (FLINT's
crossovers are $n \approx 15$ to 25).

## The fraction-free LU decomposition {#fflu}

FLINT stores the result in place, like an LU decomposition: the upper
part of $B$, from the pivot of each row on, holds the fraction-free
echelon form $U$, and below each pivot $p_k$ the entries of column $c_k$
are kept as they were when $p_k$ was chosen (step 2.3 does not change
column $c$). Let $L$ be the lower triangular matrix formed by these
columns, with the pivots $p_k$ on its diagonal.

::: proposition {#prop-fflu}
For a nonsingular $n \times n$ matrix $A$ with row permutation $P$,
$$
PA = L \, D^{-1} \, U, \qquad
D = \operatorname{diag}(p_{-1} p_0, \, p_0 p_1, \, \ldots, \, p_{n-2} p_{n-1}).
$$
:::

This is the fraction-free LU decomposition of [@NakTurWil1997; @Jef2010].
For the $4 \times 4$ example of [](fmpz_mat/rref#example), with pivots
$2, 5, 18, 82$, one checks $A = L \operatorname{diag}(2, 10, 90, 1476)^{-1} U$.

## Solving {#solve}

To solve $AX = B$ for nonsingular $A$, `fmpz_mat_solve_fflu_precomp`
applies to $PB$ the same fraction-free updates as step 2.3 (forward
substitution: the column $j$ of the right-hand side is updated like a
column of $A$), which yields the last column of the elimination of the
augmented matrix $(A \mid B)$, and then performs back substitution
$$
x_i = \frac{p_{n-1} y_i - \sum_{j > i} U_{ij} x_j}{U_{ii}},
\qquad i = n - 1, \ldots, 0,
$$
where $y$ is the transformed right-hand side. The result is
$X = p_{n-1} A^{-1} B$, which is integral by Cramer's rule since
$p_{n-1} = \pm \det A$, so all divisions are exact. `fmpz_mat_solve_fflu`
then fixes the sign, returning $\mathit{den} = \det A$. For a
rank-deficient $A$ the same procedure, with the non-pivot unknowns set to
zero, either produces a candidate or meets an inexact division and
reports failure; `fmpz_mat_can_solve_fflu` then checks the remaining
equations exactly.

## Reduced echelon form {#rref}

Let the first $r$ rows of the echelon form be $(U \mid V)$ up to a
permutation of the columns, where $U$ is the $r \times r$ upper
triangular matrix of pivot columns, with the pivots $p_0, \ldots,
p_{r-1}$ on its diagonal. The reduced echelon form over the fraction
field is $(I \mid U^{-1} V)$ in the original column order. Since
$(U \mid V)$ is obtained from the pivot rows $(A_P \mid A_N)$ of $A$ by
invertible row operations, $U^{-1} V = A_P^{-1} A_N$, and
$p_{r-1} U^{-1} V = \pm \det(A_P) A_P^{-1} A_N$ is integral by Cramer's
rule. `fmpz_mat_rref_fflu` and `fmpz_poly_mat_rref` compute it by back
substitution with exact divisions by the pivots, one non-pivot column at
a time, and set the pivot columns to $p_{r-1}$ times unit vectors. This
is fraction-free Gauss–Jordan elimination [@NakTurWil1997] restricted to
the non-pivot columns, which are all that change.

## Use in FLINT {#use}

| Ring | Functions |
|---|---|
| any (`gr`) | `gr_mat_fflu`, `gr_mat_det_fflu`, `gr_mat_rank_fflu`, `gr_mat_nonsingular_solve_fflu`, `gr_mat_rref_den_fflu` |
| $\Z$ | `fmpz_mat_fflu`, `fmpz_mat_det_bareiss`, `fmpz_mat_solve_fflu`, `fmpz_mat_solve_fflu_precomp`, `fmpz_mat_can_solve_fflu`, `fmpz_mat_rref_fflu`, `fmpz_mat_rank` (small matrices), `fmpz_mat_lll_storjohann` (on the Gram matrix) |
| $\Q$ | `fmpq_mat_solve_fraction_free`, `fmpq_mat_rref_fraction_free` (after clearing denominators) |
| $\Z[x]$ | `fmpz_poly_mat_fflu`, `fmpz_poly_mat_det_fflu`, `fmpz_poly_mat_solve_fflu`, `fmpz_poly_mat_rref`, `fmpz_poly_mat_inv`, `fmpz_poly_mat_rank` |
| $(\Z/n\Z)[x]$ | `nmod_poly_mat_fflu`, `nmod_poly_mat_det_fflu`, `nmod_poly_mat_solve_fflu`, `nmod_poly_mat_rref`, `nmod_poly_mat_inv`, `nmod_poly_mat_rank` |
