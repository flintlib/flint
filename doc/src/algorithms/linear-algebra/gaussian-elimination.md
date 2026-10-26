---
title: Gaussian elimination and LU decomposition
summary: LU decomposition with row pivoting over $\Z/p\Z$ in its classical, delayed-reduction and recursive forms, and how FLINT derives triangular solving, linear system solving, inverses, determinants, ranks, echelon forms and nullspaces from it; plus division-free determinants and the Howell form over $\Z/n\Z$.
complexity: $O(m n r^{\omega - 2})$ field operations for an $m \times n$ matrix of rank $r$, where $\omega \le 3$ is the exponent of the matrix multiplication used; $O(r^\omega)$ for solving, inversion and determinants of $r \times r$ matrices.
order: 40
---

Throughout, $p$ is a prime and $A$ an $m \times n$ matrix over the field
$\Z/p\Z$; for the last sections $N$ is an arbitrary modulus. Gaussian
elimination brings $A$ to row echelon form by row operations, and
recording the operations gives an LU decomposition, from which FLINT
computes solutions of linear systems, inverses, determinants, ranks,
reduced echelon forms and nullspaces. The recursive variant reduces all
of this to matrix multiplication
[@BH1974; @IMH1982; @JPS2013], so that the fast methods of
[](algorithms/linear-algebra/matrix-multiplication) apply.

## Generic version {#generic}

`gr_mat_lu` computes the decomposition $PA = LU$ over any field, and
over rings pretending to be fields, with the same two algorithms as
below: the classical elimination (`gr_mat_lu_classical`) and the
recursive one (`gr_mat_lu_recursive`), which splits the columns,
recurses on the left half, updates the right half with
`gr_mat_nonsingular_solve_tril` and `gr_mat_mul`, and recurses on the
lower right block. The default uses the classical algorithm when the
matrix has fewer than 5 rows or columns. Pivots are chosen by a ring
method: the first provably nonzero entry by default (`GR_UNABLE` if
none is provably nonzero but some cannot be decided, as for balls
containing zero), and the entry of largest magnitude in floating-point
rings. Over $\Z/n\Z$ the generic decomposition is not
`nmod_mat_lu`: it does without the delayed reductions of the classical
kernel, but its recursive form uses `nmod_mat_mul` and therefore runs
at the speed of matrix multiplication for large matrices. `fmpz_mod`,
`mpn_mod` and `nfloat` override the decomposition. Over integral
domains that are not fields, the generic functions use fraction-free
elimination instead ([](algorithms/linear-algebra/fraction-free)); the
choice for each operation is tabulated in
[](algorithms/linear-algebra/overview#lu).

## The decomposition {#ple}

::: definition {#def-ple}
A *PLE decomposition* of $A$ is a factorisation $A = \Pi L E$ where
$\Pi$ is an $m \times m$ permutation matrix, $L$ is an $m \times r$
lower triangular matrix with ones on the diagonal, and $E$ is an
$r \times n$ matrix in row echelon form whose $r$ rows are nonzero. The
pivot of row $i$ of $E$ is its first nonzero entry, in column $c_i$,
with $c_0 < c_1 < \cdots < c_{r-1}$.
:::

Then $r$ is the rank of $A$, since $L$ has full column rank and $E$ full
row rank. FLINT stores $L$ below the diagonal and $E$ on and above it, in
the array of $A$; the permutation is returned as an array `P` with
$(LE)_i = A_{P[i]}$, see [](nmod_mat/lu#output). Elimination chooses as
pivot of column $c$ the first row, among the rows not yet used, whose
entry in column $c$ is nonzero after the previous eliminations.

::: proposition {#prop-profile}
With this pivoting rule, column $c$ receives a pivot if and only if it is
not a linear combination of the columns to its left. The pivot columns
$c_0 < \cdots < c_{r-1}$ are therefore the column rank profile of $A$,
the lexicographically smallest set of $r$ linearly independent columns.
:::

::: proof
After the eliminations for the pivots in columns before $c$, the rows
not yet used have zeros in all those columns, and their entries in
column $c$ are zero exactly when column $c$ of $A$ lies in the span of
the earlier pivot columns, which span all earlier columns.
:::

The row permutation, in contrast, is not the row rank profile in
general when the pivot rows are brought into place by swaps, as in
`nmod_mat_lu`. Using upward rotations of rows instead reveals both
profiles, and more [@JPS2013; @DPS2017]; see [](#rank-profile). Here a
*swap* exchanges two rows and leaves the others in place, while an
*upward rotation* moves a row up to a given position and the rows it
passes over down by one position each, keeping their order; a *downward
rotation* is the inverse move.

## Classical elimination {#classical}

::: algorithm Gaussian elimination {#alg-classical}
**Input:** $A \in (\Z/p\Z)^{m \times n}$.

**Output:** the rank $r$, the permutation and $L$, $E$ stored in $A$.

1. $r \gets 0$. For $c = 0, 1, \ldots, n - 1$ while $r < m$:
   1. find the first row $i \ge r$ with $A_{ic} \ne 0$; if there is none,
      continue with the next $c$;
   2. swap rows $i$ and $r$ and record the swap;
   3. for each row $i > r$: $\ell \gets A_{ic} A_{rc}^{-1}$, subtract
      $\ell$ times row $r$ from row $i$ in columns $c + 1, \ldots, n - 1$
      (nothing to do if $\ell = 0$), and store $\ell$ in column $r$ of
      row $i$;
   4. $r \gets r + 1$.
:::

For a nonsingular $d \times d$ matrix this costs
$\sum_{i < d} (d - 1 - i)^2 \approx d^3/3$ multiply-adds, $d$ inversions and,
in `nmod_mat_lu_classical`, one reduction modulo $p$ per multiply-add.

## Rank profile matrix {#rank-profile}

::: definition {#def-rpm}
The *rank profile matrix* of $A$, of rank $r$, is the $m \times n$ matrix
$\mathcal{R}_A$ with $r$ entries equal to 1, no two in the same row or
column, and zeros elsewhere, such that for all $i \le m$ and $j \le n$ the
submatrices of $A$ and $\mathcal{R}_A$ formed by their first $i$ rows and
first $j$ columns have the same rank.
:::

It exists and is unique [@DPS2017]: writing $\rho(i, j)$ for the rank of
the leading $i \times j$ submatrix of $A$, the entry $(i, j)$ of
$\mathcal{R}_A$ is
$\rho(i+1, j+1) - \rho(i, j+1) - \rho(i+1, j) + \rho(i, j)$. Its rows
with a nonzero entry form the row rank profile of $A$, its columns with a
nonzero entry the column rank profile, and the same holds for every
leading submatrix of $A$ with the entries of $\mathcal{R}_A$ that lie in
it.

::: theorem {#thm-rotations}
Run [](#alg-classical) with step 2 replaced by the upward rotation of
rows $r, \ldots, i$: move row $i$ to position $r$ and the rows in
positions $r, \ldots, i - 1$ down by one position.
Let `P[k]` be the index in $A$ of the row in position $k$ at the end.
Then `P[r]` $< \cdots <$ `P[m-1]`, and the nonzero entries of
$\mathcal{R}_A$ are at the positions (`P[k]`, $c_k$) for $k < r$.
:::

::: proof
An upward rotation keeps the relative order of the rows it passes over,
hence preserves the following property: the rows that are not yet pivot rows are in
increasing order of their indices in $A$. This gives the first claim, and
it shows that the pivot of column $c$ is in the row of smallest index, in
$A$, among the rows that are not yet pivot rows and have a nonzero entry
in column $c$. The rows of smaller index have a zero multiplier, so every
row of $A$ is only modified by subtracting multiples of rows of smaller
index: for every $i$, the operations transform the first $i$ rows of $A$
by an invertible lower triangular matrix, which does not change the rank
of the first $j$ columns of these rows. After the steps for the columns
before $j$, these columns are zero in the rows that are not pivot rows,
and in the first $j$ columns the pivot rows form a matrix in row echelon
form, with nonzero rows the pivot rows whose pivot column is before $j$.
The rank of the leading $i \times j$ submatrix of $A$ is therefore the
number of pairs (`P[k]`, $c_k$) with `P[k]` $< i$ and $c_k < j$, which is
the definition of $\mathcal{R}_A$.
:::

`nmod_mat_ple` uses this rule. By the proof, the pivot of a column is
the row of smallest index in $A$ with a nonzero entry among the rows that
are not pivot rows yet, so that the order of the rows need not be kept
in memory: `nmod_mat_ple_classical` keeps the rows in place and records
their order in an array, an upward rotation costing $O(i - r)$
operations on indices, and `nmod_mat_ple_classical_delayed` swaps the
pivot row with the row in position $r$; both put the rows in their order
at the end, and `nmod_mat_ple_recursive` does so once for each half at
each level of the recursion. With swaps, as in
`nmod_mat_lu`, the pivot columns are the same, but a row swapped
downwards passes over rows of larger index, and the pivot rows need not
be the row rank profile. The output of `nmod_mat_ple` is unique, since
$\mathcal{R}_A$ determines the pivot rows and columns, and the other rows
are in increasing order: all its algorithms return the same
decomposition. Pivoting rules that reveal $\mathcal{R}_A$ in other
search orders, and the corresponding recursive algorithms, are studied in
[@DPS2017].

## Delayed reduction {#delayed}

The reductions dominate the cost of step 3 for small dimensions.
`nmod_mat_lu_classical_delayed` performs the same operations but lets the
entries of the active submatrix grow: a nonzero multiplier is turned
into $f = p - \ell \in [1, p - 1]$, and $f$ times the pivot row is added
without reduction, in one, two or three words. An entry is reduced
only when it is about to be used: when its column becomes the pivot
column, before the pivot is chosen, and when its row becomes the pivot
row.

::: proposition {#prop-delayed}
Let $s = \min(m, n)$. In [](#alg-classical) with delayed reduction, every
stored entry is congruent to the true entry modulo $p$ and is at most
$(p - 1) + (s - 1) (p - 1)^2 \le s (p - 1)^2$.
:::

::: proof
Adding $f$ times the pivot row, with $f \equiv -\ell$, is congruent to
subtracting $\ell$ times it, so the congruence is preserved. An entry
$(i, j)$ is updated once for each pivot that is found while it is still
active, in a row above $i$ and a column to the left of $j$; there are at
most $\min(i, j) \le s - 1$ such pivots. Each update adds $f a$ with
$f \le p - 1$ and $a \le p - 1$, a reduced entry of the pivot row, to an
initial value at most $p - 1$.
:::

FLINT chooses the number of words from the bound $s (p - 1)^2$ of
[](algorithms/linear-algebra/dot-products#limbs)
(`_nmod_vec_dot_bound_limbs`), which bounds every entry by
[](#prop-delayed). The reduction of a word uses
the Lemire–Kaser–Kurz method when the bound is below $2^{32}$ and
Barrett's otherwise ([](algorithms/arithmetic/barrett)).

`nmod_mat_ple_classical_delayed` applies the same reductions to the
elimination of [](#thm-rotations), swapping the rows of the matrix and
of the buffer as `nmod_mat_lu_classical_delayed` does; the pivot of a
column is found while the column is reduced, as the row of smallest
index in $A$ with a nonzero entry (see [](#rank-profile)).

For moduli of several limbs, `mpn_mod_mat_lu_classical_delayed` applies
the same idea with unreduced entries of $2\ell + 1$ limbs
([delayed reduction for multiword moduli](algorithms/arithmetic/multiword-delayed-reduction#lu)).

## Crout elimination {#crout}

[](#alg-classical) updates the whole active submatrix after each pivot.
The Crout schedule performs the same operations in another order,
computing each entry of $L$ and $E$ only when it is needed [@DPS2017].
When column $j$ is processed, with $r$ pivots in the columns before it,
the entry $(i, j)$ of a row $i$ without pivot is
$a_{ij} - \sum_{k < r} \ell_{ik} e_{kj}$: one dot product of the
multipliers of row $i$ with column $j$ of $E$. When row $i$ becomes the
pivot row $r$, its entries to the right of column $j$ are
$a_{it} - \sum_{k < r} \ell_{ik} e_{kt}$: one vector-matrix product of
its multipliers with the rows of $E$. The other entries are not touched
before their column is processed.

`nmod_mat_ple_crout` uses this schedule with the pivoting of
[](#thm-rotations), swapping rows as
`nmod_mat_ple_classical_delayed` does. Its dot products and
vector-matrix products delay their reductions, with SIMD code
([](algorithms/linear-algebra/dot-products),
[](algorithms/linear-algebra/modular-matrix-kernels#thin)), so that no
buffer of unreduced entries is needed, whereas the delayed algorithm
needs two or three words per entry with scalar updates when
$\min(m, n) (p - 1)^2 \ge 2^{64}$. On the other hand, the products only
skip the zero multipliers before the first and after the last nonzero
one of each row, while the classical and delayed algorithms skip every
row whose multiplier is zero, which on some structured matrices saves
most of the work.

## Recursive decomposition {#recursive}

::: algorithm Recursive PLE decomposition {#alg-recursive}
**Input:** $A$, $m \times n$.

**Output:** a PLE decomposition of $A$ in place, and its rank.

1. If $m \le 1$ or $n \le 1$, use [](#alg-classical).
2. Let $n_1 = \lfloor n/2 \rfloor$ and decompose the left $m \times n_1$
   block recursively: $A_{\text{left}} = \Pi_1 L_1 E_1$ of rank $r_1$.
3. Apply $\Pi_1^{-1}$ to the rows of the right block. Write the first
   $r_1$ rows of the result as $A_{01}$ and the others as $A_{11}$, and
   the first $r_1$ rows of $L_1$ as $L_{00}$ and the others as $L_{10}$.
4. $A_{01} \gets L_{00}^{-1} A_{01}$ (triangular solve with unit
   diagonal), and $A_{11} \gets A_{11} - L_{10} A_{01}$ (matrix
   multiplication).
5. Decompose $A_{11}$ recursively, of rank $r_2$, apply its permutation to
   the rows of $L_{10}$, and move its $L$ factor to columns
   $r_1, \ldots, r_1 + r_2 - 1$.
6. Return $r_1 + r_2$.
:::

`nmod_mat_lu_recursive` calls `nmod_mat_lu` for the two halves, which
switches to the classical or delayed algorithm below a crossover depending
on the size of $p$ (see [](nmod_mat/lu#algorithms)).

The recursion performs the operations of [](#alg-classical), with the
same pivots, in another order: the pivots of the first $n_1$ columns
depend only on these columns, step 4 applies the eliminations of the
left half to the right half, and the row moves of step 5 act on the rows
that are not pivot rows of the left half. `nmod_mat_ple_recursive` uses
the upward rotations of [](#thm-rotations) in both halves, so that it
computes the decomposition of that theorem, the same as
`nmod_mat_ple_classical`.

::: theorem {#thm-recursive-cost}
If two matrices of dimensions $a \times b$ and $b \times c$ can be
multiplied with $O(a c \, b^{\omega - 2})$ operations for $b \le a, c$,
with $2 < \omega \le 3$, then [](#alg-recursive) costs
$O(m n r^{\omega - 2} + m n \log n)$ operations.
:::

::: proof
At depth $\ell$ of the recursion there are at most $2^\ell$ subproblems,
each with at most $m$ rows and $n / 2^\ell$ columns, with ranks
$r_i \le n/2^\ell$ summing to at most $r$. Step 4 is skipped when the
left half has rank 0; otherwise the triangular solve and the product cost
$O(m \, (n/2^\ell) \, r_i^{\omega - 2})$. Let $k \le \min(2^\ell, r)$ be
the number of nonzero $r_i$. Since $x \mapsto x^{\omega - 2}$ is concave,
$\sum_i r_i^{\omega - 2} \le k^{3 - \omega} r^{\omega - 2}$. For
$2^\ell \le r$ the level therefore costs
$O(m n \, 2^{-\ell(\omega - 2)} r^{\omega - 2})$, and for $2^\ell > r$,
using $r_i^{\omega - 2} \le r_i$, it costs $O(m n r \, 2^{-\ell})$. Both
are geometric series in $\ell$, with sums $O(m n r^{\omega - 2})$ and
$O(m n)$. The row permutations of step 3 and 5 cost $O(mn)$ per level,
and there are $O(\log n)$ levels.
:::

With classical multiplication ($\omega = 3$) this is $O(mnr)$, the cost
of [](#alg-classical); with Strassen–Winograd it is $O(mnr^{0.81})$. In
practice the gain comes as much from the speed of the multiplication
kernels as from the exponent.

## Triangular systems {#triangular}

`nmod_mat_solve_tril` and `nmod_mat_solve_triu` solve $TX = B$ for a
$d \times d$ triangular $T$ and a $d \times c$ matrix $B$. For
$\min(d, c) < 64$ they use substitution, one dot product per entry of
$X$, $O(d^2 c)$ operations. Otherwise they split $T$ into blocks of
$\lfloor d/2 \rfloor$ and $\lceil d/2 \rceil$ rows and use

$$
\begin{pmatrix} T_{00} & 0 \\ T_{10} & T_{11} \end{pmatrix}^{-1}
\begin{pmatrix} B_0 \\ B_1 \end{pmatrix} =
\begin{pmatrix} X_0 \\ T_{11}^{-1} (B_1 - T_{10} X_0) \end{pmatrix},
\qquad X_0 = T_{00}^{-1} B_0,
$$

and its upper triangular analogue, which costs
$O(c \, d^{\omega - 1})$ for $c \ge d$. `nmod_mat_solve_tril_right` and
`nmod_mat_solve_triu_right` solve $XT = B$ for $B$ with $d$ columns in the
same way: by substitution along the rows of $X$, with the triangle of $T$
transposed for the dot products, or by the transposed block formulas.

## Applications {#applications}

- *Solving* $AX = B$ for nonsingular $A$: a PLE decomposition with
  early exit on rank deficiency, a permutation of $B$, and two triangular
  solves. For $d \times d$ matrices $A$ and $B$, $O(d^\omega)$.
- *Inversion*: solving $AX = I$, $O(d^\omega)$.
- *Determinant*: $\det A = \pm \prod_i E_{ii}$, the sign being that of
  the permutation; $O(d^\omega)$.
- *Rank*: the number of pivots.
- *General systems* $AX = B$: with the pivot columns $c_i$, the solution
  with zero free variables is found by a triangular solve with $L$, a
  consistency check of the rows beyond the rank (one matrix product), and
  a triangular solve with the $r \times r$ matrix of pivot columns of $E$.
- *Reduced row echelon form*: with $E = (E_1 \mid E_2)$ up to a
  permutation of the columns, where $E_1$ is the upper triangular
  $r \times r$ matrix of pivot columns, the reduced form is
  $(I \mid E_1^{-1} E_2)$ with the columns put back in place; one
  triangular solve.
- *Nullspace*: from the reduced form $R$, one basis vector per free
  column $f$, with a 1 in position $f$ and $-R_{i f}$ in the position of
  the pivot of row $i$. These entries are those of $E_1^{-1} E_2$, so the
  basis is computed from a PLE decomposition with one triangular solve,
  without forming $R$.
- *Left nullspace*: from a PLE decomposition, the rows of
  $(-L_2 L_1^{-1} \mid I)$ with their columns put in the order of the
  rows of $A$, where $L_1$ and $L_2$ are the first $r$ and the other rows
  of $L$; one triangular solve with $L_1$ on the right. With the
  decomposition of [](#thm-rotations), the identity matrix is at the rows
  of $A$ outside its row rank profile, and the basis is the transpose of
  the basis of the nullspace of $A^T$ above.

## Division-free determinants {#det}

For an arbitrary modulus $N$, elimination can meet nonzero pivots that
are zero divisors. `nmod_mat_det` avoids this: it uses cofactor expansion
for dimensions up to 4, the division-free algorithm of Berkowitz
([](algorithms/linear-algebra/charpoly#berkowitz)) up to 8, and above it
PLE decomposition if $N$ is prime and the following triangularisation by
*Euclidean row operations* otherwise.

::: lemma {#lem-unit-xgcd}
For $0 \le a, b < N$, not both zero, let $g = \gcd(a, b)$. There are
integers $s, t$ with $s a - t b = g$ and $\gcd(s, N) = 1$.
:::

::: proof
Take any $s_0, t_0$ with $s_0 a - t_0 b = g$. Then $s = s_0 + j b/g$ and
$t = t_0 + j a/g$ satisfy $sa - tb = g$ for every integer $j$, and
$\gcd(s_0, b/g) = 1$ since $s_0 (a/g) - t_0 (b/g) = 1$. By the Chinese
remainder theorem there is a $j$ with $s_0 + j b/g \not\equiv 0$ modulo
every prime factor $q$ of $N$: for $q \nmid b/g$ one value of $j$ modulo
$q$ is excluded, and for $q \mid b/g$ none is, since then
$q \nmid s_0$.
:::

To clear an entry $b$ below a pivot $a$ in the same column,
`nmod_mat_det_howell` subtracts a multiple of the pivot row if $b$ is a
multiple of $a$ modulo $N$. Otherwise it replaces the pivot row by $s$
times itself minus $t$ times row $b$, with $s, t$ from
[](#lem-unit-xgcd) (found by increasing $j$ until $\gcd(s, N) = 1$). The
new pivot $g$ divides $b$, so the entry can then be cleared. The second
operation multiplies the determinant by the unit $s$, and the product of
these units is divided out at the end. The cost is $O(d^3)$ operations
and $O(d^2)$ extended gcd computations.

## Howell form {#howell}

Over $\Z/N\Z$ for composite $N$, the row echelon form of a matrix is not
unique, and its rows need not generate all elements of the row span with
leading zeros. The *Howell form* repairs both [@How1986; @StoMul1998].

::: definition {#def-howell}
A matrix $H$ over $\Z/N\Z$ is in *Howell form* if it is in row echelon
form with its zero rows last, every pivot is a divisor of $N$, every
entry above a pivot is reduced modulo that pivot, and, for every nonzero
row $i$, every element of the row span of $H$ whose entries vanish up to
and including the pivot column of row $i$ is a combination of the rows
below row $i$.
:::

Every matrix is row equivalent to a unique matrix in Howell form, so two
matrices span the same submodule of $(\Z/N\Z)^c$ if and only if their
Howell forms are equal [@How1986]. FLINT computes the *strong echelon
form* of Fieker and Hofmann [@FieHof2014], which is the Howell form with
zero rows inserted so that each pivot lies on the diagonal, and
`nmod_mat_howell_form` then moves the zero rows to the bottom. The
algorithm:

1. Triangularise $A$ with Euclidean row operations, as for the
   determinant but with a unimodular $2 \times 2$ transformation of the
   pivot row and the other row, so that the row span is preserved.
2. For each diagonal entry $a$ from top to bottom: multiply its row by a
   unit $u$ such that $ua$ divides $N$, and reduce the entries above it
   modulo it. If $a$ is a zero divisor, the row multiplied by
   $N / \gcd(a, N)$ is a new element of the span with a zero in this
   column; it is reduced into the rows below with the same Euclidean
   operations, which is what establishes the Howell property. If the
   diagonal entry is zero, the row itself is reduced in this way.

The matrix must have at least as many rows as columns, so that each
column has a row to hold its pivot. The cost is $O(c^2 m)$ operations
plus $O(c m)$ extended gcd computations for an $m \times c$ matrix. Fast
algorithms based on matrix multiplication exist [@StoMul1998; @Sto2000]
but are not implemented.
