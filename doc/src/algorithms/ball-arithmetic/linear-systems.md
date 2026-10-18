---
title: Verified linear algebra
summary: How arb_mat and acb_mat solve linear systems, invert matrices and compute determinants with rigorous enclosures, by interval Gaussian elimination, by preconditioning with an approximate inverse, and by certifying an approximate solution with a residual bound.
complexity: $O(n^3)$ ball operations for interval LU decomposition; preconditioned solving costs an approximate LU decomposition in floating-point arithmetic, two or three ball matrix products and an interval LU decomposition or a norm bound, which is $O(n^\omega)$ with fast ball matrix multiplication.
order: 90
---

Let $A$ be an $n \times n$ matrix and $B$ an $n \times m$ matrix with
ball entries. Solving $AX = B$ rigorously means computing a ball matrix
$X$ that contains $\tilde A^{-1} \tilde B$ for every pair of matrices
$\tilde A$, $\tilde B$ with entries in the balls of $A$ and $B$, or
reporting failure if this cannot be done (in particular if the balls of
$A$ contain a singular matrix). There are two basic approaches: carry out
an exact algorithm (Gaussian elimination) in ball arithmetic, or compute
an approximate solution in floating-point arithmetic and then prove an
error bound for it [@Rum2010]. FLINT uses both, depending on the
dimension and the precision.

## Interval Gaussian elimination {#lu}

Gaussian elimination with partial pivoting in ball arithmetic
(`arb_mat_lu`) chooses as pivot in each column an entry whose ball does
not contain zero, among those the one with the largest midpoint, and
fails if all candidates contain zero. The resulting balls contain the
$LU$ factors of every matrix $P\tilde A$ in the balls of $PA$ (with the
same permutation), so triangular solving in ball arithmetic
(`arb_mat_solve_lu_precomp`) then gives a correct enclosure of the
solution. Success proves that every matrix in the balls is nonsingular.

The weakness is radius growth. Each elimination step computes $a_{jk}
- (a_{ji}/a_{ii}) a_{ik}$ and the radius of the result includes the
radii of all three inputs multiplied by the magnitudes of the others; the
errors introduced in one step are treated as independent in the next,
although they are strongly correlated. As a result, the radii can grow
exponentially with $n$ even for well-conditioned matrices, so that
$O(n)$ bits or more are lost [@HS1967]. If the precision $p$ is large
compared to $n$, a loss of $O(n)$ bits is harmless, and interval elimination is the fastest
method; this is why `arb_mat_solve` uses it when $n \le 4$ or $p > 10n$.
The recursive LU decomposition (`arb_mat_lu_recursive`, for $n \ge 8$)
performs most of the work in ball matrix products and has the same
behaviour.

## Preconditioning {#precond}

Hansen and Smith [@HS1967] observed that the radius growth disappears if
the system is first multiplied by an approximate inverse $R \approx
A^{-1}$ (a matrix with exact floating-point entries): $RA$ is close to
the identity matrix, all pivots are close to 1 and the off-diagonal
entries are small, so the elimination of $RAX = RB$ in ball arithmetic
loses only a few bits. The solution set is unchanged, since $R$ is
nonsingular whenever $RA$ is.

::: algorithm Preconditioned solving {#alg-precond}
Input: $A$, $B$ and a precision $p$.

1. Compute an approximate LU decomposition of the midpoint matrix of $A$
   in floating-point arithmetic (`arb_mat_approx_lu`) and from it an
   approximate inverse $R$ (`arb_mat_approx_solve` with $B = I$). Fail
   if a pivot is zero.
2. Compute the ball matrices $RA$ and $RB$.
3. If $m < 0.1 n + 1$, solve $(RA) X = RB$ by interval Gaussian
   elimination.
4. Otherwise, compute an upper bound $\delta$ for $\|I - RA\|_\infty$.
   If $\delta < 1$, return $X = RB$ with each entry of column $j$ enlarged
   by $\frac{\delta}{1 - \delta} \max_i |(RB)_{ij}|$ (see
   [](#thm-neumann)); otherwise solve as in step 3.
:::

::: theorem {#thm-neumann}
Let $M$ be an $n \times n$ matrix with $\|I - M\|_\infty \le \delta < 1$.
Then $M$ is nonsingular, and for every vector $y$, the solution $x$ of $M
x = y$ satisfies $\|x - y\|_\infty \le \frac{\delta}{1 - \delta}
\|y\|_\infty$.
:::

::: proof
With $E = I - M$, the Neumann series $M^{-1} = \sum_{k \ge 0} E^k$
converges since $\|E\|_\infty < 1$, and $x - y = (M^{-1} - I) y =
\sum_{k \ge 1} E^k y$ has norm at most $\sum_{k \ge 1} \delta^k
\|y\|_\infty$.
:::

Applied to each matrix $\tilde M = R \tilde A$ in the balls of $RA$ and
each column $\tilde y$ of $R \tilde B$ in the balls of $RB$, with $\delta$
bounding the norm of all $I - \tilde M$ simultaneously, the theorem shows
that the solution lies in the ball around $\tilde y$ of radius
$\frac{\delta}{1-\delta} \|\tilde y\|_\infty$, which is contained in the
returned ball. This variant costs two matrix products (forming $RA$ and
$RB$) instead of an $O(n^2 m)$ triangular solve with ball entries and is
used when there are many right-hand sides, for example by `arb_mat_inv`,
which solves $AX = I$. The cost of preconditioning is dominated by the
three matrix multiplications ($R$ from the LU decomposition, $RA$, $RB$),
which use the fast ball matrix multiplication of
[](algorithms/ball-arithmetic/matrix-multiplication).

The radius of the result is about $\kappa(A) \cdot 2^{-p}$ relative to
the solution, where $\kappa(A) = \|A\| \|A^{-1}\|$ is the condition
number, and the method succeeds as soon as the precision exceeds $\log_2
\kappa(A)$ by a few bits, independently of $n$. The interval elimination
of step 3 is rarely the bottleneck because $RA$ is almost diagonal.

## Certifying an approximate solution {#preapprox}

If an approximate solution $T$ and an approximate inverse $R$ are already
known, for example from a previous call at lower precision or from a
floating-point solver, they can be certified directly.

::: theorem {#thm-rump}
Let $A$, $R$ be $n \times n$ matrices and $T$, $B$ be $n \times m$
matrices with $\|I - RA\|_\infty \le \delta < 1$. Then $A$ is
nonsingular, and each column $x$ of $A^{-1} B$ and the corresponding
columns $t$ of $T$ and $c$ of $C = R(AT - B)$ satisfy
$$\|x - t\|_\infty \le \frac{\|c\|_\infty}{1 - \delta}.$$
:::

::: proof
$RA$ is nonsingular by [](#thm-neumann), so $A$ is. Since $RA(t - x) =
R(At - b) = c$ and $\|(RA)^{-1}\|_\infty \le \sum_k \delta^k = 1/(1 -
\delta)$, the bound follows [@Rum2010, Theorem 10.2].
:::

`arb_mat_solve_preapprox` evaluates $\delta$ and $C$ in ball arithmetic
and returns $T$ with these error bounds. The residual $AT - B$ is small
when $T$ is accurate, so the radius is of the order of the error of $T$
plus the rounding errors in the residual, and the result can be much
tighter than the one obtained by solving from scratch at the same
precision. The cost is three matrix products. This is the approach to
use for iterative refinement: compute $T$ by a fast floating-point
method, improve it if needed, and certify it at the end.

## Determinants {#det}

`arb_mat_det` uses explicit formulas for $n \le 3$, the product of the
diagonal for triangular matrices, and otherwise either interval
elimination (for $n \le 10$ or $p > 10n$; $n \le 14$ for `acb_mat_det`)
or preconditioning. The interval elimination version `arb_mat_det_lu`
multiplies the pivots; if the elimination stops at step $r$ because all
candidate pivots contain zero, the determinant of the remaining $(n - r)
\times (n - r)$ block $S$ is bounded by Hadamard's inequality,
$$|\det S| \le \prod_i \|s_i\|_2,$$
where $s_i$ are the rows of $S$, so that the result is still a valid
(if wide) enclosure rather than a failure.

::: algorithm Preconditioned determinant {#alg-det}
1. Compute an approximate decomposition $PA \approx LU$ in floating-point
   arithmetic and approximate inverses $L' \approx L^{-1}$ (exactly unit
   lower triangular) and $U' \approx U^{-1}$ (exactly upper triangular).
2. Compute the ball matrix $M = L' P A U'$, which is close to $I$, and
   $\varepsilon \ge \max_i \sum_j |M_{ij} - \delta_{ij}|$.
3. Return $\pm \det(M) / \prod_i U'_{ii}$, where $\det(M)$ is enclosed
   by the ball $1 \pm (e^{n \varepsilon} - 1)$ and the sign is that of
   the permutation $P$.
4. If the result contains zero, also compute `arb_mat_det_lu` and
   return the narrower of the two.
:::

The enclosure of $\det M$ follows from Gershgorin's theorem: every
eigenvalue $\lambda$ of $M$ satisfies $|\lambda - 1| \le
\varepsilon$, so $|\det M - 1| = |\prod_i \lambda_i - 1| \le (1 +
\varepsilon)^n - 1 \le e^{n\varepsilon} - 1$. Since $\det L' = 1$ and
$\det U' = \prod_i U'_{ii}$ exactly, $\det A = \pm \det M / \det U'$.

## Symmetric positive definite matrices {#spd}

For a symmetric positive definite matrix, `arb_mat_cho` computes the
Cholesky factor $L$ with $A = LL^T$ in ball arithmetic and fails if a
diagonal entry $a_{ii} - \sum_k l_{ik}^2$ is not provably positive;
success proves that every symmetric matrix in the balls is positive
definite. `arb_mat_ldl` computes $A = LDL^T$ without square roots.
These decompositions need no pivoting, and the radius growth is moderate
since the entries of $L$ are bounded by $\sqrt{a_{ii}}$. They are used
by `arb_mat_spd_solve` and `arb_mat_spd_inv` (and the `_precomp`
variants), which do not precondition.

## Use in FLINT {#use}

- `arb_mat_solve` and `acb_mat_solve` choose between `*_solve_lu`
  ([](#lu)) and `*_solve_precond` ([](#alg-precond)) as described above;
  `arb_mat_inv` and `acb_mat_inv` call them with $B = I$;
  `*_solve_preapprox` implements [](#thm-rump). See [](arb_mat/solve).
- `arb_mat_det`, `arb_mat_det_lu`, `arb_mat_det_precond` and the complex
  versions ([](#det)).
- The approximate functions `arb_mat_approx_lu`, `arb_mat_approx_solve`,
  `arb_mat_approx_inv` and `arb_mat_approx_mul` work with midpoints only
  and give no error bounds.
- Generic code: `gr_mat_det` over the ball contexts calls `arb_mat_det`
  and `acb_mat_det`, and `gr_mat_mul` calls the ball matrix
  multiplication, but `gr_mat_nonsingular_solve` and `gr_mat_inv` use
  generic elimination without preconditioning ([](gr_mat)).
- Eigenvalue enclosures ([](algorithms/ball-arithmetic/eigenvalues))
  use the same idea of certifying approximate results with norm bounds.
