---
title: Eigenvalue enclosures
summary: How acb_mat computes rigorous enclosures of the eigenvalues and eigenvectors of complex matrices by certifying floating-point approximations from the QR algorithm, with a global perturbation bound, a Gershgorin test after an approximate diagonalisation, and Rump's fixed-point method for eigenvalue clusters.
complexity: $O(n^3)$ floating-point operations for the approximation by Hessenberg reduction and the shifted QR algorithm, plus a few ball matrix products and one ball linear solve for certification ($O(n^\omega)$ with fast ball matrix multiplication); Rump's method costs $O(n^3)$ per eigenvalue cluster.
order: 100
---

Computing eigenvalues is a nonlinear problem, and ball arithmetic
applied directly to an eigenvalue algorithm such as the QR iteration
gives useless enclosures because the radii grow with every iteration.
The functions of [](acb_mat) therefore follow the approximate-and-certify
approach of verification methods [@Rum2010]: compute approximate
eigenvalues $\tilde\lambda_i$ and an approximate eigenvector matrix
$\tilde R$ by a floating-point algorithm without error control, and then
prove, with a short computation in ball arithmetic, that the true
eigenvalues lie within computed distances of the approximations. All
statements below hold for every matrix in the balls of $A$.

## Approximation {#approximation}

`acb_mat_approx_eig_qr` reduces the midpoint matrix to upper Hessenberg
form by Householder reflections and runs the implicitly shifted QR
algorithm with deflation to obtain a Schur form $A \approx Q T Q^*$
[@GVL1996]; the diagonal of $T$ gives the eigenvalues, and the
eigenvectors are obtained by back substitution in $T$ and multiplication
by $Q$. The implementation is adapted from mpmath's `eigen.py`. An
off-diagonal entry is treated as zero when it is smaller than the
tolerance, by default $2^{-p}/(100n)$ relative to the neighbouring
diagonal entries, and the iteration gives up if no deflation occurs
within $14n + p/10$ QR steps (by default); the function then returns 0,
but still writes its current approximations, which may be usable. All arithmetic is done on midpoints at
$p$ bits.

## A global enclosure {#global}

::: theorem {#thm-global}
Let $\tilde D = \operatorname{diag}(\tilde\lambda_1, \ldots,
\tilde\lambda_n)$, let $\tilde R$ and $Y$ be $n \times n$ matrices, and
suppose that $\|I - Y\tilde R\|_\infty \le \delta < 1$. Let
$$\varepsilon = \frac{\|Y(A\tilde R - \tilde R \tilde D)\|_\infty}{1 - \delta}.$$
Then every eigenvalue of $A$ lies in one of the disks $|z -
\tilde\lambda_i| \le \varepsilon$.
:::

::: proof
$Y\tilde R$ is nonsingular by the Neumann series, so $\tilde R$ is
nonsingular, and $\tilde R^{-1} = (Y \tilde R)^{-1} Y$ with
$\|(Y\tilde R)^{-1}\|_\infty \le 1/(1-\delta)$. The matrix $A$ is
similar to
$$\tilde R^{-1} A \tilde R = \tilde D + F, \qquad F = (Y\tilde R)^{-1} Y (A\tilde R - \tilde R\tilde D),$$
with $\|F\|_\infty \le \varepsilon$. By Gershgorin's theorem, every
eigenvalue of $\tilde D + F$ lies in a disk around $\tilde\lambda_i +
F_{ii}$ of radius $\sum_{j \ne i} |F_{ij}|$, which is contained in the
disk around $\tilde\lambda_i$ of radius $\sum_j |F_{ij}| \le
\varepsilon$.
:::

This is the bound of Miyajima [@Miy2010] in its simplest form.
`acb_mat_eig_global_enclosure` computes $Y$ as an approximate inverse of
$\tilde R$ and evaluates $\delta$ and $\varepsilon$ in ball arithmetic.
The bound holds even if $A$ has multiple eigenvalues or is not
diagonalisable (if $\tilde R$ is far from singular), but it does not say
which disk contains which eigenvalue, or how many; it is used to obtain a
quick inclusion region or as a starting point.

## Simple eigenvalues {#simple}

If the eigenvalues are simple and well separated, the approximate
diagonalisation can be certified for each eigenvalue individually.
`acb_mat_eig_simple_vdhoeven_mourrain` (called by `acb_mat_eig_simple`)
follows van der Hoeven and Mourrain [@HM2017]:

::: algorithm Certifying a diagonalisation {#alg-simple}
1. Let $T$ be the midpoint matrix of $\tilde R$. Compute the ball matrix
   $D = T^{-1} A T$ by solving $T D = AT$ with `acb_mat_solve`; fail if
   the solving fails.
2. Write $D = \Delta + H$ with $\Delta$ diagonal and $H$ zero on the
   diagonal, and compute $\eta = \|H\|_\infty$, $\mu = \|\Delta\|_\infty$,
   the separation $\sigma = \min_{i \ne j} |\Delta_{ii} - \Delta_{jj}|$
   and $\alpha = \min(\sigma/(6\mu), 1/4)$.
3. If $\eta \le \alpha\mu/8$ and $\eta \le \alpha\sigma/8$, return the
   eigenvalue enclosures $\Delta_{ii} \pm \eta$ (as boxes with real and
   imaginary radius $\eta$), provided that they are pairwise disjoint,
   and, if requested, the right eigenvectors $T(I + E)$ with entries of
   $E$ bounded by $3\eta/\sigma$ and the left eigenvectors as the inverse
   of this matrix.
:::

::: theorem {#thm-simple}
If [](#alg-simple) succeeds, each returned box contains exactly one
eigenvalue of $A$, and these are all the eigenvalues; in particular, all
eigenvalues are simple.
:::

::: proof
$D$ contains $T^{-1} \tilde A T$ for every $\tilde A$ in the balls of
$A$, which has the same eigenvalues as $\tilde A$. By Gershgorin's
theorem, its eigenvalues lie in the union of the disks around the
diagonal entries $D_{ii}$ with radii $\sum_{j \ne i} |D_{ij}| \le \eta$,
and each connected component of this union formed by $k$ disks contains
exactly $k$ eigenvalues. Each disk is contained in the corresponding
box, and the boxes are disjoint, so each disk contains exactly one
eigenvalue.
:::

The conditions on $\eta$ in step 3 are those of the theorem of van der
Hoeven and Mourrain [@HM2017], which shows that there is an exact
diagonalisation $A = (T(I+E)) \Lambda (T(I+E))^{-1}$ with $\|E\| \le
3\eta/\sigma$; they are only needed for the eigenvectors. The cost is
dominated by the ball matrix product $AT$ and the preconditioned solve
([](algorithms/ball-arithmetic/linear-systems)). The enclosures have
radius about $\|H\|$, which is of the order of the error of the
approximation times the condition number of $T$; a slightly better
enclosure could be obtained with the centre of the diagonal entries of
$D$ instead of the ball, as in [@HM2017].

## Clusters and multiple eigenvalues {#clusters}

A multiple eigenvalue cannot be certified as simple, and an eigenvalue
of a nontrivial Jordan block is perturbed by $O(\epsilon^{1/k})$ when
$A$ is perturbed by $\epsilon$, so its approximations are inaccurate. For
such clusters, Rump's method [@Rum2010, Section 13.4], also used in
INTLAB's `verifyeig`, encloses a group of $k$ eigenvalues together.

Let $\tilde\lambda$ be an approximate eigenvalue of the cluster and
$\tilde X$ an $n \times k$ matrix of approximations of a basis of the
invariant subspace. Choose $k$ row indices $v$ (the $k$ rows of $\tilde
X$ with the largest absolute value when $k = 1$, the last $k$ rows
otherwise); let $U$ select the other rows. The unknowns are a
correction $Y$ ($n \times k$) whose rows in $U$ correct $\tilde X$ and
whose rows in $v$ form a $k \times k$ matrix $M$ with $A X = X
(\tilde\lambda I + M)$ for $X = \tilde X + U U^T Y$. With an approximate
inverse $R$ of $(A - \tilde\lambda I) U U^T - \tilde X V^T$, this
equation becomes the fixed-point problem
$$Y = g(Y) = -R(A - \tilde\lambda I)\tilde X + \bigl(I - R((A - \tilde\lambda I)UU^T - \tilde X V^T)\bigr) Y + R (UU^T Y)(V^T Y).$$

::: theorem Rump {#thm-rump}
If $\mathcal{Y}$ is a ball matrix such that $g(\mathcal{Y})$ is
contained in the interior of $\mathcal{Y}$, then there is a fixed point
$Y \in \mathcal{Y}$, the matrix $\tilde\lambda I_k + V^T Y$ has $k$
eigenvalues that are eigenvalues of $A$ (counted with multiplicity), and
$X = \tilde X + UU^TY$ spans the corresponding invariant subspace.
Every eigenvalue of $\tilde\lambda I_k + V^TY$ lies in the disk around
$\tilde\lambda$ of radius $\rho(|V^T \mathcal{Y}|)$, the spectral radius
of the matrix of absolute values [@Rum2010, Section 13.4].
:::

`acb_mat_eig_enclosure_rump` implements this with
$\varepsilon$-inflation: starting from $\mathcal{Y}_0 = g(0)$, it
enlarges $\mathcal{Y}$ by $2^{-3+2j} \max|\mathcal{Y}|$ in iteration $j$
and evaluates $g(\mathcal{Y})$ in ball arithmetic, for at most $5 + b$
iterations, where $b$ is the bit length of $p$. The spectral radius in the eigenvalue bound is
bounded by the Frobenius norm and by the Collatz–Wielandt bound $\rho(N)
\le \max_i (Nw)_i / w_i$ for positive vectors $w$, using the absolute
values of approximate eigenvectors of $N$ as $w$.

`acb_mat_eig_multiple_rump` groups the approximate eigenvalues into
clusters of points within distance $2^{-3p/4} \max_i |\tilde\lambda_i|$
of the first point of the cluster (and on failure $2^{-p/2} \max_i
|\tilde\lambda_i|$), encloses each cluster of size $k$ with
`acb_mat_eig_enclosure_rump`, and succeeds if the enclosures are
pairwise disjoint. Since each enclosure contains as many eigenvalues as
its cluster has members and the cluster sizes add up to $n$, every
eigenvalue of $A$ is accounted for, and each enclosure is returned once
per member. `acb_mat_eig_multiple` first tries the simple-eigenvalue
algorithm, which is cheaper, and uses Rump's method only on failure.

As the theorem requires, `acb_mat_eig_enclosure_rump` accepts
$\mathcal{Y}$ only when every entry of $g(\mathcal{Y})$ lies in the
interior of the corresponding entry of $\mathcal{Y}$ (tested with
`acb_contains_interior`); balls that touch the boundary are rejected.

## Use in FLINT {#use}

- `acb_mat_approx_eig_qr` ([](#approximation)),
  `acb_mat_eig_global_enclosure` ([](#thm-global)),
  `acb_mat_eig_simple`, `acb_mat_eig_simple_vdhoeven_mourrain`,
  `acb_mat_eig_simple_rump` (Rump's method for each eigenvalue
  separately), `acb_mat_eig_enclosure_rump`, `acb_mat_eig_multiple`,
  `acb_mat_eig_multiple_rump`; see [](acb_mat/eigenvalues).
- There are no separate functions for real matrices; they are converted
  with `acb_mat_set_arb_mat`.
- `gr_mat_diagonalization` over `gr_ctx_init_complex_acb` calls
  `acb_mat_approx_eig_qr` and `acb_mat_eig_simple` ([](gr_mat)).
