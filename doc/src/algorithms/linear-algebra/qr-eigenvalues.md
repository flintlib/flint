---
title: The QR algorithm for eigenvalues
summary: Approximate eigenvalues and eigenvectors of a dense matrix over a floating-point ring by Householder reduction to Hessenberg form, the implicitly shifted QR iteration with deflation, back substitution in the Schur form, and Schur reordering by plane rotations.
complexity: $O(n^3)$ operations: about $\frac{10}{3}n^3$ for the Hessenberg reduction and $O(n^2)$ per QR step, with $O(n)$ steps in practice; about twice as much with eigenvectors.
order: 55
---

The eigenvalues of a matrix are the roots of its characteristic
polynomial, but computing them as such is numerically hopeless: the
roots of a polynomial are far more sensitive to its coefficients than
the eigenvalues are to the entries of the matrix. Floating-point
eigensolvers instead transform the matrix by unitary similarities,
which do not change the eigenvalues and do not amplify rounding errors,
until it is triangular. The QR algorithm of Francis and Kublanovskaya,
with the refinements of the following sections, does this in $O(n^3)$
operations and is backward stable: the computed Schur form is the exact
Schur form of a matrix $A + E$ with $\|E\| \le c\, n\, u \|A\|$ for the
unit roundoff $u$ and a modest constant $c$ [@GVL1996, §7.5]. It is the
algorithm of LAPACK's `zgeev`, and of `gr_mat_eig_qr` and
`acb_mat_approx_eig_qr` in FLINT; the implementation follows mpmath's
`eigen.py` (by Timo Hartmann), with the plane rotations of LAPACK.

Throughout, $A$ is an $n \times n$ complex matrix, $A^*$ its conjugate
transpose, and $p$ the working precision, $u = 2^{-p}$.

## Reduction to Hessenberg form {#hessenberg}

A matrix is upper Hessenberg if $a_{ij} = 0$ for $i > j + 1$. Every
matrix is unitarily similar to a Hessenberg matrix $H = Q^* A Q$, which
$n - 2$ Householder reflections $I - \tau v v^*$ compute in
$\frac{10}{3} n^3$ operations, each reflection zeroing one column below
its subdiagonal and being applied from both sides. The Hessenberg form
is preserved by the QR iteration and makes each of its steps cost
$O(n^2)$ instead of $O(n^3)$.

`gr_mat_hessenberg_reduce` works from the last row upwards and stores
each reflector in the row it zeroes (rather than in the column, as
LAPACK's `zgehd2`), so that both applications of a reflector are
unit-stride vector operations, dot products and `_gr_vec_addmul_scalar`
calls on rows, which floating-point rings implement with SIMD code. The
reflectors stay below the subdiagonal; `gr_mat_hessenberg_expand_q`
accumulates them into $Q$ when eigenvectors are wanted.

## The shifted QR iteration {#qr}

A QR step with shift $\sigma$ factors $H - \sigma I = QR$ and replaces
$H$ by $RQ + \sigma I = Q^* H Q$. If $\sigma$ is close to an eigenvalue
$\lambda$, the last subdiagonal entry $h_{n,n-1}$ is reduced by a
factor of about $|\lambda - \sigma| / |\lambda' - \sigma|$, where
$\lambda'$ is the next closest eigenvalue; with the shifts below the
convergence is quadratic (cubic for normal matrices), so that a few
steps per eigenvalue suffice.

The step is computed implicitly: the plane rotation $G_1$ that maps the
first column of $H - \sigma I$ to a multiple of $e_1$ is applied to $H$
from both sides, which creates a nonzero entry (the bulge) below the
subdiagonal, and rotations $G_2, \ldots, G_{n-1}$ chase the bulge down
and out of the matrix. By the implicit Q theorem the result equals the
explicit step. Each rotation $\begin{pmatrix} c & s \\ -\bar s & c
\end{pmatrix}$ has a real cosine $c$ and $|c|^2 + |s|^2 = 1$, as in
LAPACK's `zlartg`: the diagonal of each $2 \times 2$ transformation is
then exactly real, which a ring can exploit in its `vec_apply_2x2`
method (a quarter of the real multiplications of a general complex
rotation are saved).

The shift is the Wilkinson shift, the eigenvalue of the trailing
$2 \times 2$ block of the active part closest to its last diagonal
entry. Since it can stall on special matrices (permutation-like
structures, for which it is exactly symmetric between two eigenvalues),
exceptional shifts are used after 10, 20 and 29 steps without
deflation: the last subdiagonal entry, its absolute value and the norm
of the matrix, cycling every 30 steps.

## Deflation {#deflation}

When a subdiagonal entry $h_{k+1,k}$ becomes negligible, it is set to
zero, which splits $H$ into two Hessenberg blocks whose eigenvalues can
be computed separately. The test is
$$|h_{k+1,k}| < \varepsilon \max(|h_{kk}|, |h_{k+1,k+1}|),$$
with $\varepsilon = u/(100n)$ by default and the right-hand side
replaced by $\varepsilon \|H\|_F / n$ when the diagonal entries are
smaller than that, so that small eigenvalues do not prevent deflation.
The iteration works on the trailing block until a block of size one or
two has converged, then on the remaining leading part. Since the
ring's arithmetic only approximates the precision $p$ it reports, a
tolerance given by the caller should not be much smaller than $u$: the
iteration may otherwise fail to converge.

The eigenvalues are the diagonal entries of the resulting upper
triangular matrix $T = V^* A V$ (a complex Schur form). In a real ring
all shifts are real, so the iteration converges only if all eigenvalues
are real; the double-shift variant with $2 \times 2$ diagonal blocks for
complex conjugate pairs (LAPACK's `dhseqr`) is not implemented.

## Eigenvectors {#eigenvectors}

The right eigenvector of $T$ for the eigenvalue $t_{kk}$ has $x_k = 1$,
$x_i = 0$ for $i > k$ and
$$x_i = -\frac{1}{t_{ii} - t_{kk}} \sum_{j=i+1}^{k} t_{ij} x_j, \quad i = k-1, \ldots, 0,$$
and the eigenvector of $A$ is $V x$; left eigenvectors are computed
likewise from the rows. As in LAPACK's `ztrevc`, a difference
$t_{ii} - t_{kk}$ smaller than a threshold $s$ is replaced by $s$, so
that repeated eigenvalues give finite (though not independent)
eigenvectors, and the partial vector is scaled down whenever its
largest entry exceeds $2^{p/2}$, against overflow. The threshold is
$s = n\, 2^{-29p}$ as in `acb_mat_approx_eig_qr`, or, in a ring with a
bounded exponent range where this underflows (the formats of `dfloat`),
$n\, 2^{-p} \max |t_{ij}|$: a perturbation at the level of the backward
error.

The eigenvectors are accurate when the eigenvalues are well separated
relative to the norm; for nearly multiple eigenvalues the individual
eigenvectors are ill-conditioned, although the invariant subspace they
span may be well-conditioned.

## Schur reordering {#reordering}

Adjacent diagonal entries of a complex Schur form can be exchanged by
one plane rotation: the eigenvector of
$\begin{pmatrix} t_{11} & t_{12} \\ 0 & t_{22} \end{pmatrix}$ for
$t_{22}$ is $(t_{12}, t_{22} - t_{11})$, and the rotation whose first
column is this vector, normalised, moves $t_{22}$ to the top. The
rotation is applied to the rows and columns of $T$ and to $V$, the
exchanged diagonal entries are set exactly, and the new subdiagonal
entry, which is zero up to rounding, is set to zero. A sequence of
exchanges moves any eigenvalue to any position (LAPACK's `ztrexc`),
which serves to compute invariant subspaces for selected eigenvalues
and is the building block of aggressive early deflation in multishift
QR algorithms.

## Accuracy and precision {#accuracy}

By backward stability, the computed eigenvalues are the exact
eigenvalues of $A + E$ with $\|E\| = O(n u \|A\|)$, so a simple
eigenvalue $\lambda$ with left and right eigenvectors $y, x$ is computed
with an error of about $\|E\| \|x\| \|y\| / |y^* x|$, the
condition number times the backward error (for a normal matrix,
$|y^* x| = \|x\| \|y\|$ and the error is $O(n u \|A\|)$). An eigenvalue
of multiplicity $m$ in a Jordan block is perturbed by about
$\|E\|^{1/m}$. In arbitrary precision the remedy for ill-conditioning is
simply to raise $p$; the certification methods
([](algorithms/ball-arithmetic/eigenvalues)) tell whether it suffices.

In the multi-double formats of `dfloat` the operations are accurate to
a few units in the last place of their $53N$ bits rather than correctly
rounded, and the computed eigenvalues of random matrices are typically
one to three bits less accurate than with `acf` at the same precision;
`acb_mat_approx_eig_qr` therefore asks for a few bits of margin before
using them ([](algorithms/ball-arithmetic/linear-systems#backends)).

## Use in FLINT {#use}

- `gr_mat_eig_qr` and the functions it is built from
  (`gr_mat_hessenberg_reduce`, `gr_mat_hessenberg_expand_q`,
  `gr_mat_hessenberg_qr_step`, `gr_mat_hessenberg_qr`,
  `gr_mat_eig_triu_r`, `gr_mat_eig_triu_l`, `gr_mat_schur_reorder`),
  over any ring of approximate complex numbers; see [](gr_mat/eig).
- `acb_mat_approx_eig_qr` is `gr_mat_eig_qr` over `dfloat` or
  `nfloat_complex` when the precision and the entries allow it, and over
  `acf_t` otherwise; its results feed the certified eigenvalue functions of
  `acb_mat` ([](algorithms/ball-arithmetic/eigenvalues)).
- The plane rotations are the ring method `_gr_vec_apply_2x2`; the
  complex `dfloat` rings implement it with SIMD code, where it takes
  about three quarters of the time of the eigensolver.
