---
title: Eigenvalues and eigenvectors
functions: acb_mat_approx_eig_qr acb_mat_eig_global_enclosure acb_mat_eig_enclosure_rump acb_mat_eig_simple acb_mat_eig_multiple
---

Computing certified eigenvalues of a complex matrix $A$ in FLINT takes
two steps. First, `acb_mat_approx_eig_qr` computes floating-point
approximations of the eigenvalues and eigenvectors by the QR algorithm,
without error bounds. Then one of the certification functions turns
these approximations into rigorous enclosures:

- `acb_mat_eig_simple` proves that $A$ has $n$ distinct eigenvalues and
  returns disjoint enclosures of them, together with enclosures of the
  left and right eigenvectors if wanted;
- `acb_mat_eig_multiple` allows multiple or clustered eigenvalues and
  returns enclosures with multiplicities;
- `acb_mat_eig_global_enclosure` returns a single radius $\varepsilon$
  such that every eigenvalue lies within $\varepsilon$ of one of the
  approximations, which is the cheapest way to certify statements about
  the whole spectrum;
- `acb_mat_eig_enclosure_rump` certifies one eigenvalue (or one cluster)
  and its eigenvector (or invariant subspace).

The certification functions accept approximations from any source, and
their output is valid whatever the quality of the approximations;
poor approximations only make them fail or give wide enclosures. For a
real matrix, convert it with `acb_mat_set_arb_mat` first.

::: summary algorithms/ball-arithmetic/eigenvalues :::

## Approximation by the QR algorithm {#qr}

`acb_mat_approx_eig_qr(E, L, R, A, tol, maxiter, prec)` works on the
midpoint matrix at `prec` bits: Householder reduction to Hessenberg form,
then the implicitly shifted QR iteration until the subdiagonal entries
are below `tol` (by default $2^{-\mathit{prec}}/(100n)$ relative to
their neighbours) or `maxiter` iterations (by default $14n +
\mathit{prec}/10$) have been done, and finally back-substitution in the
triangular Schur form for the eigenvectors. The return value only says
whether the iteration converged numerically. The cost is $O(n^3)$
operations at `prec` bits.

## Certification {#certification}

**All eigenvalues, simple.** `acb_mat_eig_simple` uses by default the
algorithm of van der Hoeven and Mourrain [@HM2017]: with $T$ the
midpoint of the approximate eigenvector matrix, it computes $D = T^{-1}
A T$ in ball arithmetic (by `acb_mat_solve`), splits it into its
diagonal $\Delta$ and off-diagonal part $H$, and checks with the
separation $\sigma = \min_{i \ne j} |\delta_i - \delta_j|$ and $\mu =
\|\Delta\|_\infty$ that $\|H\|_\infty$ is small enough compared to
$\sigma$ and $\mu$ for a perturbation argument to apply. The
eigenvalues are then the diagonal entries of $D$ inflated by
$\|H\|_\infty$, and enclosures of the eigenvectors follow from a
correction of $T$. The cost is $O(n^3)$. The alternative
`acb_mat_eig_simple_rump` certifies the eigenpairs one by one with
`acb_mat_eig_enclosure_rump` at a cost of $O(n^4)$. Both fail (return
0) if two enclosures overlap, which happens for multiple eigenvalues
and for eigenvalues closer than the attainable accuracy.

**Clusters.** `acb_mat_eig_multiple` first tries
`acb_mat_eig_simple_vdhoeven_mourrain`; if that fails, it groups
approximate eigenvalues that agree to about $3/4$ (and then $1/2$) of
the precision into clusters, certifies each cluster of size $k$ with
`acb_mat_eig_enclosure_rump` (a ball containing at least $k$
eigenvalues), and succeeds if the cluster enclosures are pairwise
disjoint, since they then contain exactly $k$ eigenvalues each. The
output repeats each cluster enclosure $k$ times.

**One eigenpair or cluster.** `acb_mat_eig_enclosure_rump` follows
Section 13.4 of [@Rum2010]: the eigenproblem $AX = X\lambda$ (or
$AX = XM$ for a cluster) is written as a fixed point problem for the
correction to the approximation, in which $k$ rows of the eigenvector
are frozen to normalise it; a preconditioned fixed point iteration with
$\varepsilon$-inflation looks for a set of corrections that is mapped
into itself, which proves the existence of an eigenvalue (or of $k$
eigenvalues) in the resulting enclosure. For a cluster, the output
matrix $J$ satisfies $AR - RJ \ni 0$, so that the Jordan structure of
the cluster is that of a matrix in $J$.

**The whole spectrum.** `acb_mat_eig_global_enclosure` implements
Theorem 1 of Miyajima [@Miy2010]: with $D$ the diagonal matrix of
approximate eigenvalues, $R$ the approximate eigenvectors and $Y
\approx R^{-1}$, every eigenvalue of $A$ lies within
$$\varepsilon = \frac{\|Y(AR - RD)\|_\infty}{1 - \|YR - I\|_\infty}$$
of a diagonal entry of $D$, provided that $\|YR - I\|_\infty < 1$. This
generalises the Bauer–Fike theorem to approximate eigenvectors. It
needs only a few matrix products, but it cannot tell which disk
contains how many eigenvalues.

The algorithms and proofs are in
[](algorithms/ball-arithmetic/eigenvalues).

## Practical advice {#advice}

- Compute the approximations at the same precision as the
  certification, or a little higher; the certification cannot produce
  enclosures much tighter than the approximations are accurate.
- If certification fails, increase the precision. If it keeps failing,
  the matrix may have multiple eigenvalues (use `acb_mat_eig_multiple`)
  or be defective (non-diagonalisable), in which case only cluster
  enclosures exist, and their width is about $\varepsilon^{1/m}$ for a
  Jordan block of size $m$ and an input accuracy $\varepsilon$.
- Badly scaled matrices may need balancing first, and a similarity
  transformation by a random matrix sometimes helps when an algorithm
  fails to converge on a well-scaled matrix.
- The left eigenvectors computed by `acb_mat_eig_simple` are normalised
  so that $L = R^{-1}$, which gives a certified diagonalisation $LAR =
  D$.

## Examples {#examples}

A $4 \times 4$ nonsymmetric matrix with four real eigenvalues. The
global enclosure gives a tiny radius for the spectrum as a whole;
`acb_mat_eig_simple` isolates the eigenvalues individually (their
imaginary parts are only known to be small: certification in complex
arithmetic cannot prove that an eigenvalue is real):

```c example name=acb_mat_eig_simple
#include <flint/acb.h>
#include <flint/acb_mat.h>

int main(void)
{
    acb_mat_t A, L, R;
    acb_ptr E, F;
    mag_t eps;
    slong i, j, n = 4, prec = 64;

    acb_mat_init(A, n, n);
    acb_mat_init(L, n, n);
    acb_mat_init(R, n, n);
    E = _acb_vec_init(n);
    F = _acb_vec_init(n);
    mag_init(eps);

    /* the Hilbert matrix with the top right entry replaced by -1 */
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            acb_set_ui(acb_mat_entry(A, i, j), i + j + 1);
            acb_inv(acb_mat_entry(A, i, j), acb_mat_entry(A, i, j), prec);
        }
    acb_set_si(acb_mat_entry(A, 0, n - 1), -1);

    acb_mat_approx_eig_qr(E, NULL, R, A, NULL, 0, prec);

    acb_mat_eig_global_enclosure(eps, A, E, R, prec);
    flint_printf("global enclosure radius: ");
    mag_printd(eps, 3);
    flint_printf("\n");

    if (acb_mat_eig_simple(F, L, R, A, E, R, prec))
    {
        _acb_vec_sort_pretty(F, n);
        for (i = 0; i < n; i++)
        {
            acb_printn(F + i, 12, 0);
            flint_printf("\n");
        }
    }

    acb_mat_clear(A);
    acb_mat_clear(L);
    acb_mat_clear(R);
    _acb_vec_clear(E, n);
    _acb_vec_clear(F, n);
    mag_clear(eps);
    return 0;
}
```
```output
global enclosure radius: 3.24e-18
[-0.0392285082075 +/- 2.44e-14] + [+/- 1.32e-16]*I
[-0.00114669579895 +/- 3.57e-15] + [+/- 1.33e-16]*I
[0.558734852147 +/- 1.64e-13] + [+/- 1.36e-16]*I
[1.15783082805 +/- 1.11e-13] + [+/- 1.39e-16]*I
```

The matrix $A = T \operatorname{diag}(1, 1, 2) T^{-1}$ below has the
double eigenvalue 1. `acb_mat_eig_simple` fails, and
`acb_mat_eig_multiple` returns the cluster enclosure twice:

```c example name=acb_mat_eig_multiple
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/acb.h>
#include <flint/acb_mat.h>

int main(void)
{
    fmpq_mat_t T, Tinv, D, Q;
    acb_mat_t A, R;
    acb_ptr E, F;
    slong i, n = 3, prec = 64;

    fmpq_mat_init(T, n, n);
    fmpq_mat_init(Tinv, n, n);
    fmpq_mat_init(D, n, n);
    fmpq_mat_init(Q, n, n);
    acb_mat_init(A, n, n);
    acb_mat_init(R, n, n);
    E = _acb_vec_init(n);
    F = _acb_vec_init(n);

    /* A = T diag(1, 1, 2) T^-1 has the double eigenvalue 1 */
    fmpq_set_si(fmpq_mat_entry(T, 0, 0), 1, 1);
    fmpq_set_si(fmpq_mat_entry(T, 0, 1), 1, 1);
    fmpq_set_si(fmpq_mat_entry(T, 1, 1), 1, 1);
    fmpq_set_si(fmpq_mat_entry(T, 1, 2), 1, 1);
    fmpq_set_si(fmpq_mat_entry(T, 2, 0), 1, 1);
    fmpq_set_si(fmpq_mat_entry(T, 2, 2), 1, 1);
    fmpq_set_si(fmpq_mat_entry(D, 0, 0), 1, 1);
    fmpq_set_si(fmpq_mat_entry(D, 1, 1), 1, 1);
    fmpq_set_si(fmpq_mat_entry(D, 2, 2), 2, 1);
    fmpq_mat_inv(Tinv, T);
    fmpq_mat_mul(Q, T, D);
    fmpq_mat_mul(Q, Q, Tinv);
    acb_mat_set_fmpq_mat(A, Q, prec);

    acb_mat_approx_eig_qr(E, NULL, R, A, NULL, 0, prec);

    flint_printf("simple: %d\n", acb_mat_eig_simple(F, NULL, NULL, A, E, R, prec));
    flint_printf("multiple: %d\n", acb_mat_eig_multiple(F, A, E, R, prec));
    for (i = 0; i < n; i++)
    {
        acb_printn(F + i, 10, 0);
        flint_printf("\n");
    }

    fmpq_mat_clear(T);
    fmpq_mat_clear(Tinv);
    fmpq_mat_clear(D);
    fmpq_mat_clear(Q);
    acb_mat_clear(A);
    acb_mat_clear(R);
    _acb_vec_clear(E, n);
    _acb_vec_clear(F, n);
    return 0;
}
```
```output
simple: 0
multiple: 1
[1.000000000 +/- 4.41e-19] + [+/- 2.17e-19]*I
[1.000000000 +/- 4.41e-19] + [+/- 2.17e-19]*I
[2.000000000 +/- 3.26e-18] + [+/- 1.63e-18]*I
```

Since this matrix is diagonalisable, the cluster enclosure is about as
tight as the enclosure of a simple eigenvalue; for a Jordan block, it
would be much wider.
