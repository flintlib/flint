---
title: Numerical eigenvalues
functions: gr_mat_eig_qr gr_mat_hessenberg_reduce gr_mat_hessenberg_qr gr_mat_eig_triu_r gr_mat_schur_reorder
---

The functions on this page compute eigenvalues and eigenvectors of a
matrix of floating-point numbers, in any ring of approximate real or
complex numbers: `acf` and `nfloat_complex` at any precision, the
complex `dfloat` rings (`d1c` to `d4c`, double to quadruple-double
precision), and their real counterparts for matrices with real
eigenvalues. `acb_mat_approx_eig_qr` is a wrapper of `gr_mat_eig_qr`
that chooses the fastest of these rings. Like the other
approximate algorithms they control no errors: a result is accurate
when the problem is well conditioned, and it can be turned into
rigorous enclosures by the certification methods of `acb_mat`. More
eigensolvers (for Hermitian matrices, or multishift variants for large
matrices) may be added here.

::: summary algorithms/linear-algebra/qr-eigenvalues :::

## Eigenvalues and eigenvectors {#eig-qr}

`gr_mat_eig_qr(E, L, R, A, tol, maxiter, ctx)` is the function to call.
It reduces $A$ to upper Hessenberg form, runs the shifted QR algorithm
to a Schur form $A = Q T Q^*$, reads the eigenvalues off the diagonal of
$T$ and computes eigenvectors of $T$ by back substitution, which are
multiplied by $Q$. The eigenvalues come in the order in which they
deflate, not sorted. The right eigenvectors (the columns of `R`) and
the left eigenvectors (the rows of `L`) are scaled so that their
largest entries are of moderate size but are not normalised: for
distinct eigenvalues $LR$ is a diagonal matrix, and $R^{-1} = D^{-1} L$
for its diagonal $D$. Computing eigenvectors costs about as much again
as the eigenvalues alone, so pass `NULL` for those that are not needed.

The ring must support complex conjugation, square roots, division and
comparisons of absolute values (`gr_cmpabs`), and report a working
precision $p$ with `gr_ctx_get_real_prec`. In a real ring the algorithm
works with real shifts only, so it converges when all eigenvalues are
real (for example for a symmetric matrix) and fails otherwise
(typically with `GR_UNABLE`); convert a real matrix to the corresponding complex ring to
get complex eigenvalues. Over ball rings such as `acb` the functions
run, but the radii they produce are meaningless.

The iteration stops when every subdiagonal entry is negligible: smaller
than `tol` times the larger of the two neighbouring diagonal entries, or
times the norm of $A$ divided by $n$ when these are smaller still. With `tol` equal to `NULL`, the tolerance is
$2^{-p}/(100n)$. If some eigenvalue has not converged after `maxiter`
QR steps since the previous one ($14n + p/10$ if `maxiter` is not
positive), the function returns `GR_UNABLE` and writes nothing.

```c example name=gr_mat_eig_qr_rings
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_vec.h>
#include <flint/gr_mat.h>
#include <flint/acb.h>
#include <flint/nfloat.h>

/* the same code over any complex floating-point ring */
static void
eig_demo(const char * name, gr_ctx_t ctx)
{
    slong n = 4, i, j, sz = ctx->sizeof_elem;
    gr_mat_t A, R, AR;
    gr_ptr E, t, res, eps;
    acb_t z;
    gr_ctx_t CC;
    int cmp;

    gr_mat_init(A, n, n, ctx);
    gr_mat_init(R, n, n, ctx);
    gr_mat_init(AR, n, n, ctx);
    E = gr_heap_init_vec(n, ctx);
    t = gr_heap_init(ctx);
    res = gr_heap_init(ctx);
    eps = gr_heap_init(ctx);

    GR_MUST_SUCCEED(gr_mat_set_str(A,
        "[[1, 2, 0, 0], [-2, 1, 0, 0], [1, 0, 3, 1], [0, 1, 0, 2]]", 0, ctx));

    /* eigenvalues and right eigenvectors, default tolerance */
    GR_MUST_SUCCEED(gr_mat_eig_qr(E, NULL, R, A, NULL, 0, ctx));

    /* print the eigenvalues to 10 digits through acb_t */
    gr_ctx_init_complex_acb(CC, 64);
    acb_init(z);
    flint_printf("%s:\n", name);
    for (i = 0; i < n; i++)
    {
        GR_MUST_SUCCEED(gr_set_other(z, GR_ENTRY(E, i, sz), ctx, CC));
        flint_printf("  ");
        acb_printn(z, 10, ARB_STR_NO_RADIUS);
        flint_printf("\n");
    }
    acb_clear(z);
    gr_ctx_clear(CC);

    /* sum of |A R - R diag(E)| over the entries, against 2^-100 */
    GR_MUST_SUCCEED(gr_mat_mul(AR, A, R, ctx));
    GR_MUST_SUCCEED(gr_zero(res, ctx));
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            GR_MUST_SUCCEED(gr_mul(t, gr_mat_entry_srcptr(R, i, j, ctx), GR_ENTRY(E, j, sz), ctx));
            GR_MUST_SUCCEED(gr_sub(t, gr_mat_entry_srcptr(AR, i, j, ctx), t, ctx));
            GR_MUST_SUCCEED(gr_abs(t, t, ctx));
            GR_MUST_SUCCEED(gr_add(res, res, t, ctx));
        }
    GR_MUST_SUCCEED(gr_one(eps, ctx));
    GR_MUST_SUCCEED(gr_mul_2exp_si(eps, eps, -100, ctx));
    GR_MUST_SUCCEED(gr_cmp(&cmp, res, eps, ctx));
    flint_printf("  residual below 2^-100: %s\n", cmp < 0 ? "yes" : "no");

    gr_mat_clear(A, ctx);
    gr_mat_clear(R, ctx);
    gr_mat_clear(AR, ctx);
    gr_heap_clear_vec(E, n, ctx);
    gr_heap_clear(t, ctx);
    gr_heap_clear(res, ctx);
    gr_heap_clear(eps, ctx);
}

int main(void)
{
    gr_ctx_t ctx;

    gr_ctx_init_complex_float_acf(ctx, 128);
    eig_demo("acf, 128 bits", ctx);
    gr_ctx_clear(ctx);

    GR_MUST_SUCCEED(nfloat_complex_ctx_init(ctx, 128, 0));
    eig_demo("nfloat_complex, 128 bits", ctx);
    gr_ctx_clear(ctx);

    flint_cleanup();
    return 0;
}
```

```output
acf, 128 bits:
  3.000000000
  1.000000000 - 2.000000000*I
  1.000000000 + 2.000000000*I
  2.000000000
  residual below 2^-100: yes
nfloat_complex, 128 bits:
  3.000000000
  1.000000000 - 2.000000000*I
  1.000000000 + 2.000000000*I
  2.000000000
  residual below 2^-100: yes
```

With `gr_ctx_init_dfloat(ctx, 3, DFLOAT_COMPLEX)` the same function
gives the eigenvalues in triple-double precision (159 bits), several
times faster than `nfloat_complex` at 192 bits ([](dfloat/performance#linear-algebra)).

## Building blocks {#blocks}

The steps of `gr_mat_eig_qr` are available separately, for computing
a Schur form or for building other eigensolvers.
`gr_mat_hessenberg_reduce` overwrites $A$ with its Hessenberg form,
keeping the Householder reflectors below the subdiagonal and in `tau`
(as LAPACK's `zgehd2` does, but stored by rows); the caller zeroes the
entries below the subdiagonal to get a clean Hessenberg matrix, after
copying them for `gr_mat_hessenberg_expand_q` if $Q$ is needed.
`gr_mat_hessenberg_qr` then iterates `gr_mat_hessenberg_qr_step` with
Wilkinson shifts (and exceptional shifts after 10, 20 and 29 steps
without deflation) and deflation, accumulating the plane rotations into
`Q`; on success $A$ is upper triangular, and if `Q` started as the
matrix of the reduction, the original matrix is $Q T Q^*$.
`gr_mat_eig_triu_r` and `gr_mat_eig_triu_l` compute eigenvectors of
the triangular factor $T$, which `gr_mat_eig_qr` multiplies by $Q$ and
$Q^*$.

`gr_mat_schur_reorder` moves an eigenvalue along the diagonal of a
Schur form by a sequence of plane rotations, each exchanging two
adjacent diagonal entries (like LAPACK's `ztrexc`); applied repeatedly,
it groups selected eigenvalues at the top left, so that the first
columns of the updated $V$ span their invariant subspace. The diagonal
entries are permuted exactly and the part below the diagonal stays
zero.

The plane rotations are applied with `gr_mat_apply_2x2_rows` and
`gr_mat_apply_2x2_cols`, which call the ring method
`_gr_vec_apply_2x2`; the rotations have real cosines, which the complex
`dfloat` rings exploit to save a quarter of the multiplications.

## Certified eigenvalues {#certified}

For rigorous results, convert the matrix to `acb_mat_t` and use the
approximations as the input of `acb_mat_eig_simple` (distinct
eigenvalues) or `acb_mat_eig_multiple` (clusters), which certify them
with error bounds, or call `gr_mat_diagonalization` over an `acb` ring,
which does both
([](algorithms/ball-arithmetic/eigenvalues)). `acb_mat_approx_eig_qr`
itself converts its input to the fastest of these rings that holds the
working precision and the exponents of the entries
([](algorithms/ball-arithmetic/linear-systems#backends)).
