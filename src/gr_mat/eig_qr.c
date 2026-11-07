/*
    Copyright (C) 2013 Timo Hartmann
    Copyright (C) 2018, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"

/*
    Approximate eigenvalues and eigenvectors of a general square matrix
    by the shifted QR algorithm.

    E must point to a vector of length n, which receives the
    eigenvalues.  If R is not NULL it receives a matrix of right
    eigenvectors (as columns); if L is not NULL it receives a matrix of
    left eigenvectors (as rows).  Either may be NULL.

    It is intended for approximate rings without error tracking (acf,
    nfloat_complex, dfloat).  Over a ball ring such as acb it will run,
    but the enclosures will not be useful, since no attempt is made to
    keep the radii small; use acb_mat_eig_enclosure_rump with these
    approximations instead.

    With converged not NULL, a QR iteration that does not converge is not
    an error: the outputs are computed from its current state (as
    acb_mat_approx_eig_qr does) and *converged is set to 0.
*/
int
_gr_mat_eig_qr(gr_ptr E, gr_mat_t L, gr_mat_t R, const gr_mat_t A, gr_srcptr tol, slong maxiter, int * converged, gr_ctx_t ctx)
{
    slong i, j, n;
    slong sz = ctx->sizeof_elem;
    gr_ptr tau;
    gr_mat_t Acopy, Q;
    int status = GR_SUCCESS;

    n = gr_mat_nrows(A, ctx);

    if (converged != NULL)
        *converged = 1;

    if (n != gr_mat_ncols(A, ctx))
        return GR_DOMAIN;

    if (n == 0)
        return GR_SUCCESS;

    GR_TMP_INIT_VEC(tau, n, ctx);
    gr_mat_init(Acopy, n, n, ctx);

    status |= gr_mat_set(Acopy, A, ctx);
    status |= gr_mat_hessenberg_reduce(Acopy, tau, ctx);

    if (status != GR_SUCCESS)
        goto cleanup1;

    if (L != NULL || R != NULL)
    {
        gr_mat_init(Q, n, n, ctx);
        status |= gr_mat_set(Q, Acopy, ctx);
        status |= gr_mat_hessenberg_expand_q(Q, tau, ctx);
    }

    /* Discard the reflectors, leaving a clean Hessenberg matrix. */
    for (i = 0; i < n; i++)
        for (j = i + 2; j < n; j++)
            status |= gr_zero(GR_MAT_ENTRY(Acopy, j, i, sz), ctx);

    if (status != GR_SUCCESS)
        goto cleanup2;

    status = gr_mat_hessenberg_qr(Acopy, (L != NULL || R != NULL) ? Q : NULL,
                                   tol, maxiter, ctx);

    if (status == GR_UNABLE && converged != NULL)
    {
        *converged = 0;
        status = GR_SUCCESS;
    }

    if (status != GR_SUCCESS)
        goto cleanup2;

    for (i = 0; i < n; i++)
        status |= gr_set(GR_ENTRY(E, i, sz), GR_MAT_ENTRY(Acopy, i, i, sz), ctx);

    if (R != NULL)
    {
        gr_mat_t ER;
        gr_mat_init(ER, n, n, ctx);
        status |= gr_mat_eig_triu_r(ER, Acopy, ctx);
        status |= gr_mat_mul(R, Q, ER, ctx);
        gr_mat_clear(ER, ctx);
    }

    if (L != NULL)
    {
        gr_mat_t EL;
        gr_mat_init(EL, n, n, ctx);
        status |= gr_mat_eig_triu_l(EL, Acopy, ctx);
        status |= gr_mat_conjugate_transpose(Q, Q, ctx);
        status |= gr_mat_mul(L, EL, Q, ctx);
        gr_mat_clear(EL, ctx);
    }

cleanup2:
    if (L != NULL || R != NULL)
        gr_mat_clear(Q, ctx);

cleanup1:
    gr_mat_clear(Acopy, ctx);
    GR_TMP_CLEAR_VEC(tau, n, ctx);

    return status;
}

int
gr_mat_eig_qr(gr_ptr E, gr_mat_t L, gr_mat_t R, const gr_mat_t A, gr_srcptr tol, slong maxiter, gr_ctx_t ctx)
{
    return _gr_mat_eig_qr(E, L, R, A, tol, maxiter, NULL, ctx);
}
