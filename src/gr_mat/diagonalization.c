/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_mat.h"

int
gr_mat_diagonalization_precomp(gr_vec_t D, gr_mat_t L, gr_mat_t R, const gr_mat_t A, const gr_vec_t eigenvalues, const gr_vec_t mult, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_mat_t AIe, b;
    const fmpz * am = mult->entries;
    slong i, j, k, n;
    slong nullity, added, deg;

    if (gr_mat_is_square(A, ctx) != T_TRUE)
        return GR_DOMAIN;

    n = gr_mat_nrows(A, ctx);

    if (L != NULL && R == NULL)
    {
        gr_mat_t tmpR;
        gr_mat_init(tmpR, n, n, ctx);
        status = gr_mat_diagonalization_precomp(D, L, tmpR, A, eigenvalues, mult, ctx);
        gr_mat_clear(tmpR, ctx);
        return status;
    }

    gr_vec_set_length(D, n, ctx);
    /* todo: resize L, R if needed */


    deg = 0;
    for (i = 0; i < eigenvalues->length; i++)
        deg += am[i];

    if (deg != n)
        return GR_DOMAIN;

    gr_mat_init(AIe, n, n, ctx);
    gr_mat_init(b, 0, 0, ctx);

    status |= _gr_vec_zero(D->entries, n, ctx);

    added = 0;
    for (i = 0; i < eigenvalues->length; i++)
    {
        status |= gr_mat_set(AIe, A, ctx);
        for (j = 0; j < n; j++)
            status |= gr_sub(gr_mat_entry_ptr(AIe, j, j, ctx), gr_mat_entry_ptr(AIe, j, j, ctx), gr_vec_entry_ptr((gr_vec_struct *) eigenvalues, i, ctx), ctx);

        status |= gr_mat_nullspace(b, AIe, ctx);
        if (status != GR_SUCCESS)
        {
            status = GR_UNABLE;
            break;
        }

        nullity = gr_mat_ncols(b, ctx);

        if (nullity != am[i])
        {
            status = GR_DOMAIN;
            break;
        }

        for (j = 0; j < am[i]; j++)
        {
            status |= gr_set(gr_vec_entry_ptr(D, added + j, ctx), gr_vec_entry_ptr((gr_vec_struct *) eigenvalues, i, ctx), ctx);

            if (R != NULL)
            {
                for (k = 0; k < n; k++)
                    status |= gr_set(gr_mat_entry_ptr(R, k, added + j, ctx), gr_mat_entry_ptr(b, k, j, ctx), ctx);
            }
        }

        added += am[i];
    }

    if (status == GR_SUCCESS && L != NULL)
        status = gr_mat_inv(L, R, ctx);

    gr_mat_clear(AIe, ctx);
    gr_mat_clear(b, ctx);

    return status;
}

int
gr_mat_diagonalization_generic(gr_vec_t D, gr_mat_t L, gr_mat_t R, const gr_mat_t A, int flags, gr_ctx_t ctx)
{
    int status;
    gr_ctx_t ZZ;
    gr_vec_t eigenvalues, mult;

    if (gr_mat_is_square(A, ctx) != T_TRUE)
        return GR_DOMAIN;

    gr_ctx_init_fmpz(ZZ);
    gr_vec_init(eigenvalues, 0, ctx);
    gr_vec_init(mult, 0, ZZ);

    if (gr_mat_eigenvalues(eigenvalues, mult, A, flags, ctx) == GR_SUCCESS)
    {
        status = gr_mat_diagonalization_precomp(D, L, R, A, eigenvalues, mult, ctx);
    }
    else
    {
        status = GR_UNABLE;
    }

    gr_vec_clear(eigenvalues, ctx);
    gr_vec_clear(mult, ZZ);
    gr_ctx_clear(ZZ);

    return status;
}

int
gr_mat_diagonalization(gr_vec_t D, gr_mat_t L, gr_mat_t R, const gr_mat_t A, int flags, gr_ctx_t ctx)
{
    return GR_MAT_DIAGONALIZATION_OP(ctx, MAT_DIAGONALIZATION)(D, L, R, A, flags, ctx);
}
