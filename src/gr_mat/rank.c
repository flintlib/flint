/*
    Copyright (C) 2022, 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "perm.h"
#include "gr.h"
#include "gr_mat.h"

int
gr_mat_rank_fflu(slong * rank, const gr_mat_t A, gr_ctx_t ctx)
{
    slong n, m;
    slong * P;
    int status;
    gr_mat_t T;
    gr_ptr den;

    n = gr_mat_nrows(A, ctx);
    m = gr_mat_ncols(A, ctx);

    if (n == 0 || m == 0)
    {
        *rank = 0;
        return GR_SUCCESS;
    }
    else
    {
        GR_TMP_INIT(den, ctx);

        gr_mat_init(T, n, m, ctx);
        P = _perm_init(n);

        status = gr_mat_fflu(rank, P, T, den, A, 0, ctx);

        gr_mat_clear(T, ctx);
        _perm_clear(P);

        GR_TMP_CLEAR(den, ctx);

        if (status != GR_SUCCESS)
            status |= GR_UNABLE;

        return status;
    }
}

int
gr_mat_rank_lu(slong * rank, const gr_mat_t A, gr_ctx_t ctx)
{
    slong n, m;
    slong * P;
    int status;
    gr_mat_t T;

    n = gr_mat_nrows(A, ctx);
    m = gr_mat_ncols(A, ctx);

    if (n == 0 || m == 0)
    {
        *rank = 0;
        return GR_SUCCESS;
    }
    else
    {
        gr_mat_init(T, n, m, ctx);
        P = _perm_init(n);

        status = gr_mat_lu(rank, P, T, A, 0, ctx);

        gr_mat_clear(T, ctx);
        _perm_clear(P);

        if (status != GR_SUCCESS)
            status |= GR_UNABLE;

        return status;
    }
}

int
gr_mat_rank(slong * rank, const gr_mat_t A, gr_ctx_t ctx)
{
    truth_t dom;

    /* Sensible definition over any ring. */
    if (gr_mat_nrows(A, ctx) == 0 || gr_mat_ncols(A, ctx) == 0)
    {
        *rank = 0;
        return GR_SUCCESS;
    }

    dom = gr_ctx_is_field(ctx);

    /* Prefer standard LU only over finite fields */
    /* Prefer FFLU over non-finite fields as it often results in
       smaller coefficients. TODO: this choice surely wants
       tuning, e.g. when we have fast matrix multiplication.
       For example, ca_mat currently uses LU for number fields. */
    if (dom == T_TRUE && gr_ctx_is_finite(ctx) == T_TRUE)
        return gr_mat_rank_lu(rank, A, ctx);

    if (dom != T_TRUE)
        dom = gr_ctx_is_integral_domain(ctx);

    if (dom == T_TRUE)
        return gr_mat_rank_fflu(rank, A, ctx);

    /* There are ways to generalize the notation of rank to
       non-integral domains, but we do not currently implement them;
       for now we define the rank as undefined. */
    if (dom == T_FALSE)
        return GR_DOMAIN;

    return GR_UNABLE;
}
