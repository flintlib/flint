/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"
#include "gr_mat.h"

void
ca_mat_solve_triu_classical(ca_mat_t X, const ca_mat_t U,
    const ca_mat_t B, int unit, ca_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_ca_from_ref(gr_ctx, GR_CTX_CC_CA, ctx);
    GR_MUST_SUCCEED(gr_mat_nonsingular_solve_triu_classical((gr_mat_struct *) X, (const gr_mat_struct *) U, (const gr_mat_struct *) B, unit, gr_ctx));
}

void
ca_mat_solve_triu_recursive(ca_mat_t X,
        const ca_mat_t U, const ca_mat_t B, int unit, ca_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_ca_from_ref(gr_ctx, GR_CTX_CC_CA, ctx);
    GR_MUST_SUCCEED(gr_mat_nonsingular_solve_triu_recursive((gr_mat_struct *) X, (const gr_mat_struct *) U, (const gr_mat_struct *) B, unit, gr_ctx));}

void
ca_mat_solve_triu(ca_mat_t X, const ca_mat_t U,
                                    const ca_mat_t B, int unit, ca_ctx_t ctx)
{
    /* todo: tune thresholds */
    if (B->r < 10 || B->c < 10)
        ca_mat_solve_triu_classical(X, U, B, unit, ctx);
    else
        ca_mat_solve_triu_recursive(X, U, B, unit, ctx);
}
