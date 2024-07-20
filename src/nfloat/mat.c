/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "nfloat.h"
#include "gr_mat.h"

int
nfloat_mat_nonsingular_solve_tril(gr_mat_t X, const gr_mat_t L,
                                    const gr_mat_t B, int unit, gr_ctx_t ctx)
{
    slong cutoff, prec = NFLOAT_CTX_PREC(ctx);

    if (prec <= 256)
        cutoff = 96;
    else if (prec <= 512)
        cutoff = 16;
    else if (prec <= 576)
        cutoff = 32;
    else if (prec <= 1536)
        cutoff = 8;
    else if (prec <= 2176)
        cutoff = 7;
    else
        cutoff = 6;

    if (B->r < cutoff || B->c < cutoff)
        return gr_mat_nonsingular_solve_tril_classical(X, L, B, unit, ctx);
    else
        return gr_mat_nonsingular_solve_tril_recursive(X, L, B, unit, ctx);
}

int
nfloat_mat_nonsingular_solve_triu(gr_mat_t X, const gr_mat_t L,
                                    const gr_mat_t B, int unit, gr_ctx_t ctx)
{
    slong cutoff, prec = NFLOAT_CTX_PREC(ctx);

    if (prec <= 256)
        cutoff = 96;
    else if (prec <= 512)
        cutoff = 16;
    else if (prec <= 576)
        cutoff = 32;
    else if (prec <= 1536)
        cutoff = 8;
    else if (prec <= 2176)
        cutoff = 7;
    else
        cutoff = 6;

    if (B->r < cutoff || B->c < cutoff)
        return gr_mat_nonsingular_solve_triu_classical(X, L, B, unit, ctx);
    else
        return gr_mat_nonsingular_solve_triu_recursive(X, L, B, unit, ctx);
}

int
nfloat_complex_mat_nonsingular_solve_tril(gr_mat_t X, const gr_mat_t L,
                                    const gr_mat_t B, int unit, gr_ctx_t ctx)
{
    slong cutoff, prec = NFLOAT_CTX_PREC(ctx);

    if (prec <= 192)
        cutoff = 64;
    else if (prec <= 256)
        cutoff = 16;
    else if (prec <= 384)
        cutoff = 7;
    else if (prec == 576)
        cutoff = 16;
    else
        cutoff = 6;

    if (B->r < cutoff || B->c < cutoff)
        return gr_mat_nonsingular_solve_tril_classical(X, L, B, unit, ctx);
    else
        return gr_mat_nonsingular_solve_tril_recursive(X, L, B, unit, ctx);
}

int
nfloat_complex_mat_nonsingular_solve_triu(gr_mat_t X, const gr_mat_t L,
                                    const gr_mat_t B, int unit, gr_ctx_t ctx)
{
    slong cutoff, prec = NFLOAT_CTX_PREC(ctx);

    if (prec <= 192)
        cutoff = 64;
    else if (prec <= 256)
        cutoff = 16;
    else if (prec <= 384)
        cutoff = 7;
    else if (prec == 576)
        cutoff = 16;
    else
        cutoff = 6;

    if (B->r < cutoff || B->c < cutoff)
        return gr_mat_nonsingular_solve_triu_classical(X, L, B, unit, ctx);
    else
        return gr_mat_nonsingular_solve_triu_recursive(X, L, B, unit, ctx);
}

int
nfloat_mat_lu(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)
{
    slong cutoff, prec = NFLOAT_CTX_PREC(ctx);

    if (prec <= 256)
        cutoff = 32;
    else if (prec <= 576)
        cutoff = 28;
    else if (prec <= 768)
        cutoff = 16;
    else if (prec <= 1536)
        cutoff = 12;
    else if (prec <= 2560)
        cutoff = 8;
    else
        cutoff = 7;

    if (A->r < cutoff || A->c < cutoff)
        return gr_mat_lu_classical(rank, P, LU, A, rank_check, ctx);
    else
        return gr_mat_lu_recursive(rank, P, LU, A, rank_check, ctx);
}

int
nfloat_complex_mat_lu(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)
{
    slong cutoff, prec = NFLOAT_CTX_PREC(ctx);

    if (prec <= 256)
        cutoff = 12;
    else if (prec <= 512)
        cutoff = 8;
    else if (prec <= 576)
        cutoff = 16;
    else if (prec <= 1024)
        cutoff = 7;
    else
        cutoff = 6;

    if (A->r < cutoff || A->c < cutoff)
        return gr_mat_lu_classical(rank, P, LU, A, rank_check, ctx);
    else
        return gr_mat_lu_recursive(rank, P, LU, A, rank_check, ctx);
}
