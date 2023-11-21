/*
    Copyright (C) 2008, Martin Albrecht
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010, 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

#include "fmpz_mat.h"

/* todo: optimize for small matrices */
/* todo: bodrato squaring */
/* todo: use fused add-mul operations when supported by
         the matrix interface in the future */

int gr_mat_mul_strassen(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    slong ar, ac, br, bc;
    slong anr, anc, bnr, bnc;
    int status = GR_SUCCESS;

    gr_mat_t A11, A12, A21, A22;
    gr_mat_t B11, B12, B21, B22;
    gr_mat_t C11, C12, C21, C22;
    gr_mat_t X1, X2;

    ar = A->r;
    ac = A->c;
    br = B->r;
    bc = B->c;

    if (ar <= 1 || ac <= 1 || bc <= 1)
    {
        return gr_mat_mul_classical(C, A, B, ctx);
    }

    if (ac != br || ar != C->r || bc != C->c)
    {
        return GR_DOMAIN;
    }

    if (A == C || B == C)
    {
        gr_mat_t T;
        gr_mat_init(T, ar, bc, ctx);
        status |= gr_mat_mul_strassen(T, A, B, ctx);
        status |= gr_mat_swap_entrywise(T, C, ctx);
        gr_mat_clear(T, ctx);
        return status;
    }

    anr = ar / 2;
    anc = ac / 2;
    bnr = anc;
    bnc = bc / 2;

    gr_mat_window_init(A11, A, 0, 0, anr, anc, ctx);
    gr_mat_window_init(A12, A, 0, anc, anr, 2 * anc, ctx);
    gr_mat_window_init(A21, A, anr, 0, 2 * anr, anc, ctx);
    gr_mat_window_init(A22, A, anr, anc, 2 * anr, 2 * anc, ctx);

    gr_mat_window_init(B11, B, 0, 0, bnr, bnc, ctx);
    gr_mat_window_init(B12, B, 0, bnc, bnr, 2 * bnc, ctx);
    gr_mat_window_init(B21, B, bnr, 0, 2 * bnr, bnc, ctx);
    gr_mat_window_init(B22, B, bnr, bnc, 2 * bnr, 2 * bnc, ctx);

    gr_mat_window_init(C11, C, 0, 0, anr, bnc, ctx);
    gr_mat_window_init(C12, C, 0, bnc, anr, 2 * bnc, ctx);
    gr_mat_window_init(C21, C, anr, 0, 2 * anr, bnc, ctx);
    gr_mat_window_init(C22, C, anr, bnc, 2 * anr, 2 * bnc, ctx);

    gr_mat_init(X1, anr, FLINT_MAX(bnc, anc), ctx);
    gr_mat_init(X2, anc, bnc, ctx);

    X1->c = anc;

    status |= gr_mat_sub(X1, A11, A21, ctx);
    status |= gr_mat_sub(X2, B22, B12, ctx);
    status |= gr_mat_mul(C21, X1, X2, ctx);

    status |= gr_mat_add(X1, A21, A22, ctx);
    status |= gr_mat_sub(X2, B12, B11, ctx);
    status |= gr_mat_mul(C22, X1, X2, ctx);

    status |= gr_mat_sub(X1, X1, A11, ctx);
    status |= gr_mat_sub(X2, B22, X2, ctx);
    status |= gr_mat_mul(C12, X1, X2, ctx);

    status |= gr_mat_sub(X1, A12, X1, ctx);
    status |= gr_mat_mul(C11, X1, B22, ctx);

    X1->c = bnc;
    status |= gr_mat_mul(X1, A11, B11, ctx);
    status |= gr_mat_add(C12, X1, C12, ctx);
    status |= gr_mat_add(C21, C12, C21, ctx);
    status |= gr_mat_add(C12, C12, C22, ctx);
    status |= gr_mat_add(C22, C21, C22, ctx);
    status |= gr_mat_add(C12, C12, C11, ctx);
    status |= gr_mat_sub(X2, X2, B21, ctx);
    status |= gr_mat_mul(C11, A22, X2, ctx);

    gr_mat_clear(X2, ctx);

    status |= gr_mat_sub(C21, C21, C11, ctx);
    status |= gr_mat_mul(C11, A12, B21, ctx);

    status |= gr_mat_add(C11, X1, C11, ctx);

    X1->c = FLINT_MAX(bnc, anc);
    gr_mat_clear(X1, ctx);

    gr_mat_window_clear(A11, ctx);
    gr_mat_window_clear(A12, ctx);
    gr_mat_window_clear(A21, ctx);
    gr_mat_window_clear(A22, ctx);

    gr_mat_window_clear(B11, ctx);
    gr_mat_window_clear(B12, ctx);
    gr_mat_window_clear(B21, ctx);
    gr_mat_window_clear(B22, ctx);

    gr_mat_window_clear(C11, ctx);
    gr_mat_window_clear(C12, ctx);
    gr_mat_window_clear(C21, ctx);
    gr_mat_window_clear(C22, ctx);

    if (bc > 2 * bnc)
    {
        gr_mat_t Bc, Cc;
        gr_mat_window_init(Bc, B, 0, 2 * bnc, ac, bc, ctx);
        gr_mat_window_init(Cc, C, 0, 2 * bnc, ar, bc, ctx);
        status |= gr_mat_mul(Cc, A, Bc, ctx);
        gr_mat_window_clear(Bc, ctx);
        gr_mat_window_clear(Cc, ctx);
    }

    if (ar > 2 * anr)
    {
        gr_mat_t Ar, Cr;
        gr_mat_window_init(Ar, A, 2 * anr, 0, ar, ac, ctx);
        gr_mat_window_init(Cr, C, 2 * anr, 0, ar, bc, ctx);
        status |= gr_mat_mul(Cr, Ar, B, ctx);
        gr_mat_window_clear(Ar, ctx);
        gr_mat_window_clear(Cr, ctx);
    }

    if (ac > 2 * anc)
    {
        gr_mat_t Ac, Br, Cb, tmp;
        slong mt, nt;

        gr_mat_window_init(Ac, A, 0, 2 * anc, 2 * anr, ac, ctx);
        gr_mat_window_init(Br, B, 2 * bnr, 0, ac, 2 * bnc, ctx);
        gr_mat_window_init(Cb, C, 0, 0, 2 * anr, 2 * bnc, ctx);

        mt = Ac->r;
        nt = Br->c;

        gr_mat_init(tmp, mt, nt, ctx);
        status |= gr_mat_mul(tmp, Ac, Br, ctx);
        status |= gr_mat_add(Cb, Cb, tmp, ctx);
        gr_mat_clear(tmp, ctx);
        gr_mat_window_clear(Ac, ctx);
        gr_mat_window_clear(Br, ctx);
        gr_mat_window_clear(Cb, ctx);
    }

    return status;
}
