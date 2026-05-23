/*
    Copyright (C) 2026 Maria Neagoie

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "gr.h"
#include "gr_ore_poly.h"
#include "gr_vec.h"

// Returns the unique pair (Q, R) such that U = QV + R and ord(R) < ord(V)
int _gr_ore_poly_divrem(gr_ptr Q, gr_ptr R, gr_srcptr U, slong lenU, gr_srcptr V, slong lenV, gr_ore_poly_ctx_t ctx)
{
    gr_ctx_struct *cctx = GR_ORE_POLY_ELEM_CTX(ctx);
    slong el_size = gr_ctx_sizeof_elem(cctx);
    slong lenQ, lenR, ordV, ordR;
    int status = GR_SUCCESS;

    lenQ = lenU - lenV + 1;
    lenR = lenU;
    ordV = lenV - 1;
    ordR = lenR - 1;

    // Set Q to 0
    status |= _gr_vec_zero(Q, lenQ, cctx);

    // Set R to U
    status |= _gr_vec_set(R, U, lenU, cctx);

    gr_ptr c;
    GR_TMP_INIT(c, cctx);

    gr_ptr A, B;
    GR_TMP_INIT_VEC(A, lenQ, cctx);
    GR_TMP_INIT_VEC(B, lenU, cctx);

    gr_ptr sigma_pows;
    GR_TMP_INIT_VEC(sigma_pows, lenQ, cctx);
    status |= gr_set(GR_ENTRY(sigma_pows, 0, el_size), GR_ENTRY(V, ordV, el_size), cctx);
    for (slong i = 1; i < lenQ; i++)
        status |= gr_ore_poly_sigma(GR_ENTRY(sigma_pows, i, el_size), GR_ENTRY(sigma_pows, i - 1, el_size), ctx);

    while (ordR > ordV)
    {
        slong k = ordR - ordV;

        // c = lc(R) / sigma ^ k (lc(V)) using sigma_pows computed above
        status |= gr_div(c, GR_ENTRY(R, ordR, el_size), GR_ENTRY(sigma_pows, k, el_size), cctx);

        // R -= c * x^k * V. We compute A =  c * x^k, then B = A * V
        // A = c * x^k, so A[k] = c, rest 0
        slong lenA = k + 1;
        status |= _gr_vec_zero(A, lenA, cctx);
        status |= gr_set(GR_ENTRY(A, k, el_size), c, cctx);

        // B = A * V
        status |= _gr_ore_poly_mul(B, A, lenA, V, lenV, ctx);

        // R -= B
        slong lenB = lenA + lenV - 1;
        status |= _gr_vec_sub(R, R, B, lenB, cctx);

        // Q += c * x^k, so Q[k] += c
        status |= gr_add(GR_ENTRY(Q, k, el_size), GR_ENTRY(Q, k, el_size), c, cctx);

        ordR--;
    }

    GR_TMP_CLEAR_VEC(sigma_pows, lenQ, cctx);
    GR_TMP_CLEAR(c, cctx);

    GR_TMP_CLEAR_VEC(A, lenQ, cctx);
    GR_TMP_CLEAR_VEC(B, lenU, cctx);

    return status;
}

int gr_ore_poly_divrem(gr_ore_poly_t Q, gr_ore_poly_t R, const gr_ore_poly_t U, gr_ore_poly_t V, gr_ore_poly_ctx_t ctx)
{
    int status = GR_SUCCESS;
    const slong lenU = U->length;
    const slong lenV = V->length;

    if (lenV == 0) // division by zero polynomial
        return GR_DOMAIN;

    if (lenU == 0) // result is 0
    {
        status |= gr_ore_poly_zero(Q, ctx);
        status |= gr_ore_poly_zero(R, ctx);
        return status;
    }

    if (lenU < lenV)
    {
        status |= gr_ore_poly_set(R, U, ctx);
        status |= gr_ore_poly_zero(Q, ctx);
        return status;
    }

    slong lenQ = lenU - lenV + 1;

    if (Q == U || Q == V || R == U || R == V || Q == R) // treat aliasing case separately
    {
        gr_ore_poly_t tQ, tR;
        gr_ore_poly_init(tQ, ctx);
        gr_ore_poly_init(tR, ctx);

        gr_ore_poly_fit_length(tQ, lenQ, ctx);
        gr_ore_poly_fit_length(tR, lenU, ctx);
        status = _gr_ore_poly_divrem(tQ->coeffs, tR->coeffs, U->coeffs, lenU, V->coeffs, lenV, ctx);
        _gr_ore_poly_set_length(tQ, lenQ, ctx);
        _gr_ore_poly_set_length(tR, lenU, ctx);
        _gr_ore_poly_normalise(tQ, ctx);
        _gr_ore_poly_normalise(tR, ctx);

        gr_ore_poly_swap(Q, tQ, ctx);
        gr_ore_poly_swap(R, tR, ctx);
        gr_ore_poly_clear(tQ, ctx);
        gr_ore_poly_clear(tR, ctx);
    }
    else
    {
        gr_ore_poly_fit_length(Q, lenQ, ctx);
        gr_ore_poly_fit_length(R, lenU, ctx);
        status = _gr_ore_poly_divrem(Q->coeffs, R->coeffs, U->coeffs, lenU, V->coeffs, lenV, ctx);
        _gr_ore_poly_set_length(Q, lenQ, ctx);
        _gr_ore_poly_set_length(R, lenU, ctx);
        _gr_ore_poly_normalise(Q, ctx);
        _gr_ore_poly_normalise(R, ctx);
    }
    return status;
}

int gr_ore_poly_div(gr_ore_poly_t Q, const gr_ore_poly_t U, gr_ore_poly_t V, gr_ore_poly_ctx_t ctx)
{
    gr_ore_poly_t R;
    int status;
    gr_ore_poly_init(R, ctx);
    status = gr_ore_poly_divrem(Q, R, U, V, ctx);
    gr_ore_poly_clear(R, ctx);
    return status;
}

int gr_ore_poly_rem(gr_ore_poly_t R, const gr_ore_poly_t U, gr_ore_poly_t V, gr_ore_poly_ctx_t ctx)
{
    gr_ore_poly_t Q;
    int status;
    gr_ore_poly_init(Q, ctx);
    status = gr_ore_poly_divrem(Q, R, U, V, ctx);
    gr_ore_poly_clear(Q, ctx);
    return status;
}
