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

// Returns the unique pair (Q, R) such that U = QV + R and ord(R) < ord(V)
int _gr_ore_poly_divrem(gr_ptr Q, gr_ptr R, gr_srcptr U, slong lenU, gr_srcptr V, slong lenV, gr_ore_poly_ctx_t ctx)
{
    gr_ctx_struct *cctx = GR_ORE_POLY_ELEM_CTX(ctx);
    slong el_size = gr_ctx_sizeof_elem(cctx);
    slong lenQ, lenR, ordV, ordR;
    int status = GR_SUCCESS;

    if (GR_ORE_POLY_CTX(ctx)->sigma_delta == NULL)
        return GR_UNABLE;

    lenQ = lenU - lenV + 1;
    lenR = lenU;
    ordV = lenV - 1;
    ordR = lenR - 1;

    // Set Q to 0
    for (slong k = 0; k < lenQ; k++)
        status |= gr_zero(GR_ENTRY(Q, k, el_size), cctx);

    // Set R to U
    for (slong k = 0; k < lenU; k++)
        status |= gr_set(GR_ENTRY(R, k, el_size), GR_ENTRY(U, k, el_size), cctx);

    gr_ptr lcR, lcV, denominator, c;
    GR_TMP_INIT(lcR, cctx);
    GR_TMP_INIT(lcV, cctx);
    GR_TMP_INIT(denominator, cctx);
    GR_TMP_INIT(c, cctx);

    gr_ptr A = flint_malloc(lenQ * el_size);
    gr_ptr B = flint_malloc(lenU * el_size);
    _gr_vec_init(A, lenQ, cctx);
    _gr_vec_init(B, lenU, cctx);

    while (ordR > ordV)
    {
        slong k = ordR - ordV;

        status |= gr_set(lcR, GR_ENTRY(R, ordR, el_size), cctx);
        status |= gr_set(lcV, GR_ENTRY(V, ordV, el_size), cctx);

        // Compute denominator = sigma ^ k (lc(V))
        status |= gr_set(denominator, lcV, cctx);
        for (slong i = 0; i < k; i++)
            status |= gr_ore_poly_sigma(denominator, denominator, ctx);

        // c = lc(R) / denominator
        status |= gr_div(c, lcR, denominator, cctx);

        // R -= c * x^k * V. We compute A =  c * x^k, then B = A * V
        // A = c * x^k, so A[k] = c, rest 0
        slong lenA = k + 1;
        for (slong i = 0; i < lenA; i++)
            status |= gr_zero(GR_ENTRY(A, i, el_size), cctx);
        status |= gr_set(GR_ENTRY(A, k, el_size), c, cctx);

        // B = A * V
        status |= _gr_ore_poly_mul(B, A, lenA, V, lenV, ctx);

        // R -= B
        slong lenB = lenA + lenV - 1;
        for (slong i = 0; i < lenB; i++)
            status |= gr_sub(GR_ENTRY(R, i, el_size), GR_ENTRY(R, i, el_size), GR_ENTRY(B, i, el_size), cctx);

        // Q += c * x^k, so Q[k] += c
        status |= gr_add(GR_ENTRY(Q, k, el_size), GR_ENTRY(Q, k, el_size), c, cctx);

        ordR--;
    }

    gr_clear(lcR, cctx);
    gr_clear(lcV, cctx);
    gr_clear(denominator, cctx);
    gr_clear(c, cctx);

    _gr_vec_clear(A, lenQ, cctx);
    _gr_vec_clear(B, lenU, cctx);
    flint_free(A);
    flint_free(B);

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
        gr_ore_poly_t tU;
        // take care of aliasing case with temp
        gr_ore_poly_init(tU, ctx);
        status |= gr_ore_poly_set(tU, U, ctx);
        status |= gr_ore_poly_zero(Q, ctx);
        status |= gr_ore_poly_set(R, tU, ctx);
        gr_ore_poly_clear(tU, ctx);
        return status;
    }

    slong lenQ = lenU - lenV + 1;

    if (Q == U || Q == V || R == U || R == V) // treat aliasing case separately
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
