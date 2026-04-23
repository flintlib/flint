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

// Implements left multiplication: D * Q = sigma(Q) * D + delta(Q)
int _gr_ore_poly_lmul_gen(gr_ptr res, gr_srcptr poly, slong len, gr_ore_poly_ctx_t ctx)
{
    slong el_size = gr_ctx_sizeof_elem(GR_ORE_POLY_ELEM_CTX(ctx));
    int status = GR_SUCCESS;

    // Set the new coefficient (L) to zero
    status |= gr_zero(GR_ENTRY(res, len, el_size), GR_ORE_POLY_ELEM_CTX(ctx));

    gr_ptr temps;
    gr_ptr tempd;
    GR_TMP_INIT(temps, GR_ORE_POLY_ELEM_CTX(ctx));
    GR_TMP_INIT(tempd, GR_ORE_POLY_ELEM_CTX(ctx));

    for (slong j = len - 1; j >= 0; j--)
    {
        gr_srcptr qj = GR_ENTRY(poly, j, el_size);

        // Compute both sigma and delta at once
        status |= gr_ore_poly_sigma_delta(temps, tempd, qj, ctx);

        // res[j+1] += sigma(qj)
        status |= gr_add(GR_ENTRY(res, j+1, el_size), GR_ENTRY(res, j+1, el_size), temps, GR_ORE_POLY_ELEM_CTX(ctx));

        // res[j] = delta(qj) (overwrite)
        status |= gr_set(GR_ENTRY(res, j, el_size), tempd, GR_ORE_POLY_ELEM_CTX(ctx));
    }

    GR_TMP_CLEAR(temps, GR_ORE_POLY_ELEM_CTX(ctx));
    GR_TMP_CLEAR(tempd, GR_ORE_POLY_ELEM_CTX(ctx));

    return status;
}

int gr_ore_poly_lmul_gen(gr_ore_poly_t res, const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
{
    const slong len = poly->length;
    if (len == 0)
        return gr_ore_poly_zero(res, ctx);

    // By multiplication, length can increase by 1
    gr_ore_poly_fit_length(res, len + 1, ctx);
    int status = _gr_ore_poly_lmul_gen(res->coeffs, poly->coeffs, len, ctx);
    _gr_ore_poly_set_length(res, len + 1, ctx);
    _gr_ore_poly_normalise(res, ctx);
    return status;
}

// Assumes: res != a, res != b (no aliasing), a->length > 0, b->length > 0
int _gr_ore_poly_mul(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ore_poly_ctx_t ctx)
{
    slong el_size = gr_ctx_sizeof_elem(GR_ORE_POLY_ELEM_CTX(ctx));
    int status = GR_SUCCESS;

    const slong lres = len1 + len2 - 1;

    // Output initialized at 0
    for (slong k = 0; k < lres; k++)
        status |= gr_zero(GR_ENTRY(res, k, el_size), GR_ORE_POLY_ELEM_CTX(ctx));

    // Initialize Q = b, then Q = D^i * b.
    gr_ore_poly_t Q;
    gr_ore_poly_init2(Q, len2, ctx);
    for (slong i = 0; i < len2; i++)
        status |= gr_set(GR_ENTRY(Q->coeffs, i, el_size), GR_ENTRY(poly2, i, el_size), GR_ORE_POLY_ELEM_CTX(ctx));
    _gr_ore_poly_set_length(Q, len2, ctx);
    _gr_ore_poly_normalise(Q, ctx);

    gr_ptr temp_prod;
    GR_TMP_INIT(temp_prod, GR_ORE_POLY_ELEM_CTX(ctx));

    for (slong i = 0; i < len1; i++)
    {
        gr_srcptr ai = GR_ENTRY(poly1, i, el_size);

        const slong LQ = Q->length;
        for (slong j = 0; j < LQ; j++)
        {
            gr_srcptr qj = GR_ENTRY(Q->coeffs, j, el_size);

            // Compute a_i * q_j
            status |= gr_mul(temp_prod, ai, qj, GR_ORE_POLY_ELEM_CTX(ctx));

            // Add it at degree j of res, since Q = D^i * b
            status |= gr_add(GR_ENTRY(res, j, el_size), GR_ENTRY(res, j, el_size), temp_prod, GR_ORE_POLY_ELEM_CTX(ctx));
        }

        // Update Q = Q * D so that Q = D^i * b is also true at the next step
        if (i + 1 < len1)
        {
            status |= gr_ore_poly_lmul_gen(Q, Q, ctx);
        }
    }

    GR_TMP_CLEAR(temp_prod, GR_ORE_POLY_ELEM_CTX(ctx));

    gr_ore_poly_clear(Q, ctx);

    return status;
}

// High-level multiplication of a(D) * b(D)
int gr_ore_poly_mul(gr_ore_poly_t res, const gr_ore_poly_t poly1, const gr_ore_poly_t poly2, gr_ore_poly_ctx_t ctx)
{
    const slong len1 = poly1->length;
    const slong len2 = poly2->length;

    if (len1 == 0 || len2 == 0)
        return gr_ore_poly_zero(res, ctx);

    if (len1 + len2 - 1 > ctx->size_limit)
        return GR_UNABLE;

    // Aliasing
    if (res == poly1 || res == poly2)
    {
        gr_ore_poly_t t;
        gr_ore_poly_init2(t, len1 + len2 - 1, ctx);
        int status = _gr_ore_poly_mul(t->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, ctx);
        _gr_ore_poly_set_length(t, len1 + len2 - 1, ctx);
        _gr_ore_poly_normalise(t, ctx);
        gr_ore_poly_swap(res, t, ctx);
        gr_ore_poly_clear(t, ctx);
        return status;
    }

    gr_ore_poly_fit_length(res, len1 + len2 - 1, ctx);
    int status = _gr_ore_poly_mul(res->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, ctx);
    _gr_ore_poly_set_length(res, len1 + len2 - 1, ctx);
    _gr_ore_poly_normalise(res, ctx);
    return status;
}
