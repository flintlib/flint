/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "ulong_extras.h"
#include "gr_vec.h"
#include "gr_poly.h"

ulong
gr_poly_deflation(const gr_poly_t poly, gr_ctx_t ctx)
{
    slong i, sz = ctx->sizeof_elem;
    ulong deflation;

    if (poly->length <= 1)
        return poly->length;

    /* The exponents of the nonzero terms. A term whose coefficient cannot
       be verified to be zero is treated as nonzero, so the deflation is
       always valid (possibly not optimal). */
    deflation = 0;
    for (i = 1; i < poly->length - 1 && deflation != 1; i++)
        if (gr_is_zero(GR_ENTRY(poly->coeffs, i, sz), ctx) != T_TRUE)
            deflation = n_gcd(deflation, i);

    return n_gcd(deflation, poly->length - 1);
}

int
gr_poly_deflate(gr_poly_t res, const gr_poly_t poly, ulong deflation, gr_ctx_t ctx)
{
    slong i, len, sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (deflation == 0)
        return GR_DOMAIN;

    if (poly->length <= 1 || deflation == 1)
        return gr_poly_set(res, poly, ctx);

    len = (poly->length - 1) / deflation + 1;

    if (res == poly)
    {
        for (i = 1; i < len; i++)
            gr_swap(GR_ENTRY(res->coeffs, i, sz), GR_ENTRY(res->coeffs, i * deflation, sz), ctx);
    }
    else
    {
        gr_method_unary_op set = GR_UNARY_OP(ctx, SET);
        gr_ptr r;
        gr_srcptr p;

        gr_poly_fit_length(res, len, ctx);
        r = res->coeffs;
        p = poly->coeffs;
        for (i = 0; i < len; i++)
            status |= set(GR_ENTRY(r, i, sz), GR_ENTRY(p, i * deflation, sz), ctx);
    }

    _gr_poly_set_length_normalise(res, len, ctx);

    return status;
}

FLINT_FORCE_INLINE void
_gr_poly_swap_entries_bytewise(char * a, char * b, slong sz)
{
    if (sz == sizeof(ulong))
    {
        ulong t;
        memcpy(&t, a, sizeof(ulong));
        memcpy(a, b, sizeof(ulong));
        memcpy(b, &t, sizeof(ulong));
    }
    else
    {
        slong k;
        char t;

        for (k = 0; k < sz; k++)
        {
            t = a[k];
            a[k] = b[k];
            b[k] = t;
        }
    }
}

/* poly(x^0) = poly(1); kept out of line so that the main path of
   gr_poly_inflate needs no stack temporaries */
static int
_gr_poly_inflate_zero(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx)
{
    gr_ptr t;
    int status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);
    status |= gr_one(t, ctx);
    status |= gr_poly_evaluate(t, poly, t, ctx);
    status |= gr_poly_set_scalar(res, t, ctx);
    GR_TMP_CLEAR(t, ctx);

    return status;
}

int
gr_poly_inflate(gr_poly_t res, const gr_poly_t poly, ulong inflation, gr_ctx_t ctx)
{
    slong len;
    int status = GR_SUCCESS;

    if (inflation == 0)
        return _gr_poly_inflate_zero(res, poly, ctx);

    if (poly->length <= 1 || inflation == 1)
        return gr_poly_set(res, poly, ctx);

    len = (poly->length - 1) * inflation + 1;

    if (res == poly)
    {
        slong i, sz = ctx->sizeof_elem, plen = poly->length;
        char * r;

        gr_poly_fit_length(res, len, ctx);
        r = res->coeffs;

        /* With zeros after the coefficients, the transpositions
           (i, i * inflation) for decreasing i move the coefficients to
           their final positions and leave zeros in the gaps (this is the
           second pass of the non-aliased case below). */
        status |= _gr_vec_zero(r + plen * sz, len - plen, ctx);
        for (i = plen - 1; i >= 1; i--)
            _gr_poly_swap_entries_bytewise(r + i * sz, r + i * inflation * sz, sz);
    }
    else if (poly->length < 32)
    {
        gr_method_unary_op set = GR_UNARY_OP(ctx, SET);
        gr_method_vec_constant_op vec_zero = GR_VEC_CONSTANT_OP(ctx, VEC_ZERO);
        slong i, sz = ctx->sizeof_elem;
        gr_ptr r;
        gr_srcptr p;

        /* write each coefficient directly to its final position (reusing
           any memory already allocated for the entries of res) */
        gr_poly_fit_length(res, len, ctx);
        r = res->coeffs;
        p = poly->coeffs;
        status |= set(r, p, ctx);
        for (i = 1; i < poly->length; i++)
        {
            status |= vec_zero(GR_ENTRY(r, (i - 1) * inflation + 1, sz), inflation - 1, ctx);
            status |= set(GR_ENTRY(r, i * inflation, sz), GR_ENTRY(p, i, sz), ctx);
        }
    }
    else
    {
        slong i, sz = ctx->sizeof_elem, plen = poly->length;
        char * r;

        gr_poly_fit_length(res, len, ctx);
        r = res->coeffs;

        /* Same result with O(1) instead of O(len) method calls: move the
           entries at the target positions i * inflation to the front (the
           transpositions (i, i * inflation) for increasing i), so that the
           coefficients can be copied with a single vector operation
           (reusing any memory already allocated for these entries) and the
           gaps can be zeroed with a single vector operation; then undo the
           permutation. Elements are relocatable, so they can be moved
           bytewise. */
        for (i = 1; i < plen; i++)
            _gr_poly_swap_entries_bytewise(r + i * sz, r + i * inflation * sz, sz);

        status |= _gr_vec_set(r, poly->coeffs, plen, ctx);
        status |= _gr_vec_zero(r + plen * sz, len - plen, ctx);

        for (i = plen - 1; i >= 1; i--)
            _gr_poly_swap_entries_bytewise(r + i * sz, r + i * inflation * sz, sz);
    }

    _gr_poly_set_length_normalise(res, len, ctx);

    return status;
}
