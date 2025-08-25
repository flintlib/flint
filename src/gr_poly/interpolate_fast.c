/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "gr_vec.h"
#include "gr_poly.h"

/* todo */
#define _gr_poly_mul_monic _gr_poly_mul

int
_gr_poly_interpolation_weights(gr_ptr w,
    const gr_ptr * tree, slong len, gr_ctx_t ctx)
{
    gr_ptr tmp;
    slong i, n, height;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (len == 0)
        return status;

    if (len == 1)
        return gr_one(w, ctx);

    GR_TMP_INIT_VEC(tmp, len + 1, ctx);
    height = FLINT_CLOG2(len);
    n = WORD(1) << (height - 1);

    status |= _gr_poly_mul_monic(tmp, tree[height-1], n + 1,
                        GR_ENTRY(tree[height-1], (n + 1), sz), (len - n + 1), ctx);

    status |= _gr_poly_derivative(tmp, tmp, len + 1, ctx);
    status |= _gr_poly_evaluate_vec_fast_precomp(w, tmp, len, tree, len, ctx);

    for (i = 0; i < len; i++)
        status |= gr_inv(GR_ENTRY(w, i, sz), GR_ENTRY(w, i, sz), ctx);

    GR_TMP_CLEAR_VEC(tmp, len + 1, ctx);

    return status;
}

int
_gr_poly_interpolate_fast_precomp(gr_ptr poly,
    gr_srcptr ys, const gr_ptr * tree, gr_srcptr weights,
    slong len, gr_ctx_t ctx)
{
    gr_ptr t, u, pa, pb;
    slong i, pow, left;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (len == 0)
        return status;

    GR_TMP_INIT_VEC(t, 2 * len, ctx);
    u = GR_ENTRY(t, len, sz);

    status |= _gr_vec_mul(poly, weights, ys, len, ctx);

    for (i = 0; i < FLINT_CLOG2(len); i++)
    {
        pow = (WORD(1) << i);
        pa = tree[i];
        pb = poly;
        left = len;

        while (left >= 2 * pow)
        {
            status |= _gr_poly_mul(t, pa, pow + 1, GR_ENTRY(pb, pow, sz), pow, ctx);
            status |= _gr_poly_mul(u, GR_ENTRY(pa, pow + 1, sz), pow + 1, pb, pow, ctx);
            status |= _gr_vec_add(pb, t, u, 2 * pow, ctx);

            left -= 2 * pow;
            pa = GR_ENTRY(pa, 2 * pow + 2, sz);
            pb = GR_ENTRY(pb, 2 * pow, sz);
        }

        if (left > pow)
        {
            status |= _gr_poly_mul(t, pa, pow + 1, GR_ENTRY(pb, pow, sz), left - pow, ctx);
            status |= _gr_poly_mul(u, pb, pow, GR_ENTRY(pa, pow + 1, sz), left - pow + 1, ctx);
            status |= _gr_vec_add(pb, t, u, left, ctx);
        }
    }

    GR_TMP_CLEAR_VEC(t, 2 * len, ctx);

    return status;
}

int
_gr_poly_interpolate_fast(gr_ptr poly,
    gr_srcptr xs, gr_srcptr ys, slong len, gr_ctx_t ctx)
{
    gr_ptr * tree;
    gr_ptr w;
    int status = GR_SUCCESS;

    tree = _gr_poly_tree_alloc(len, ctx);

    GR_TMP_INIT_VEC(w, len, ctx);

    status |= _gr_poly_tree_build(tree, xs, len, ctx);
    status |= _gr_poly_interpolation_weights(w, tree, len, ctx);
    status |= _gr_poly_interpolate_fast_precomp(poly, ys, tree, w, len, ctx);

    GR_TMP_CLEAR_VEC(w, len, ctx);
    _gr_poly_tree_free(tree, len, ctx);

    return status;
}

int
gr_poly_interpolate_fast(gr_poly_t poly, const gr_vec_t xs, const gr_vec_t ys, gr_ctx_t ctx)
{
    int status;
    slong n = xs->length;

    if (n != ys->length)
        return GR_DOMAIN;

    gr_poly_fit_length(poly, n, ctx);
    status = _gr_poly_interpolate_fast(poly->coeffs, xs->entries, ys->entries, n, ctx);
    _gr_poly_set_length(poly, n, ctx);
    _gr_poly_normalise(poly, ctx);
    return status;
}

