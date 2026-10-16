/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "gr_poly.h"

void
_arb_poly_interpolation_weights(arb_ptr w,
    arb_ptr * tree, slong len, slong prec)
{
    gr_ctx_t ctx;

    gr_ctx_init_real_arb(ctx, prec);

    /* Weights that are not invertible give GR_UNABLE, in which case
       gr_inv has already written the same (non-finite) value as arb_inv.
       Exact zeros give GR_DOMAIN and are not written; this can only
       happen for degenerate input, where we simply redo the computation
       with arb_inv to get the same output as arb_inv. */
    if (_gr_poly_interpolation_weights(w, (const gr_ptr *) tree, len, ctx) & GR_DOMAIN)
    {
        arb_ptr tmp;
        slong i, n, height;

        tmp = _arb_vec_init(len + 1);
        height = FLINT_CLOG2(len);
        n = WORD(1) << (height - 1);

        _arb_poly_mul_monic(tmp, tree[height-1], n + 1,
                            tree[height-1] + (n + 1), (len - n + 1), prec);
        _arb_poly_derivative(tmp, tmp, len + 1, prec);
        _arb_poly_evaluate_vec_fast_precomp(w, tmp, len, tree, len, prec);

        for (i = 0; i < len; i++)
            arb_inv(w + i, w + i, prec);

        _arb_vec_clear(tmp, len + 1);
    }
}

void
_arb_poly_interpolate_fast_precomp(arb_ptr poly,
    arb_srcptr ys, arb_ptr * tree, arb_srcptr weights,
    slong len, slong prec)
{
    gr_ctx_t ctx;
    gr_ctx_init_real_arb(ctx, prec);
    GR_MUST_SUCCEED(_gr_poly_interpolate_fast_precomp(poly, ys,
        (const gr_ptr *) tree, weights, len, ctx));
}

void
_arb_poly_interpolate_fast(arb_ptr poly,
    arb_srcptr xs, arb_srcptr ys, slong len, slong prec)
{
    arb_ptr * tree;
    arb_ptr w;

    tree = _arb_poly_tree_alloc(len);
    _arb_poly_tree_build(tree, xs, len, prec);

    w = _arb_vec_init(len);
    _arb_poly_interpolation_weights(w, tree, len, prec);

    _arb_poly_interpolate_fast_precomp(poly, ys, tree, w, len, prec);

    _arb_vec_clear(w, len);
    _arb_poly_tree_free(tree, len);
}

void
arb_poly_interpolate_fast(arb_poly_t poly,
        arb_srcptr xs, arb_srcptr ys, slong n, slong prec)
{
    if (n == 0)
    {
        arb_poly_zero(poly);
    }
    else
    {
        arb_poly_fit_length(poly, n);
        _arb_poly_set_length(poly, n);
        _arb_poly_interpolate_fast(poly->coeffs, xs, ys, n, prec);
        _arb_poly_normalise(poly);
    }
}
