/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "gr_poly.h"

void
_acb_poly_interpolation_weights(acb_ptr w,
    acb_ptr * tree, slong len, slong prec)
{
    gr_ctx_t ctx;

    gr_ctx_init_complex_acb(ctx, prec);

    /* Weights that are not invertible give GR_UNABLE, in which case
       gr_inv has already written the same (non-finite) value as acb_inv.
       Exact zeros give GR_DOMAIN and are not written; this can only
       happen for degenerate input, where we simply redo the computation
       with acb_inv to get the same output as acb_inv. */
    if (_gr_poly_interpolation_weights(w, (const gr_ptr *) tree, len, ctx) & GR_DOMAIN)
    {
        acb_ptr tmp;
        slong i, n, height;

        tmp = _acb_vec_init(len + 1);
        height = FLINT_CLOG2(len);
        n = WORD(1) << (height - 1);

        _acb_poly_mul_monic(tmp, tree[height-1], n + 1,
                            tree[height-1] + (n + 1), (len - n + 1), prec);
        _acb_poly_derivative(tmp, tmp, len + 1, prec);
        _acb_poly_evaluate_vec_fast_precomp(w, tmp, len, tree, len, prec);

        for (i = 0; i < len; i++)
            acb_inv(w + i, w + i, prec);

        _acb_vec_clear(tmp, len + 1);
    }
}

void
_acb_poly_interpolate_fast_precomp(acb_ptr poly,
    acb_srcptr ys, acb_ptr * tree, acb_srcptr weights,
    slong len, slong prec)
{
    gr_ctx_t ctx;
    gr_ctx_init_complex_acb(ctx, prec);
    GR_MUST_SUCCEED(_gr_poly_interpolate_fast_precomp(poly, ys,
        (const gr_ptr *) tree, weights, len, ctx));
}

void
_acb_poly_interpolate_fast(acb_ptr poly,
    acb_srcptr xs, acb_srcptr ys, slong len, slong prec)
{
    acb_ptr * tree;
    acb_ptr w;

    tree = _acb_poly_tree_alloc(len);
    _acb_poly_tree_build(tree, xs, len, prec);

    w = _acb_vec_init(len);
    _acb_poly_interpolation_weights(w, tree, len, prec);

    _acb_poly_interpolate_fast_precomp(poly, ys, tree, w, len, prec);

    _acb_vec_clear(w, len);
    _acb_poly_tree_free(tree, len);
}

void
acb_poly_interpolate_fast(acb_poly_t poly,
        acb_srcptr xs, acb_srcptr ys, slong n, slong prec)
{
    if (n == 0)
    {
        acb_poly_zero(poly);
    }
    else
    {
        acb_poly_fit_length(poly, n);
        _acb_poly_set_length(poly, n);
        _acb_poly_interpolate_fast(poly->coeffs, xs, ys, n, prec);
        _acb_poly_normalise(poly);
    }
}
