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
_acb_poly_evaluate_vec_fast_precomp(acb_ptr vs, acb_srcptr poly,
    slong plen, acb_ptr * tree, slong len, slong prec)
{
    gr_ctx_t ctx;

    /* _acb_poly_evaluate may differ from the generic Horner scheme */
    if (len == 1)
    {
        acb_t tmp;
        acb_init(tmp);
        acb_neg(tmp, tree[0] + 0);
        _acb_poly_evaluate(vs + 0, poly, plen, tmp, prec);
        acb_clear(tmp);
        return;
    }

    /* The tree polynomials are monic, so the division can not fail. */
    gr_ctx_init_complex_acb(ctx, prec);
    GR_MUST_SUCCEED(_gr_poly_evaluate_vec_fast_precomp(vs, poly, plen,
        (const gr_ptr *) tree, len, ctx));
}

void _acb_poly_evaluate_vec_fast(acb_ptr ys, acb_srcptr poly, slong plen,
    acb_srcptr xs, slong n, slong prec)
{
    acb_ptr * tree;

    tree = _acb_poly_tree_alloc(n);
    _acb_poly_tree_build(tree, xs, n, prec);
    _acb_poly_evaluate_vec_fast_precomp(ys, poly, plen, tree, n, prec);
    _acb_poly_tree_free(tree, n);
}

void
acb_poly_evaluate_vec_fast(acb_ptr ys,
        const acb_poly_t poly, acb_srcptr xs, slong n, slong prec)
{
    _acb_poly_evaluate_vec_fast(ys, poly->coeffs,
                                        poly->length, xs, n, prec);
}
