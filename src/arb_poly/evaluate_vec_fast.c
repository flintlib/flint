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
_arb_poly_evaluate_vec_fast_precomp(arb_ptr vs, arb_srcptr poly,
    slong plen, arb_ptr * tree, slong len, slong prec)
{
    gr_ctx_t ctx;

    /* _arb_poly_evaluate may differ from the generic Horner scheme */
    if (len == 1)
    {
        arb_t tmp;
        arb_init(tmp);
        arb_neg(tmp, tree[0] + 0);
        _arb_poly_evaluate(vs + 0, poly, plen, tmp, prec);
        arb_clear(tmp);
        return;
    }

    /* The tree polynomials are monic, so the division can not fail. */
    gr_ctx_init_real_arb(ctx, prec);
    GR_MUST_SUCCEED(_gr_poly_evaluate_vec_fast_precomp(vs, poly, plen,
        (const gr_ptr *) tree, len, ctx));
}

void _arb_poly_evaluate_vec_fast(arb_ptr ys, arb_srcptr poly, slong plen,
    arb_srcptr xs, slong n, slong prec)
{
    arb_ptr * tree;

    tree = _arb_poly_tree_alloc(n);
    _arb_poly_tree_build(tree, xs, n, prec);
    _arb_poly_evaluate_vec_fast_precomp(ys, poly, plen, tree, n, prec);
    _arb_poly_tree_free(tree, n);
}

void
arb_poly_evaluate_vec_fast(arb_ptr ys,
        const arb_poly_t poly, arb_srcptr xs, slong n, slong prec)
{
    _arb_poly_evaluate_vec_fast(ys, poly->coeffs,
                                        poly->length, xs, n, prec);
}
