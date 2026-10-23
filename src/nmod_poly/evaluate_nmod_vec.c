/*
    Copyright (C) 2011, 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "gr_poly.h"

void
_nmod_poly_evaluate_nmod_vec(nn_ptr ys, nn_srcptr coeffs, slong len,
    nn_srcptr xs, slong n, nmod_t mod)
{
    if (len < 32)
        _nmod_poly_evaluate_nmod_vec_iter(ys, coeffs, len, xs, n, mod);
    else
        _nmod_poly_evaluate_nmod_vec_fast(ys, coeffs, len, xs, n, mod);
}

void
nmod_poly_evaluate_nmod_vec(nn_ptr ys,
        const nmod_poly_t poly, nn_srcptr xs, slong n)
{
    _nmod_poly_evaluate_nmod_vec(ys, poly->coeffs,
                                        poly->length, xs, n, poly->mod);
}

void
_nmod_poly_evaluate_nmod_vec_fast_precomp(nn_ptr vs, nn_srcptr poly,
    slong plen, const nn_ptr * tree, slong len, nmod_t mod)
{
    gr_ctx_t ctx;
    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(_gr_poly_evaluate_vec_fast_precomp(vs, poly, plen, (const gr_ptr *) tree, len, ctx));
}

void _nmod_poly_evaluate_nmod_vec_fast(nn_ptr ys, nn_srcptr poly, slong plen,
    nn_srcptr xs, slong n, nmod_t mod)
{
    nn_ptr * tree;

    tree = _nmod_poly_tree_alloc(n);
    _nmod_poly_tree_build(tree, xs, n, mod);
    _nmod_poly_evaluate_nmod_vec_fast_precomp(ys, poly, plen, tree, n, mod);
    _nmod_poly_tree_free(tree, n);
}

void
nmod_poly_evaluate_nmod_vec_fast(nn_ptr ys,
        const nmod_poly_t poly, nn_srcptr xs, slong n)
{
    _nmod_poly_evaluate_nmod_vec_fast(ys, poly->coeffs,
                                        poly->length, xs, n, poly->mod);
}

void
_nmod_poly_evaluate_nmod_vec_iter(nn_ptr ys, nn_srcptr coeffs, slong len,
    nn_srcptr xs, slong n, nmod_t mod)
{
    slong i;
    for (i = 0; i < n; i++)
        ys[i] = _nmod_poly_evaluate_nmod(coeffs, len, xs[i], mod);
}

void
nmod_poly_evaluate_nmod_vec_iter(nn_ptr ys,
        const nmod_poly_t poly, nn_srcptr xs, slong n)
{
    _nmod_poly_evaluate_nmod_vec_iter(ys, poly->coeffs,
                                        poly->length, xs, n, poly->mod);
}
