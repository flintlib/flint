/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_interpolation_weights(acb_ptr w,
    acb_ptr * tree, slong len, slong prec)
{
    acb_ptr tmp;
    slong i, n, height;

    if (len == 0)
        return;

    if (len == 1)
    {
        acb_one(w);
        return;
    }

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

void
_acb_poly_interpolate_fast_precomp(acb_ptr poly,
    acb_srcptr ys, acb_ptr * tree, acb_srcptr weights,
    slong len, slong prec)
{
    acb_ptr t, u, pa, pb;
    slong i, pow, left;

    if (len == 0)
        return;

    t = _acb_vec_init(len);
    u = _acb_vec_init(len);

    for (i = 0; i < len; i++)
        acb_mul(poly + i, weights + i, ys + i, prec);

    for (i = 0; i < FLINT_CLOG2(len); i++)
    {
        pow = (WORD(1) << i);
        pa = tree[i];
        pb = poly;
        left = len;

        while (left >= 2 * pow)
        {
            _acb_poly_mul(t, pa, pow + 1, pb + pow, pow, prec);
            _acb_poly_mul(u, pa + pow + 1, pow + 1, pb, pow, prec);
            _acb_vec_add(pb, t, u, 2 * pow, prec);

            left -= 2 * pow;
            pa += 2 * pow + 2;
            pb += 2 * pow;
        }

        if (left > pow)
        {
            _acb_poly_mul(t, pa, pow + 1, pb + pow, left - pow, prec);
            _acb_poly_mul(u, pb, pow, pa + pow + 1, left - pow + 1, prec);
            _acb_vec_add(pb, t, u, left, prec);
        }
    }

    _acb_vec_clear(t, len);
    _acb_vec_clear(u, len);
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
