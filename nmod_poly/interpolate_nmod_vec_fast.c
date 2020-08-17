/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void
_nmod_poly_interpolation_weights(mp_ptr w, const mp_ptr * tree, slong len, nmod_t mod)
{
    mp_ptr tmp;
    slong i, n, height;

    if (len == 0)
        return;

    if (len == 1)
    {
        w[0] = 1;
        return;
    }

    tmp = _nmod_vec_init(len + 1);
    height = FLINT_CLOG2(len);
    n = WORD(1) << (height - 1);

    _nmod_poly_mul(tmp, tree[height-1], n + 1,
                        tree[height-1] + (n + 1), (len - n + 1), mod);

    _nmod_poly_derivative(tmp, tmp, len + 1, mod);
    _nmod_poly_evaluate_nmod_vec_fast_precomp(w, tmp, len, tree, len, mod);

    for (i = 0; i < len; i++)
        w[i] = n_invmod(w[i], mod.n);

    _nmod_vec_clear(tmp);
}

void
_nmod_poly_interpolate_nmod_vec_fast_precomp(mp_ptr poly, mp_srcptr ys,
    const mp_ptr * tree, mp_srcptr weights, slong len, nmod_t mod)
{
    mp_ptr t, u, pa, pb;
    slong i, pow, left;

    if (len == 0)
        return;

    t = _nmod_vec_init(len);
    u = _nmod_vec_init(len);

    for (i = 0; i < len; i++)
        poly[i] = nmod_mul(weights[i], ys[i], mod);

    for (i = 0; i < FLINT_CLOG2(len); i++)
    {
        pow = (WORD(1) << i);
        pa = tree[i];
        pb = poly;
        left = len;

        while (left >= 2 * pow)
        {
            _nmod_poly_mul(t, pa, pow + 1, pb + pow, pow, mod);
            _nmod_poly_mul(u, pa + pow + 1, pow + 1, pb, pow, mod);
            _nmod_vec_add(pb, t, u, 2 * pow, mod);

            left -= 2 * pow;
            pa += 2 * pow + 2;
            pb += 2 * pow;
        }

        if (left > pow)
        {
            _nmod_poly_mul(t, pa, pow + 1, pb + pow, left - pow, mod);
            _nmod_poly_mul(u, pb, pow, pa + pow + 1, left - pow + 1, mod);
            _nmod_vec_add(pb, t, u, left, mod);
        }
    }

    _nmod_vec_clear(t);
    _nmod_vec_clear(u);
}


void
_nmod_poly_interpolate_nmod_vec_fast(mp_ptr poly,
                            mp_srcptr xs, mp_srcptr ys, slong len, nmod_t mod)
{
    mp_ptr * tree;
    mp_ptr w;

    tree = _nmod_poly_tree_alloc(len);
    _nmod_poly_tree_build(tree, xs, len, mod);

    w = _nmod_vec_init(len);
    _nmod_poly_interpolation_weights(w, tree, len, mod);

    _nmod_poly_interpolate_nmod_vec_fast_precomp(poly, ys, tree, w, len, mod);

    _nmod_vec_clear(w);
    _nmod_poly_tree_free(tree, len);
}

void
nmod_poly_interpolate_nmod_vec_fast(nmod_poly_t poly,
                                    mp_srcptr xs, mp_srcptr ys, slong n)
{
    if (n == 0)
    {
        nmod_poly_zero(poly);
    }
    else
    {
        nmod_poly_fit_length(poly, n);
        poly->length = n;
        _nmod_poly_interpolate_nmod_vec_fast(poly->coeffs,
            xs, ys, n, poly->mod);
        _nmod_poly_normalise(poly);
    }
}
