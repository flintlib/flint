/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

/* This gives some speedup for small lengths. */
static inline void
_acb_poly_rem_2(acb_ptr r, acb_srcptr a, slong al,
    acb_srcptr b, slong bl, slong prec)
{
    if (al == 2)
    {
        acb_mul(r + 0, a + 1, b + 0, prec);
        acb_sub(r + 0, a + 0, r + 0, prec);
    }
    else
    {
        _acb_poly_rem(r, a, al, b, bl, prec);
    }
}

void
_acb_poly_evaluate_vec_fast_precomp(acb_ptr vs, acb_srcptr poly,
    slong plen, acb_ptr * tree, slong len, slong prec)
{
    slong height, i, j, pow, left;
    slong tree_height;
    slong tlen;
    acb_ptr t, u, swap, pa, pb, pc;

    /* avoid worrying about some degenerate cases */
    if (len < 2 || plen < 2)
    {
        if (len == 1)
        {
            acb_t tmp;
            acb_init(tmp);
            acb_neg(tmp, tree[0] + 0);
            _acb_poly_evaluate(vs + 0, poly, plen, tmp, prec);
            acb_clear(tmp);
        }
        else if (len != 0 && plen == 0)
        {
            _acb_vec_zero(vs, len);
        }
        else if (len != 0 && plen == 1)
        {
            for (i = 0; i < len; i++)
                acb_set(vs + i, poly + 0);
        }
        return;
    }

    t = _acb_vec_init(len);
    u = _acb_vec_init(len);

    left = len;

    /* Initial reduction. We allow the polynomial to be larger
        or smaller than the number of points. */
    height = FLINT_BIT_COUNT(plen - 1) - 1;
    tree_height = FLINT_CLOG2(len);
    while (height >= tree_height)
        height--;
    pow = WORD(1) << height;

    for (i = j = 0; i < len; i += pow, j += (pow + 1))
    {
        tlen = ((i + pow) <= len) ? pow : len % pow;
        _acb_poly_rem(t + i, poly, plen, tree[height] + j, tlen + 1, prec);
    }

    for (i = height - 1; i >= 0; i--)
    {
        pow = WORD(1) << i;
        left = len;
        pa = tree[i];
        pb = t;
        pc = u;

        while (left >= 2 * pow)
        {
            _acb_poly_rem_2(pc, pb, 2 * pow, pa, pow + 1, prec);
            _acb_poly_rem_2(pc + pow, pb, 2 * pow, pa + pow + 1, pow + 1, prec);

            pa += 2 * pow + 2;
            pb += 2 * pow;
            pc += 2 * pow;
            left -= 2 * pow;
        }

        if (left > pow)
        {
            _acb_poly_rem(pc, pb, left, pa, pow + 1, prec);
            _acb_poly_rem(pc + pow, pb, left, pa + pow + 1, left - pow + 1, prec);
        }
        else if (left > 0)
            _acb_vec_set(pc, pb, left);

        swap = t;
        t = u;
        u = swap;
    }

    _acb_vec_set(vs, t, len);
    _acb_vec_clear(t, len);
    _acb_vec_clear(u, len);
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
