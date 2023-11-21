/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

/* This gives some speedup for small lengths. */
static inline void
_arb_poly_rem_2(arb_ptr r, arb_srcptr a, slong al,
    arb_srcptr b, slong bl, slong prec)
{
    if (al == 2)
    {
        arb_mul(r + 0, a + 1, b + 0, prec);
        arb_sub(r + 0, a + 0, r + 0, prec);
    }
    else
    {
        _arb_poly_rem(r, a, al, b, bl, prec);
    }
}

void
_arb_poly_evaluate_vec_fast_precomp(arb_ptr vs, arb_srcptr poly,
    slong plen, arb_ptr * tree, slong len, slong prec)
{
    slong height, i, j, pow, left;
    slong tree_height;
    slong tlen;
    arb_ptr t, u, swap, pa, pb, pc;

    /* avoid worrying about some degenerate cases */
    if (len < 2 || plen < 2)
    {
        if (len == 1)
        {
            arb_t tmp;
            arb_init(tmp);
            arb_neg(tmp, tree[0] + 0);
            _arb_poly_evaluate(vs + 0, poly, plen, tmp, prec);
            arb_clear(tmp);
        }
        else if (len != 0 && plen == 0)
        {
            _arb_vec_zero(vs, len);
        }
        else if (len != 0 && plen == 1)
        {
            for (i = 0; i < len; i++)
                arb_set(vs + i, poly + 0);
        }
        return;
    }

    t = _arb_vec_init(len);
    u = _arb_vec_init(len);

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
        _arb_poly_rem(t + i, poly, plen, tree[height] + j, tlen + 1, prec);
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
            _arb_poly_rem_2(pc, pb, 2 * pow, pa, pow + 1, prec);
            _arb_poly_rem_2(pc + pow, pb, 2 * pow, pa + pow + 1, pow + 1, prec);

            pa += 2 * pow + 2;
            pb += 2 * pow;
            pc += 2 * pow;
            left -= 2 * pow;
        }

        if (left > pow)
        {
            _arb_poly_rem(pc, pb, left, pa, pow + 1, prec);
            _arb_poly_rem(pc + pow, pb, left, pa + pow + 1, left - pow + 1, prec);
        }
        else if (left > 0)
            _arb_vec_set(pc, pb, left);

        swap = t;
        t = u;
        u = swap;
    }

    _arb_vec_set(vs, t, len);
    _arb_vec_clear(t, len);
    _arb_vec_clear(u, len);
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
