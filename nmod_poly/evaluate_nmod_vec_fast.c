/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_poly.h"
#include "flint-impl.h"

/* This gives some speedup for small lengths. */
static __inline__ void _nmod_poly_rem_2(ulong_ptr r, ulong_srcptr a, slong al,
    ulong_srcptr b, slong bl, nmod_t mod)
{
    if (al == 2)
        r[0] = nmod_sub(a[0], nmod_mul(a[1], b[0], mod), mod);
    else
        _nmod_poly_rem(r, a, al, b, bl, mod);
}

void
_nmod_poly_evaluate_nmod_vec_fast_precomp(ulong_ptr vs, ulong_srcptr poly,
    slong plen, const ulong_ptr * tree, slong len, nmod_t mod)
{
    slong height, i, j, pow, left;
    slong tree_height;
    slong tlen;
    ulong_ptr t, u, swap, pa, pb, pc;

    /* avoid worrying about some degenerate cases */
    if (len < 2 || plen < 2)
    {
        if (len == 1)
            vs[0] = _nmod_poly_evaluate_nmod(poly, plen,
                nmod_neg(tree[0][0], mod), mod);
        else if (len != 0 && plen == 0)
            _NMOD_VEC_ZERO(vs, len);
        else if (len != 0 && plen == 1)
            for (i = 0; i < len; i++)
                vs[i] = poly[0];
        return;
    }

    t = _nmod_vec_init(len);
    u = _nmod_vec_init(len);

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
        _nmod_poly_rem(t + i, poly, plen, tree[height] + j, tlen + 1, mod);
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
            _nmod_poly_rem_2(pc, pb, 2 * pow, pa, pow + 1, mod);
            _nmod_poly_rem_2(pc + pow, pb, 2 * pow, pa + pow + 1, pow + 1, mod);

            pa += 2 * pow + 2;
            pb += 2 * pow;
            pc += 2 * pow;
            left -= 2 * pow;
        }

        if (left > pow)
        {
            _nmod_poly_rem(pc, pb, left, pa, pow + 1, mod);
            _nmod_poly_rem(pc + pow, pb, left, pa + pow + 1, left - pow + 1, mod);
        }
        else if (left > 0)
            _NMOD_VEC_SET(pc, pb, left);

        swap = t;
        t = u;
        u = swap;
    }

    _NMOD_VEC_SET(vs, t, len);
    _nmod_vec_clear(t);
    _nmod_vec_clear(u);
}

void _nmod_poly_evaluate_nmod_vec_fast(ulong_ptr ys, ulong_srcptr poly, slong plen,
    ulong_srcptr xs, slong n, nmod_t mod)
{
    ulong_ptr * tree;

    tree = _nmod_poly_tree_alloc(n);
    _nmod_poly_tree_build(tree, xs, n, mod);
    _nmod_poly_evaluate_nmod_vec_fast_precomp(ys, poly, plen, tree, n, mod);
    _nmod_poly_tree_free(tree, n);
}

void
nmod_poly_evaluate_nmod_vec_fast(ulong_ptr ys,
        const nmod_poly_t poly, ulong_srcptr xs, slong n)
{
    _nmod_poly_evaluate_nmod_vec_fast(ys, poly->coeffs,
                                        poly->length, xs, n, poly->mod);
}
