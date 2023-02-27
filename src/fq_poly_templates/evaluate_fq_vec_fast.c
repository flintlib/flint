/*
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2012 William Hart
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
_TEMPLATE4(T, poly_evaluate, T, vec_fast_precomp)
    (TEMPLATE(T, struct) * vs,
     const TEMPLATE(T, struct) * poly, slong plen,
     TEMPLATE(T, poly_struct) * const * tree, slong len,
     const TEMPLATE(T, ctx_t) ctx)
{
    slong height, i, j, pow, left;
    slong tree_height;
    TEMPLATE(T, t) temp, inv;
    TEMPLATE(T, struct) * t, * u, * pb, * pc, * swap;
    TEMPLATE(T, poly_struct) * pa;

    TEMPLATE(T, init)(temp, ctx);
    TEMPLATE(T, init)(inv, ctx);

    /* avoid worrying about some degenerate cases */
    if (len < 2 || plen < 2)
    {
        if (len == 1)
        {
            TEMPLATE(T, neg)(temp, tree[0]->coeffs, ctx);
            _TEMPLATE3(T, poly_evaluate, T)(vs, poly, plen, temp, ctx);
        } else if (len != 0 && plen == 0)
            _TEMPLATE(T, vec_zero)(vs, len, ctx);
        else if (len != 0 && plen == 1)
            for (i = 0; i < len; i++)
                TEMPLATE(T, set)(vs + i, poly, ctx);
        
        TEMPLATE(T, clear)(temp, ctx);
        TEMPLATE(T, clear)(inv, ctx);
        return;
    }

    t = _TEMPLATE(T, vec_init)(2*len, ctx);
    u = _TEMPLATE(T, vec_init)(2*len, ctx);

    left = len;

    /* Initial reduction. We allow the polynomial to be larger
       or smaller than the number of points. */
    height = FLINT_BIT_COUNT(plen - 1) - 1;
    tree_height = FLINT_CLOG2(len);
    while (height >= tree_height)
        height--;
    pow = WORD(1) << height;

    for (i = j = 0; i < len; i += pow, j++)
    {
        pa = tree[height] + j;
        TEMPLATE(T, inv)(inv, pa->coeffs + pa->length - 1, ctx);
        _TEMPLATE(T, poly_rem)(t + i, poly, plen, pa->coeffs, pa->length, inv, ctx);
    }

    for (i = height - 1; i >= 0; i--)
    {
        pow = WORD(1) << i;
        left = len;
        pa = tree[i];
        pb = t;
        pc = u;

        left = len;
        while (left >= 2 * pow)
        {
            TEMPLATE(T, inv)(inv, pa->coeffs + pa->length - 1, ctx);
            _TEMPLATE(T, poly_rem)(pc, pb, 2 * pow, pa->coeffs, pa->length, inv, ctx);
            
            pa++;
            TEMPLATE(T, inv)(inv, pa->coeffs + pa->length - 1, ctx);
            _TEMPLATE(T, poly_rem)(pc + pow, pb, 2 * pow, pa->coeffs, pa->length, inv, ctx);
            
            pa++;
            pb += 2 * pow;
            pc += 2 * pow;
            left -= 2 * pow;
        }
        
        if (left > pow)
        {
            TEMPLATE(T, inv)(inv, pa->coeffs + pa->length - 1, ctx);
            _TEMPLATE(T, poly_rem)(pc, pb, left, pa->coeffs, pa->length, inv, ctx);
            
            pa ++;
            TEMPLATE(T, inv)(inv, pa->coeffs + pa->length - 1, ctx);
            _TEMPLATE(T, poly_rem)(pc + pow, pb, left, pa->coeffs, pa->length, inv, ctx);
        }
        else if (left > 0)
            _TEMPLATE(T, vec_set)(pc, pb, left, ctx);

        swap = t;
        t = u;
        u = swap;
    }

    TEMPLATE(T, clear)(temp, ctx);
    TEMPLATE(T, clear)(inv, ctx);

    _TEMPLATE(T, vec_set)(vs, t, len, ctx);

    _TEMPLATE(T, vec_clear)(t, 2*len, ctx);
    _TEMPLATE(T, vec_clear)(u, 2*len, ctx);
}

void
_TEMPLATE4(T, poly_evaluate, T, vec_fast)(TEMPLATE(T, struct) * ys,
                                          const TEMPLATE(T, struct) * poly, slong plen,
                                          const TEMPLATE(T, struct) * xs, slong n,
                                          const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_struct) ** tree;

    tree = _TEMPLATE(T, poly_tree_alloc)(n, ctx);
    _TEMPLATE(T, poly_tree_build)(tree, xs, n, ctx);
    _TEMPLATE4(T, poly_evaluate, T, vec_fast_precomp)(ys, poly, plen, tree, n, ctx);
    _TEMPLATE(T, poly_tree_free)(tree, n, ctx);
}

void
TEMPLATE4(T, poly_evaluate, T, vec_fast)(TEMPLATE(T, struct) * ys,
                                         const TEMPLATE(T, poly_t) poly,
                                         const TEMPLATE(T, struct) *xs, slong n,
                                         const TEMPLATE(T, ctx_t) ctx)
{
    _TEMPLATE4(T, poly_evaluate, T, vec_fast)(ys, poly->coeffs, poly->length,
                                              xs, n, ctx);
}


#endif
