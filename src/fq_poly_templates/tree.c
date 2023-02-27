/*
    Copyright (C) 2011 Fredrik Johansson
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

TEMPLATE(T, poly_struct) **
_TEMPLATE(T, poly_tree_alloc)(slong len, const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_struct) ** tree = NULL;

    if (len)
    {
        slong i, j, height = FLINT_CLOG2(len);

        tree = flint_malloc(sizeof(TEMPLATE(T, poly_struct) *) * (height + 1));
        for (i = 0; i <= height; i++, len = (len + 1)/2)
        {
            tree[i] = flint_malloc(sizeof(TEMPLATE(T, poly_struct)) * len);
            for (j = 0; j < len; j++)
                TEMPLATE(T, poly_init)(tree[i] + j, ctx);
        }
            
    }

    return tree;
}

void
_TEMPLATE(T, poly_tree_free)(TEMPLATE(T, poly_struct) ** tree, slong len,
                             const TEMPLATE(T, ctx_t) ctx)
{
    if (len)
    {
        slong i, j, height = FLINT_CLOG2(len);

        for (i = 0; i <= height; i++, len = (len + 1)/2)
        {
           for (j = 0; j < len; j++)
               TEMPLATE(T, poly_clear)(tree[i] + j, ctx);
           flint_free(tree[i]);
        }

        flint_free(tree);
    }
}

void
_TEMPLATE(T, poly_tree_build)(TEMPLATE(T, poly_struct) ** tree,
                              const TEMPLATE(T, struct) * roots,
                              slong len,
                              const TEMPLATE(T, ctx_t) ctx)
{
    slong height, pow, left, i;
    TEMPLATE(T, poly_struct) * pa, * pb;

    if (len == 0)
        return;

    height = FLINT_CLOG2(len);

    /* zeroth level, (x-a) */
    for (i = 0; i < len; i++)
    {
        TEMPLATE(T, poly_gen)(tree[0] + i, ctx);
        TEMPLATE(T, neg)((tree[0] + i)->coeffs, roots + i, ctx);
    }

    for (i = 0; i < height - 1; i++)
    {
        left = len;
        pow = WORD(1) << i;
        pa = tree[i];
        pb = tree[i + 1];

        while (left >= 2 * pow)
        {
            TEMPLATE(T, poly_fit_length)(pb, pa->length + (pa + 1)->length - 1,
                                         ctx);
            _TEMPLATE(T, poly_mul)(pb->coeffs,
                                   pa->coeffs, pa->length,
                                   (pa + 1)->coeffs, (pa + 1)->length,
                                   ctx);
            _TEMPLATE(T, poly_set_length)(pb, pa->length + (pa + 1)->length - 1,
                                          ctx);
            left -= 2 * pow;
            pa += 2;
            pb += 1;
        }

        if (left > pow)
        {
            TEMPLATE(T, poly_fit_length)(pb, pa->length + (pa + 1)->length - 1,
                                         ctx);
            _TEMPLATE(T, poly_mul)(pb->coeffs,
                                   pa->coeffs, pa->length,
                                   (pa + 1)->coeffs, (pa + 1)->length,
                                   ctx);
            _TEMPLATE(T, poly_set_length)(pb, pa->length + (pa + 1)->length - 1,
                                          ctx);
        } else if (left > 0)
            TEMPLATE(T, poly_set)(pb, pa, ctx);
    }
}


#endif
