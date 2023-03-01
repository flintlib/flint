/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

void
_ca_poly_compose_horner(ca_ptr res,
    ca_srcptr poly1, slong len1,
    ca_srcptr poly2, slong len2, ca_ctx_t ctx)
{
    if (len1 == 1)
    {
        ca_set(res, poly1, ctx);
    }
    else if (len2 == 1)
    {
        _ca_poly_evaluate(res, poly1, len1, poly2, ctx);
    }
    else if (len1 == 2)
    {
        _ca_vec_scalar_mul_ca(res, poly2, len2, poly1 + 1, ctx);
        ca_add(res, res, poly1, ctx);
    }
    else
    {
        const slong alloc = (len1 - 1) * (len2 - 1) + 1;
        slong i = len1 - 1, lenr = len2;
        ca_ptr t, t1, t2;
        t = _ca_vec_init(alloc, ctx);

        if (len1 % 2 == 0)
        {
            t1 = res;
            t2 = t;
        }
        else
        {
            t1 = t;
            t2 = res;
        }

        /* Perform the first two steps as one,
            "res = a(m) * poly2 + a(m-1)". */
        {
            _ca_vec_scalar_mul_ca(t1, poly2, len2, poly1 + i, ctx);
            i--;
            ca_add(t1 + 0, t1 + 0, poly1 + i, ctx);
        }
        while (i--)
        {
            _ca_poly_mul(t2, t1, lenr, poly2, len2, ctx);
            lenr += len2 - 1;
            {
                void *t_ = t1;
                t1 = t2;
                t2 = t_;
            }
            ca_add(t1 + 0, t1 + 0, poly1 + i, ctx);
        }
        _ca_vec_clear(t, alloc, ctx);
    }
}

void ca_poly_compose_horner(ca_poly_t res,
              const ca_poly_t poly1, const ca_poly_t poly2, ca_ctx_t ctx)
{
    const slong len1 = poly1->length;
    const slong len2 = poly2->length;
    
    if (len1 == 0)
    {
        ca_poly_zero(res, ctx);
    }
    else if (len1 == 1 || len2 == 0)
    {
        ca_poly_set_ca(res, poly1->coeffs, ctx);
    }
    else
    {
        const slong lenr = (len1 - 1) * (len2 - 1) + 1;
        
        if (res != poly1 && res != poly2)
        {
            ca_poly_fit_length(res, lenr, ctx);
            _ca_poly_compose_horner(res->coeffs, poly1->coeffs, len1,
                                                   poly2->coeffs, len2, ctx);
        }
        else
        {
            ca_poly_t t;
            ca_poly_init2(t, lenr, ctx);
            _ca_poly_compose_horner(t->coeffs, poly1->coeffs, len1,
                                                 poly2->coeffs, len2, ctx);
            ca_poly_swap(res, t, ctx);
            ca_poly_clear(t, ctx);
        }

        _ca_poly_set_length(res, lenr, ctx);
        _ca_poly_normalise(res, ctx);
    }
}
