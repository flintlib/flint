/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2012, 2016, 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

/* todo: implement taylor shift */
#if 0
/* compose by poly2 = a*x^n + c, no aliasing; n >= 1 */
void
_ca_poly_compose_axnc(ca_ptr res, ca_srcptr poly1, slong len1,
    const ca_t c, const ca_t a, slong n, ca_ctx_t ctx)
{
    slong i;

    _ca_vec_set(res, poly1, len1, ctx);

    /* shift by c (c = 0 case will be fast) */
    _ca_poly_taylor_shift(res, c, len1, ctx);

    /* multiply by powers of a */
    if (!ca_is_one(a, ctx))
    {
        if (ca_equal_si(a, -1, ctx))
        {
            for (i = 1; i < len1; i += 2)
                ca_neg(res + i, res + i, ctx);
        }
        else if (len1 == 2)
        {
            ca_mul(res + 1, res + 1, a, ctx);
        }
        else
        {
            ca_t t;
            ca_init(t, ctx);
            ca_set(t, a, ctx);

            for (i = 1; i < len1; i++)
            {
                ca_mul(res + i, res + i, t, ctx);
                if (i + 1 < len1)
                    ca_mul(t, t, a, ctx);
            }

            ca_clear(t, ctx);
        }
    }

    /* stretch */
    for (i = len1 - 1; i >= 1 && n > 1; i--)
    {
        ca_swap(res + i * n, res + i, ctx);
        _ca_vec_zero(res + (i - 1) * n + 1, n - 1, ctx);
    }
}
#endif

void
_ca_poly_compose(ca_ptr res,
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
#if 0
    else if (_ca_vec_is_zero(poly2 + 1, len2 - 2))
    {
        _ca_poly_compose_axnc(res, poly1, len1, poly2, poly2 + len2 - 1, len2 - 1, ctx);
    }
#endif
    else if (len1 <= 7)
    {
        _ca_poly_compose_horner(res, poly1, len1, poly2, len2, ctx);
    }
    else
    {
        _ca_poly_compose_divconquer(res, poly1, len1, poly2, len2, ctx);
    }
}

void ca_poly_compose(ca_poly_t res,
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
            _ca_poly_compose(res->coeffs, poly1->coeffs, len1,
                                                   poly2->coeffs, len2, ctx);
        }
        else
        {
            ca_poly_t t;
            ca_poly_init2(t, lenr, ctx);
            _ca_poly_compose(t->coeffs, poly1->coeffs, len1,
                                                 poly2->coeffs, len2, ctx);
            ca_poly_swap(res, t, ctx);
            ca_poly_clear(t, ctx);
        }

        _ca_poly_set_length(res, lenr, ctx);
        _ca_poly_normalise(res, ctx);
    }
}
