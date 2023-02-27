/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2012, 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

/* compose by poly2 = a*x^n + c, no aliasing; n >= 1 */
void
_arb_poly_compose_axnc(arb_ptr res, arb_srcptr poly1, slong len1,
    const arb_t c, const arb_t a, slong n, slong prec)
{
    slong i;

    _arb_vec_set_round(res, poly1, len1, prec);
    /* shift by c (c = 0 case will be fast) */
    _arb_poly_taylor_shift(res, c, len1, prec);

    /* multiply by powers of a */
    if (!arb_is_one(a))
    {
        if (arb_equal_si(a, -1))
        {
            for (i = 1; i < len1; i += 2)
                arb_neg(res + i, res + i);
        }
        else if (len1 == 2)
        {
            arb_mul(res + 1, res + 1, a, prec);
        }
        else
        {
            arb_t t;
            arb_init(t);
            arb_set(t, a);

            for (i = 1; i < len1; i++)
            {
                arb_mul(res + i, res + i, t, prec);
                if (i + 1 < len1)
                    arb_mul(t, t, a, prec);
            }

            arb_clear(t);
        }
    }

    /* stretch */
    for (i = len1 - 1; i >= 1 && n > 1; i--)
    {
        arb_swap(res + i * n, res + i);
        _arb_vec_zero(res + (i - 1) * n + 1, n - 1);
    }
}

void
_arb_poly_compose(arb_ptr res,
    arb_srcptr poly1, slong len1,
    arb_srcptr poly2, slong len2, slong prec)
{
    if (len1 == 1)
    {
        arb_set_round(res, poly1, prec);
    }
    else if (len2 == 1)
    {
        _arb_poly_evaluate(res, poly1, len1, poly2, prec);
    }
    else if (_arb_vec_is_zero(poly2 + 1, len2 - 2))
    {
        _arb_poly_compose_axnc(res, poly1, len1, poly2, poly2 + len2 - 1, len2 - 1, prec);
    }
    else if (len1 <= 7)
    {
        _arb_poly_compose_horner(res, poly1, len1, poly2, len2, prec);
    }
    else
    {
        _arb_poly_compose_divconquer(res, poly1, len1, poly2, len2, prec);
    }
}

void arb_poly_compose(arb_poly_t res,
              const arb_poly_t poly1, const arb_poly_t poly2, slong prec)
{
    const slong len1 = poly1->length;
    const slong len2 = poly2->length;
    
    if (len1 == 0)
    {
        arb_poly_zero(res);
    }
    else if (len1 == 1 || len2 == 0)
    {
        arb_poly_set_arb(res, poly1->coeffs);
    }
    else
    {
        const slong lenr = (len1 - 1) * (len2 - 1) + 1;
        
        if (res != poly1 && res != poly2)
        {
            arb_poly_fit_length(res, lenr);
            _arb_poly_compose(res->coeffs, poly1->coeffs, len1,
                                                   poly2->coeffs, len2, prec);
        }
        else
        {
            arb_poly_t t;
            arb_poly_init2(t, lenr);
            _arb_poly_compose(t->coeffs, poly1->coeffs, len1,
                                                 poly2->coeffs, len2, prec);
            arb_poly_swap(res, t);
            arb_poly_clear(t);
        }

        _arb_poly_set_length(res, lenr);
        _arb_poly_normalise(res);
    }
}

