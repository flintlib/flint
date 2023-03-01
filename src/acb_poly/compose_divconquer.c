/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"
#include "acb_poly.h"

void
_acb_poly_compose_divconquer(acb_ptr res, acb_srcptr poly1, slong len1,
                                          acb_srcptr poly2, slong len2, slong prec)
{
    if (len1 == 1)
    {
        acb_set(res, poly1);
    }
    else if (len2 == 1)
    {
        _acb_poly_evaluate(res, poly1, len1, poly2, prec);
    }
    else if (len1 == 2)
    {
        _acb_poly_compose_horner(res, poly1, len1, poly2, len2, prec);
    }
    else if (!_acb_vec_is_finite(poly1, len1) || !_acb_vec_is_finite(poly2, len2))
    {
        _acb_vec_indeterminate(res, (len1 - 1) * (len2 - 1) + 1);
    }
    else
    {
        gr_ctx_t ctx;
        gr_ctx_init_complex_acb(ctx, prec);
        GR_MUST_SUCCEED(_gr_poly_compose_divconquer(res, poly1, len1, poly2, len2, ctx));
    }
}

void
acb_poly_compose_divconquer(acb_poly_t res,
                             const acb_poly_t poly1, const acb_poly_t poly2, slong prec)
{
    const slong len1 = poly1->length;
    const slong len2 = poly2->length;
    
    if (len1 == 0)
    {
        acb_poly_zero(res);
    }
    else if (len1 == 1 || len2 == 0)
    {
        acb_poly_set_acb(res, poly1->coeffs);
    }
    else
    {
        const slong lenr = (len1 - 1) * (len2 - 1) + 1;
        
        if (res != poly1 && res != poly2)
        {
            acb_poly_fit_length(res, lenr);
            _acb_poly_compose_divconquer(res->coeffs, poly1->coeffs, len1,
                                                   poly2->coeffs, len2, prec);
        }
        else
        {
            acb_poly_t t;
            acb_poly_init2(t, lenr);
            _acb_poly_compose_divconquer(t->coeffs, poly1->coeffs, len1,
                                                 poly2->coeffs, len2, prec);
            acb_poly_swap(res, t);
            acb_poly_clear(t);
        }

        _acb_poly_set_length(res, lenr);
        _acb_poly_normalise(res);
    }
}
