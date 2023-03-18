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
#include "arb_poly.h"

void
_arb_poly_compose_divconquer(arb_ptr res, arb_srcptr poly1, slong len1,
                                          arb_srcptr poly2, slong len2, slong prec)
{
    if (len1 == 1)
    {
        arb_set(res, poly1);
    }
    else if (len2 == 1)
    {
        _arb_poly_evaluate(res, poly1, len1, poly2, prec);
    }
    else if (len1 == 2)
    {
        _arb_poly_compose_horner(res, poly1, len1, poly2, len2, prec);
    }
    else if (!_arb_vec_is_finite(poly1, len1) || !_arb_vec_is_finite(poly2, len2))
    {
        _arb_vec_indeterminate(res, (len1 - 1) * (len2 - 1) + 1);
    }
    else
    {
        gr_ctx_t ctx;
        gr_ctx_init_real_arb(ctx, prec);
        GR_MUST_SUCCEED(_gr_poly_compose_divconquer(res, poly1, len1, poly2, len2, ctx));
    }
}

void
arb_poly_compose_divconquer(arb_poly_t res,
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
            _arb_poly_compose_divconquer(res->coeffs, poly1->coeffs, len1,
                                                   poly2->coeffs, len2, prec);
        }
        else
        {
            arb_poly_t t;
            arb_poly_init2(t, lenr);
            _arb_poly_compose_divconquer(t->coeffs, poly1->coeffs, len1,
                                                 poly2->coeffs, len2, prec);
            arb_poly_swap(res, t);
            arb_poly_clear(t);
        }

        _arb_poly_set_length(res, lenr);
        _arb_poly_normalise(res);
    }
}
