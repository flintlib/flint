/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void _acb_poly_mul(acb_ptr C,
    acb_srcptr A, slong lenA,
    acb_srcptr B, slong lenB, slong prec)
{
    _acb_poly_mullow(C, A, lenA, B, lenB, lenA + lenB - 1, prec);
}

void
acb_poly_mul(acb_poly_t res, const acb_poly_t poly1,
              const acb_poly_t poly2, slong prec)
{
    slong len_out;

    if ((poly1->length == 0) || (poly2->length == 0))
    {
        acb_poly_zero(res);
        return;
    }

    len_out = poly1->length + poly2->length - 1;

    if (res == poly1 || res == poly2)
    {
        acb_poly_t temp;
        acb_poly_init2(temp, len_out);
        _acb_poly_mul(temp->coeffs, poly1->coeffs, poly1->length,
                                 poly2->coeffs, poly2->length, prec);
        acb_poly_swap(res, temp);
        acb_poly_clear(temp);
    }
    else
    {
        acb_poly_fit_length(res, len_out);
        _acb_poly_mul(res->coeffs, poly1->coeffs, poly1->length,
                                 poly2->coeffs, poly2->length, prec);
    }

    _acb_poly_set_length(res, len_out);
    _acb_poly_normalise(res);
}

