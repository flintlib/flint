/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

void _ca_poly_mul(ca_ptr C,
    ca_srcptr A, slong lenA,
    ca_srcptr B, slong lenB, ca_ctx_t ctx)
{
    _ca_poly_mullow(C, A, lenA, B, lenB, lenA + lenB - 1, ctx);
}

void
ca_poly_mul(ca_poly_t res, const ca_poly_t poly1,
              const ca_poly_t poly2, ca_ctx_t ctx)
{
    slong len_out;

    if ((poly1->length == 0) || (poly2->length == 0))
    {
        ca_poly_zero(res, ctx);
        return;
    }

    len_out = poly1->length + poly2->length - 1;

    if (res == poly1 || res == poly2)
    {
        ca_poly_t temp;
        ca_poly_init2(temp, len_out, ctx);
        _ca_poly_mul(temp->coeffs, poly1->coeffs, poly1->length,
                                 poly2->coeffs, poly2->length, ctx);
        ca_poly_swap(res, temp, ctx);
        ca_poly_clear(temp, ctx);
    }
    else
    {
        ca_poly_fit_length(res, len_out, ctx);
        _ca_poly_mul(res->coeffs, poly1->coeffs, poly1->length,
                                 poly2->coeffs, poly2->length, ctx);
    }

    _ca_poly_set_length(res, len_out, ctx);
}
