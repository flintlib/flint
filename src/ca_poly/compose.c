/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2012, 2016, 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"
#include "gr_poly.h"

void
_ca_poly_compose(ca_ptr res,
    ca_srcptr poly1, slong len1,
    ca_srcptr poly2, slong len2, ca_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_ca_from_ref(gr_ctx, GR_CTX_CC_CA, ctx);
    GR_MUST_SUCCEED(_gr_poly_compose(res, poly1, len1, poly2, len2, gr_ctx));
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
