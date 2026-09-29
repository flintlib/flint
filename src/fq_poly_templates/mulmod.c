/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"
#include "gr_poly.h"

void
_TEMPLATE(T, poly_mulmod) (TEMPLATE(T, struct) * res,
                           const TEMPLATE(T, struct) * poly1, slong len1,
                           const TEMPLATE(T, struct) * poly2, slong len2,
                           const TEMPLATE(T, struct) * f, slong lenf,
                           const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    TEMPLATE3(_gr_ctx_init, T, from_ref)(gr_ctx, ctx);
    GR_MUST_SUCCEED(_gr_poly_mulmod(res, poly1, len1, poly2, len2,
                                    f, lenf, gr_ctx));
}

void
TEMPLATE(T, poly_mulmod) (TEMPLATE(T, poly_t) res,
                          const TEMPLATE(T, poly_t) poly1,
                          const TEMPLATE(T, poly_t) poly2,
                          const TEMPLATE(T, poly_t) f,
                          const TEMPLATE(T, ctx_t) ctx)
{
    slong len1, len2, lenf;
    TEMPLATE(T, struct) * fcoeffs;

    lenf = f->length;
    len1 = poly1->length;
    len2 = poly2->length;

    if (lenf == 0)
    {
        flint_throw(FLINT_DIVZERO, "(%s): Division by zero\n", __func__);
    }

    if (lenf == 1 || len1 == 0 || len2 == 0)
    {
        TEMPLATE(T, poly_zero) (res, ctx);
        return;
    }

    if (len1 + len2 - lenf > 0)
    {
        if (f == res)
        {
            fcoeffs = _TEMPLATE(T, vec_init) (lenf, ctx);
            _TEMPLATE(T, vec_set) (fcoeffs, f->coeffs, lenf, ctx);
        }
        else
            fcoeffs = f->coeffs;

        TEMPLATE(T, poly_fit_length) (res, len1 + len2 - 1, ctx);
        _TEMPLATE(T, poly_mulmod) (res->coeffs,
                                   poly1->coeffs, len1,
                                   poly2->coeffs, len2, fcoeffs, lenf, ctx);
        if (f == res)
            _TEMPLATE(T, vec_clear) (fcoeffs, lenf, ctx);

        _TEMPLATE(T, poly_set_length) (res, lenf - 1, ctx);
        _TEMPLATE(T, poly_normalise) (res, ctx);
    }
    else
    {
        TEMPLATE(T, poly_mul) (res, poly1, poly2, ctx);
    }
}


#endif
