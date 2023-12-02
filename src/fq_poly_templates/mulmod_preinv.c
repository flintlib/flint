/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Martin Lee
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
_TEMPLATE(T, poly_mulmod_preinv) (
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * poly1, slong len1,
    const TEMPLATE(T, struct) * poly2, slong len2,
    const TEMPLATE(T, struct) * f, slong lenf,
    const TEMPLATE(T, struct) * finv, slong lenfinv,
    const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, struct) * T, *Q;
    slong lenT, lenQ;

    lenT = len1 + len2 - 1;
    lenQ = lenT - lenf + 1;

    if (len1 + len2 > lenf) /* reduction necessary */
    {
        T = _TEMPLATE(T, vec_init) (lenT + lenQ, ctx);
        Q = T + lenT;

        if (len1 >= len2)
            _TEMPLATE(T, poly_mul) (T, poly1, len1, poly2, len2, ctx);
        else
            _TEMPLATE(T, poly_mul) (T, poly2, len2, poly1, len1, ctx);

        _TEMPLATE(T, poly_divrem_newton_n_preinv) (Q, res, T, lenT, f, lenf,
                                               finv, lenfinv, ctx);
        _TEMPLATE(T, vec_clear) (T, lenT + lenQ, ctx);
    } else /* just use mul */
    {
        if (len1 >= len2)
            _TEMPLATE(T, poly_mul) (res, poly1, len1, poly2, len2, ctx);
	else
	    _TEMPLATE(T, poly_mul) (res, poly2, len2, poly1, len1, ctx);

	if (lenT < lenf - 1)
            _TEMPLATE(T, vec_zero) (res + lenT, lenf - lenT - 1, ctx);
    }
}

void
TEMPLATE(T, poly_mulmod_preinv) (TEMPLATE(T, poly_t) res,
                                 const TEMPLATE(T, poly_t) poly1,
                                 const TEMPLATE(T, poly_t) poly2,
                                 const TEMPLATE(T, poly_t) f,
                                 const TEMPLATE(T, poly_t) finv,
                                 const TEMPLATE(T, ctx_t) ctx)
{
    slong len1, len2, lenf;
    TEMPLATE(T, struct) * fcoeffs, * coeffs1, * coeffs2;

    lenf = f->length;
    len1 = poly1->length;
    len2 = poly2->length;

    if (lenf == 0)
    {
        flint_throw(FLINT_ERROR, "(%s): Division by zero\n", __func__);
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

        if (poly1 == res)
        {
            coeffs1 = _TEMPLATE(T, vec_init) (len1, ctx);
            _TEMPLATE(T, vec_set) (coeffs1, poly1->coeffs, len1, ctx);
        }
        else
            coeffs1 = poly1->coeffs;

        if (poly2 == res)                                                                        {
            coeffs2 = _TEMPLATE(T, vec_init) (len2, ctx);
            _TEMPLATE(T, vec_set) (coeffs2, poly2->coeffs, len2, ctx);
        }
        else
            coeffs2 = poly2->coeffs;

         TEMPLATE(T, poly_fit_length) (res, lenf - 1, ctx);
        _TEMPLATE(T, poly_mulmod_preinv) (res->coeffs, coeffs1, len1,
                                          coeffs2, len2,
                                          fcoeffs, lenf, finv->coeffs,
                                          finv->length, ctx);
        if (f == res)
            _TEMPLATE(T, vec_clear) (fcoeffs, lenf, ctx);

        if (poly1 == res)
            _TEMPLATE(T, vec_clear) (coeffs1, len1, ctx);

        if (poly2 == res)
            _TEMPLATE(T, vec_clear) (coeffs2, len2, ctx);

        res->length = lenf - 1;
        _TEMPLATE(T, poly_normalise) (res, ctx);
    }
    else
    {
        TEMPLATE(T, poly_mul) (res, poly1, poly2, ctx);
    }
}


#endif
