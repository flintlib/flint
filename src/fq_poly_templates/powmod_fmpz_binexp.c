/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include "flint.h"
void
_TEMPLATE(T, poly_powmod_fmpz_binexp) (
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * poly,
    const fmpz_t e,
    const TEMPLATE(T, struct) * f, slong lenf,
    const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, struct) * T, *Q;
    TEMPLATE(T, t) invf;
    slong lenT, lenQ;
    slong i;

    if (lenf == 2)
    {
        TEMPLATE(T, pow) (res, poly, e, ctx);
        return;
    }

    lenT = 2 * lenf - 3;
    lenQ = lenT - lenf + 1;

    T = _TEMPLATE(T, vec_init) (lenT + lenQ, ctx);
    Q = T + lenT;

    TEMPLATE(T, init) (invf, ctx);
    TEMPLATE(T, inv) (invf, f + lenf - 1, ctx);

    _TEMPLATE(T, vec_set) (res, poly, lenf - 1, ctx);

    for (i = fmpz_sizeinbase(e, 2) - 2; i >= 0; i--)
    {
        _TEMPLATE(T, poly_sqr) (T, res, lenf - 1, ctx);
        _TEMPLATE(T, poly_divrem) (Q, res, T, 2 * lenf - 3, f, lenf, invf,
                                   ctx);

        if (fmpz_tstbit(e, i))
        {
            _TEMPLATE(T, poly_mul) (T, res, lenf - 1, poly, lenf - 1, ctx);
            _TEMPLATE(T, poly_divrem) (Q, res, T, 2 * lenf - 3, f, lenf, invf,
                                       ctx);
        }
    }

    TEMPLATE(T, clear) (invf, ctx);
    _TEMPLATE(T, vec_clear) (T, lenT + lenQ, ctx);
}


void
TEMPLATE(T, poly_powmod_fmpz_binexp) (TEMPLATE(T, poly_t) res,
                                      const TEMPLATE(T, poly_t) poly,
                                      const fmpz_t e,
                                      const TEMPLATE(T, poly_t) f,
                                      const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, struct) * q;
    slong len = poly->length;
    slong lenf = f->length;
    slong trunc = lenf - 1;
    int qcopy = 0;

    if (lenf == 0)
    {
        flint_throw(FLINT_ERROR, "(%s): Division by zero\n", __func__);
    }

    if (fmpz_sgn(e) < 0)
    {
        flint_throw(FLINT_ERROR, "(%s): Negative exp not implemented\n", __func__);
    }

    if (len >= lenf)
    {
        TEMPLATE(T, poly_t) t, r;
        TEMPLATE(T, poly_init) (t, ctx);
        TEMPLATE(T, poly_init) (r, ctx);
        TEMPLATE(T, poly_divrem) (t, r, poly, f, ctx);
        TEMPLATE(T, poly_powmod_fmpz_binexp) (res, r, e, f, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
        TEMPLATE(T, poly_clear) (r, ctx);
        return;
    }

    if (fmpz_abs_fits_ui(e))
    {
        ulong exp = fmpz_get_ui(e);

        if (exp <= 2)
        {
            if (exp == UWORD(0))
            {
                TEMPLATE(T, poly_fit_length) (res, 1, ctx);
                TEMPLATE(T, one) (res->coeffs, ctx);
                _TEMPLATE(T, poly_set_length) (res, 1, ctx);
            }
            else if (exp == UWORD(1))
            {
                TEMPLATE(T, poly_set) (res, poly, ctx);
            }
            else
                TEMPLATE(T, poly_mulmod) (res, poly, poly, f, ctx);
            return;
        }
    }

    if (lenf == 1 || len == 0)
    {
        TEMPLATE(T, poly_zero) (res, ctx);
        return;
    }

    if (poly->length < trunc)
    {
        q = _TEMPLATE(T, vec_init) (trunc, ctx);
        _TEMPLATE(T, vec_set) (q, poly->coeffs, len, ctx);
        _TEMPLATE(T, vec_zero) (q + len, trunc - len, ctx);
        qcopy = 1;
    }
    else
    {
        q = poly->coeffs;
    }

    if ((res == poly && !qcopy) || (res == f))
    {
        TEMPLATE(T, poly_t) t;
        TEMPLATE(T, poly_init2) (t, 2 * lenf - 3, ctx);
        _TEMPLATE(T, poly_powmod_fmpz_binexp) (t->coeffs, q, e, f->coeffs,
                                               lenf, ctx);
        TEMPLATE(T, poly_swap) (res, t, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length) (res, 2 * lenf - 3, ctx);
        _TEMPLATE(T, poly_powmod_fmpz_binexp) (res->coeffs, q, e, f->coeffs,
                                               lenf, ctx);
    }

    if (qcopy)
        _TEMPLATE(T, vec_clear) (q, trunc, ctx);

    _TEMPLATE(T, poly_set_length) (res, trunc, ctx);
    _TEMPLATE(T, poly_normalise) (res, ctx);
}


#endif
