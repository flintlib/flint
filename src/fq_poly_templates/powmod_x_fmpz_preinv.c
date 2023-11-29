/*
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
#include "long_extras.h"

void
_TEMPLATE(T, poly_powmod_x_fmpz_preinv) (
    TEMPLATE(T, struct) * res,
    const fmpz_t e,
    const TEMPLATE(T, struct) * f, slong lenf,
    const TEMPLATE(T, struct) * finv, slong lenfinv,
    const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, struct) * T, *Q;
    slong lenT, lenQ;
    slong i, window, l, c;

    lenT = 2 * lenf - 3;
    lenQ = lenT - lenf + 1;

    T = _TEMPLATE(T, vec_init) (lenT + lenQ, ctx);
    Q = T + lenT;

    TEMPLATE(T, one) (res, ctx);
    _TEMPLATE(T, vec_zero) (res + 1, lenf - 2, ctx);

    l = z_sizeinbase(lenf - 1, 2) - 2;
    window = 0;
    window = (1 << l);
    c = l;
    i = fmpz_sizeinbase(e, 2) - 2;

    if (i <= l)
    {
        window = 0;
        window = (1 << i);
        c = i;
        l = i;
    }

    if (c == 0)
    {
        _TEMPLATE(T, poly_shift_left) (T, res, lenf - 1, window, ctx);
        _TEMPLATE(T, poly_divrem_newton_n_preinv) (Q, res, T,
                                                   lenf - 1 + window, f, lenf,
                                                   finv, lenfinv, ctx);
        c = l + 1;
        window = 0;
    }

    for (; i >= 0; i--)
    {
        _TEMPLATE(T, poly_sqr) (T, res, lenf - 1, ctx);
        _TEMPLATE(T, poly_divrem_newton_n_preinv) (Q, res, T, 2 * lenf - 3, f,
                                                   lenf, finv, lenfinv, ctx);

        c--;
        if (fmpz_tstbit(e, i))
        {
            if (window == 0 && i <= l - 1)
                c = i;
            if (c >= 0)
                window = window | (1 << c);
        }
        else if (window == 0)
        {
            c = l + 1;
        }
        if (c == 0)
        {
            _TEMPLATE(T, poly_shift_left) (T, res, lenf - 1, window, ctx);

            _TEMPLATE(T, poly_divrem_newton_n_preinv) (Q, res, T,
                                                       lenf - 1 + window, f,
                                                       lenf, finv, lenfinv,
                                                       ctx);
            c = l + 1;
            window = 0;
        }
    }

    _TEMPLATE(T, vec_clear) (T, lenT + lenQ, ctx);
}


void
TEMPLATE(T, poly_powmod_x_fmpz_preinv) (TEMPLATE(T, poly_t) res,
                                        const fmpz_t e,
                                        const TEMPLATE(T, poly_t) f,
                                        const TEMPLATE(T, poly_t) finv,
                                        const TEMPLATE(T, ctx_t) ctx)
{
    slong lenf = f->length;
    slong trunc = lenf - 1;
    TEMPLATE(T, poly_t) tmp;

    if (lenf == 0)
    {
        flint_throw(FLINT_ERROR, "Exception: " TEMPLATE_STR(T) "_poly_powmod_x_preinv: divide by zero\n");
    }

    if (fmpz_sgn(e) < 0)
    {
        flint_throw(FLINT_ERROR, "Exception: " TEMPLATE_STR(T) "_poly_powmod_x_preinv: negative exp not implemented\n");
    }

    if (lenf == 1)
    {
        TEMPLATE(T, poly_zero) (res, ctx);
        return;
    }

    if (lenf == 2)
    {
        TEMPLATE(T, poly_t) r, poly;
        TEMPLATE(T, poly_init) (tmp, ctx);
        TEMPLATE(T, poly_init) (r, ctx);
        TEMPLATE(T, poly_init2) (poly, 2, ctx);
        TEMPLATE(T, poly_gen) (poly, ctx);
        TEMPLATE(T, poly_divrem) (tmp, r, poly, f, ctx);
        TEMPLATE(T, poly_powmod_fmpz_binexp_preinv) (res, r, e, f, finv, ctx);
        TEMPLATE(T, poly_clear) (tmp, ctx);
        TEMPLATE(T, poly_clear) (r, ctx);
        TEMPLATE(T, poly_clear) (poly, ctx);
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
                TEMPLATE(T, poly_t) r;
                TEMPLATE(T, poly_init2) (r, 2, ctx);
                TEMPLATE(T, poly_gen) (r, ctx);
                TEMPLATE(T, poly_init) (tmp, ctx);
                TEMPLATE(T, poly_divrem) (tmp, res, r, f, ctx);
                TEMPLATE(T, poly_clear) (tmp, ctx);
                TEMPLATE(T, poly_clear) (r, ctx);
            }
            else
            {
                TEMPLATE(T, poly_init2) (tmp, 3, ctx);
                TEMPLATE(T, poly_gen) (tmp, ctx);
                TEMPLATE(T, poly_mulmod) (res, tmp, tmp, f, ctx);
                TEMPLATE(T, poly_clear) (tmp, ctx);
            }
            return;
        }
    }

    if ((res == f) || (res == finv))
    {
        TEMPLATE(T, poly_init2) (tmp, trunc, ctx);
        _TEMPLATE(T, poly_powmod_x_fmpz_preinv) (tmp->coeffs, e, f->coeffs,
                                                 lenf, finv->coeffs,
                                                 finv->length, ctx);
        TEMPLATE(T, poly_swap) (res, tmp, ctx);
        TEMPLATE(T, poly_clear) (tmp, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length) (res, trunc, ctx);
        _TEMPLATE(T, poly_powmod_x_fmpz_preinv) (res->coeffs, e, f->coeffs,
                                                 lenf, finv->coeffs,
                                                 finv->length, ctx);
    }

    _TEMPLATE(T, poly_set_length) (res, trunc, ctx);
    _TEMPLATE(T, poly_normalise) (res, ctx);
}

#endif
