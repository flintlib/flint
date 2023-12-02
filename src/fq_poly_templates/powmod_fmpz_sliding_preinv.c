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

void
_TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * poly,
    const fmpz_t e, ulong k,
    const TEMPLATE(T, struct) * f, slong lenf,
    const TEMPLATE(T, struct) * finv, slong lenfinv,
    const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, struct) * T, *Q;
    TEMPLATE(T, poly_struct) * precomp;
    TEMPLATE(T, poly_t) poly_squared;
    ulong twokm1;
    slong lenT, lenQ;
    slong i, l, j;
    int index;

    if (lenf == 2)
    {
        TEMPLATE(T, pow) (res, poly, e, ctx);
        return;
    }

    lenT = 2 * lenf - 3;
    lenQ = lenT - lenf + 1;

    T = _TEMPLATE(T, vec_init) (lenT + lenQ, ctx);
    Q = T + lenT;

    /* Precomputation */
    twokm1 = n_pow(2, k - 1);
    precomp = flint_malloc(twokm1 * sizeof(TEMPLATE(T, poly_struct)));
    TEMPLATE(T, poly_init) (precomp, ctx);
    TEMPLATE(T, poly_fit_length) (precomp, lenf - 1, ctx);
    _TEMPLATE(T, vec_set) (precomp->coeffs, poly, lenf - 1, ctx);

    TEMPLATE(T, poly_init) (poly_squared, ctx);
    if (k > 1)
    {
        TEMPLATE(T, poly_fit_length) (poly_squared, lenf - 1, ctx);
        _TEMPLATE(T, poly_mul) (T, poly, lenf - 1, poly, lenf - 1, ctx);
        _TEMPLATE(T, poly_divrem_newton_n_preinv) (Q, poly_squared->coeffs, T,
                                                   2 * lenf - 3, f, lenf, finv,
                                                   lenfinv, ctx);
    }
    for (i = 1; i < twokm1; i++)
    {
        TEMPLATE(T, poly_init) (precomp + i, ctx);
        TEMPLATE(T, poly_fit_length) (precomp + i, lenf - 1, ctx);
        _TEMPLATE(T, poly_mul) (T, (precomp + i - 1)->coeffs, lenf - 1,
                                poly_squared->coeffs, lenf - 1, ctx);
        _TEMPLATE(T, poly_divrem_newton_n_preinv) (Q, (precomp + i)->coeffs, T,
                                                   2 * lenf - 3, f, lenf, finv,
                                                   lenfinv, ctx);
    }

    _TEMPLATE(T, vec_set) (res, poly, lenf - 1, ctx);

    i = fmpz_sizeinbase(e, 2) - 2;
    while (i >= 0)
    {
        if (fmpz_tstbit(e, i) == 0)
        {
            _TEMPLATE(T, poly_sqr) (T, res, lenf - 1, ctx);
            _TEMPLATE(T, poly_divrem_newton_n_preinv) (Q, res, T, 2 * lenf - 3,
                                                       f, lenf, finv, lenfinv,
                                                       ctx);
            i -= 1;
        }
        else
        {
            l = FLINT_MAX(i - k + 1, 0);
            while (fmpz_tstbit(e, l) == 0)
            {
                l += 1;
            }
            for (j = 0; j < i - l + 1; j++)
            {
                _TEMPLATE(T, poly_sqr) (T, res, lenf - 1, ctx);
                _TEMPLATE(T, poly_divrem_newton_n_preinv) (Q, res, T,
                                                           2 * lenf - 3, f,
                                                           lenf, finv, lenfinv,
                                                           ctx);
            }

            index = fmpz_tstbit(e, i);
            for (j = i - 1; j >= l; j--)
            {
                index = index << 1;
                index += fmpz_tstbit(e, j);
            }
            index = (index - 1) / 2;

            _TEMPLATE(T, poly_mul) (T, res, lenf - 1,
                                    (precomp + index)->coeffs, lenf - 1, ctx);
            _TEMPLATE(T, poly_divrem_newton_n_preinv) (Q, res, T, 2 * lenf - 3,
                                                       f, lenf, finv, lenfinv,
                                                       ctx);
            i = l - 1;
        }
    }

    for (j = 0; j < twokm1; j++)
    {
        TEMPLATE(T, poly_clear) (precomp + j, ctx);
    }
    flint_free(precomp);
    TEMPLATE(T, poly_clear) (poly_squared, ctx);
    _TEMPLATE(T, vec_clear) (T, lenT + lenQ, ctx);
}


void
TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (TEMPLATE(T, poly_t) res,
                                              const TEMPLATE(T, poly_t) poly,
                                              const fmpz_t e, ulong k,
                                              const TEMPLATE(T, poly_t) f,
                                              const TEMPLATE(T, poly_t) finv,
                                              const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, struct) * q;
    slong len = poly->length;
    slong lenf = f->length;
    slong trunc = lenf - 1;
    int qcopy = 0;
    flint_bitcnt_t bits;

    if (lenf == 0)
    {
        flint_throw(FLINT_ERROR, "Exception: " TEMPLATE_STR(T) "_poly_powmod_fmpz_sliding_preinv: divide by zero\n");
    }

    if (fmpz_sgn(e) < 0)
    {
        flint_throw(FLINT_ERROR, "Exception: " TEMPLATE_STR(T) "_poly_powmod_fmpz_sliding_preinv: negative exp not implemented\n");
    }

    if (len >= lenf)
    {
        TEMPLATE(T, poly_t) t, r;
        TEMPLATE(T, poly_init) (t, ctx);
        TEMPLATE(T, poly_init) (r, ctx);
        TEMPLATE(T, poly_divrem) (t, r, poly, f, ctx);
        TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (res, r, e, k, f, finv,
                                                      ctx);
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
                TEMPLATE(T, poly_mulmod_preinv) (res, poly, poly, f, finv,
                                                 ctx);
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

    /* Determine "optimum" sliding window size */
    if (k == 0)
    {
        bits = fmpz_bits(e);
        if (bits < 9)
            k = 1;
        else if (bits < 15)
            k = 2;
        else if (bits < 62)
            k = 3;
        else if (bits < 203)
            k = 4;
        else if (bits < 587)
            k = 5;
        else if (bits < 1560)
            k = 6;
        else
            k = 7;
    }

    if ((res == poly && !qcopy) || (res == f))
    {
        TEMPLATE(T, poly_t) t;
        TEMPLATE(T, poly_init2) (t, 2 * lenf - 3, ctx);
        _TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (t->coeffs, q, e, k,
                                                       f->coeffs, lenf,
                                                       finv->coeffs,
                                                       finv->length, ctx);
        TEMPLATE(T, poly_swap) (res, t, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length) (res, 2 * lenf - 3, ctx);
        _TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (res->coeffs, q, e, k,
                                                       f->coeffs, lenf,
                                                       finv->coeffs,
                                                       finv->length, ctx);
    }

    if (qcopy)
        _TEMPLATE(T, vec_clear) (q, trunc, ctx);

    _TEMPLATE(T, poly_set_length) (res, trunc, ctx);
    _TEMPLATE(T, poly_normalise) (res, ctx);
}


#endif
