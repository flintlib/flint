/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_mulmod(fmpz * res, const fmpz * poly1, slong len1,
                           const fmpz * poly2, slong len2, const fmpz * f,
                           slong lenf, const fmpz_mod_ctx_t ctx)
{
    fmpz * T, * Q;
    fmpz_t invf;
    slong lenT, lenQ;

    lenT = len1 + len2 - 1;
    lenQ = lenT - lenf + 1;

    T = _fmpz_vec_init(lenT + lenQ);
    Q = T + lenT;

    if (len1 >= len2)
        _fmpz_mod_poly_mul(T, poly1, len1, poly2, len2, ctx);
    else
        _fmpz_mod_poly_mul(T, poly2, len2, poly1, len1, ctx);

    fmpz_init(invf);
    fmpz_mod_inv(invf, f + lenf - 1, ctx);

    _fmpz_mod_poly_divrem(Q, res, T, lenT, f, lenf, invf, ctx);

    _fmpz_vec_clear(T, lenT + lenQ);
    fmpz_clear(invf);
}

void
fmpz_mod_poly_mulmod(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1,
               const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t f,
                                                      const fmpz_mod_ctx_t ctx)
{
    slong len1, len2, lenf;
    fmpz * fcoeffs;

    lenf = f->length;
    len1 = poly1->length;
    len2 = poly2->length;

    if (lenf == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mod_poly_mulmod). Divide by zero\n");
    }

    if (lenf == 1 || len1 == 0 || len2 == 0)
    {
        fmpz_mod_poly_zero(res, ctx);
        return;
    }

    if (len1 + len2 - lenf > 0)
    {
        if (f == res)
        {
            fcoeffs = _fmpz_vec_init(lenf);
            _fmpz_vec_set(fcoeffs, f->coeffs, lenf);
        }
        else
        {
            fcoeffs = f->coeffs;
        }

        fmpz_mod_poly_fit_length(res, len1 + len2 - 1, ctx);
        _fmpz_mod_poly_mulmod(res->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, fcoeffs, lenf, ctx);

        if (f == res)
            _fmpz_vec_clear(fcoeffs, lenf);

        _fmpz_mod_poly_set_length(res, lenf - 1);
        _fmpz_mod_poly_normalise(res);
    }
    else
    {
        fmpz_mod_poly_mul(res, poly1, poly2, ctx);
    }
}
