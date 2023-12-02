/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void
_fmpz_poly_CRT_ui_precomp(fmpz * res, const fmpz * poly1, slong len1,
               const fmpz_t m1, mp_srcptr poly2, slong len2, mp_limb_t m2,
                mp_limb_t m2inv, fmpz_t m1m2, mp_limb_t c, int sign)
{
    slong i;

    for (i = 0; i < FLINT_MIN(len1, len2); i++)
    {
        _fmpz_CRT_ui_precomp(res + i, poly1 + i, m1,
                                poly2[i], m2, m2inv, m1m2, c, sign);
    }

    if (len2 > len1)
    {
        fmpz_t zero;
        fmpz_init(zero);
        for (i = len1; i < len2; i++)
        {
            _fmpz_CRT_ui_precomp(res + i, zero, m1,
                                    poly2[i], m2, m2inv, m1m2, c, sign);
        }
        fmpz_clear(zero);
    }

    for (i = len2; i < len1; i++)
    {
        _fmpz_CRT_ui_precomp(res + i, res + i, m1,
                            0, m2, m2inv, m1m2, c, sign);
    }
}

void
_fmpz_poly_CRT_ui(fmpz * res, const fmpz * poly1, slong len1,
               const fmpz_t m1, mp_srcptr poly2, slong len2, mp_limb_t m2,
                mp_limb_t m2inv, int sign)
{
    mp_limb_t c;
    fmpz_t m1m2;

    c = fmpz_fdiv_ui(m1, m2);
    c = n_invmod(c, m2);

    if (c == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (_fmpz_poly_CRT_ui): m1 not invertible modulo m2.\n");
    }

    fmpz_init(m1m2);
    fmpz_mul_ui(m1m2, m1, m2);

    _fmpz_poly_CRT_ui_precomp(res, poly1, len1, m1,
                                   poly2, len2, m2, m2inv, m1m2, c, sign);

    fmpz_clear(m1m2);
}

void
fmpz_poly_CRT_ui(fmpz_poly_t res, const fmpz_poly_t poly1,
                        const fmpz_t m1, const nmod_poly_t poly2, int sign)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len = FLINT_MAX(len1, len2);

    if (len == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    fmpz_poly_fit_length(res, len);

    _fmpz_poly_CRT_ui(res->coeffs, poly1->coeffs, poly1->length, m1,
        poly2->coeffs, poly2->length, poly2->mod.n, poly2->mod.ninv, sign);

    _fmpz_poly_set_length(res, len);
    _fmpz_poly_normalise(res);
}
