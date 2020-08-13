/*
    Copyright (C) 2008-2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "fft.h"
#include "fft_tuning.h"

void _fmpz_poly_mul_SS(fmpz *output, const fmpz *input1, slong len1, 
                       const fmpz *input2, slong len2)
{
    const slong rlen = len1 + len2 - 1;

    _fmpz_poly_mullow_SS(output, input1, len1, input2, len2, rlen);
}

void
fmpz_poly_mul_SS(fmpz_poly_t res,
                 const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    const slong len1 = poly1->length, len2 = poly2->length;
    slong rlen;

    if (len1 == 0 || len2 == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    if (len1 <= 2 || len2 <= 2)
    {
        fmpz_poly_mul_classical(res, poly1, poly2);
        return;
    }

    rlen = len1 + len2 - 1;

    fmpz_poly_fit_length(res, rlen);
    if (len1 >= len2)
        _fmpz_poly_mul_SS(res->coeffs, poly1->coeffs, len1,
            poly2->coeffs, len2);
    else
        _fmpz_poly_mul_SS(res->coeffs, poly2->coeffs, len2,
                          poly1->coeffs, len1);
    _fmpz_poly_set_length(res, rlen);
}
