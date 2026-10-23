/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"


void _fmpz_poly_sqr(fmpz * res, const fmpz * poly, slong len)
{
    if (len == 1)
    {
        fmpz_mul(res, poly, poly);
        return;
    }

    _fmpz_poly_mulmid(res, poly, len, poly, len, 0, 2 * len - 1);
}

void fmpz_poly_sqr(fmpz_poly_t res, const fmpz_poly_t poly)
{
    slong len = poly->length;
    slong rlen;

    if (len == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    rlen = 2 * len - 1;

    if (res == poly)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, rlen);
        _fmpz_poly_sqr(t->coeffs, poly->coeffs, len);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
    else
    {
        fmpz_poly_fit_length(res, rlen);
        _fmpz_poly_sqr(res->coeffs, poly->coeffs, len);
    }

    _fmpz_poly_set_length(res, rlen);
}
