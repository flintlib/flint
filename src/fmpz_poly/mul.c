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


void
_fmpz_poly_mul(fmpz * res, const fmpz * poly1,
                         slong len1, const fmpz * poly2, slong len2)
{
    /* multiplication by a scalar: avoid the dispatch overhead entirely
       (Z[x] operations on scalar values are common) */
    if (len2 == 1)
    {
        _fmpz_vec_scalar_mul_fmpz(res, poly1, len1, poly2);
        return;
    }
    if (len1 == 1)
    {
        _fmpz_vec_scalar_mul_fmpz(res, poly2, len2, poly1);
        return;
    }

    /* all the remaining algorithm selection happens in _fmpz_poly_mulmid */
    _fmpz_poly_mulmid(res, poly1, len1, poly2, len2, 0, len1 + len2 - 1);
}

void
fmpz_poly_mul(fmpz_poly_t res,
              const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong rlen;

    if (len1 == 0 || len2 == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    rlen = len1 + len2 - 1;

    if (res == poly1 || res == poly2)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, rlen);
        _fmpz_poly_mul(t->coeffs, poly1->coeffs, len1, poly2->coeffs, len2);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
    else
    {
        fmpz_poly_fit_length(res, rlen);
        _fmpz_poly_mul(res->coeffs, poly1->coeffs, len1, poly2->coeffs, len2);
    }

    _fmpz_poly_set_length(res, rlen);
}
