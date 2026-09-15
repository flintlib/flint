/*
    Copyright (C) 2008, 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
fmpz_poly_mulhigh_n(fmpz_poly_t res, const fmpz_poly_t poly1,
                    const fmpz_poly_t poly2, slong n)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong lenfull, lo;

    if (n == 0 || len1 == 0 || len2 == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    if (res == poly1 || res == poly2)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, 2 * n - 1);
        fmpz_poly_mulhigh_n(t, poly1, poly2, n);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
        return;
    }

    /* the top n coefficients of the product (which has lenfull
       coefficients; the top coefficient is nonzero); anything below is
       left arbitrary, matching the documented interface. The operands
       are multiplied at their full lengths even beyond n, as some
       callers (fmpz_poly_divhigh_smodp) rely on this. */
    lenfull = len1 + len2 - 1;
    lo = FLINT_MAX(lenfull - n, 0);

    fmpz_poly_fit_length(res, lenfull);
    _fmpz_poly_mulmid(res->coeffs + lo, poly1->coeffs, len1,
                      poly2->coeffs, len2, lo, lenfull);
    _fmpz_poly_set_length(res, lenfull);
}
