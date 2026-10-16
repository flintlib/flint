/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

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

#if FLINT_HAVE_FFT_SMALL
#include "fft_small.h"
#endif

void
_fmpz_poly_mullow(fmpz * res, const fmpz * poly1, slong len1,
                                const fmpz * poly2, slong len2, slong n)
{
    if (len2 == 1)
    {
        _fmpz_vec_scalar_mul_fmpz(res, poly1, FLINT_MIN(len1, n), poly2);
        return;
    }
    if (len1 == 1)
    {
        _fmpz_vec_scalar_mul_fmpz(res, poly2, FLINT_MIN(len2, n), poly1);
        return;
    }

    _fmpz_poly_mulmid(res, poly1, len1, poly2, len2, 0, n);
}

void
fmpz_poly_mullow(fmpz_poly_t res,
                   const fmpz_poly_t poly1, const fmpz_poly_t poly2,
                   slong n)
{
    const slong len1 = poly1->length;
    const slong len2 = poly2->length;

    if (len1 == 0 || len2 == 0 || n == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    if (res == poly1 || res == poly2)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        fmpz_poly_mullow(t, poly1, poly2, n);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
        return;
    }

    n = FLINT_MIN(n, len1 + len2 - 1);
    fmpz_poly_fit_length(res, n);
    _fmpz_poly_mullow(res->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, n);
    _fmpz_poly_set_length(res, n);
    _fmpz_poly_normalise(res);
}
