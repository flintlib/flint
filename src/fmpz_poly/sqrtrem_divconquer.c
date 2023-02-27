/*
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2018, 2019 William Hart
    
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

/*
    Implements the square root with remainder function of "speeding up the
    division and square root of power series", by Hanrot, Quercia and
    Zimmermann (see https://hal.inria.fr/inria-00072675/document).

    TODO: implement middle product so that we can implement their Sqrt function
    (their specialised square root Divide function requires it).
*/
int
_fmpz_poly_sqrtrem_divconquer(fmpz * res, fmpz * r, const fmpz * poly, slong len, fmpz * temp)
{
    slong i, n, n2;

    if (len < FMPZ_POLY_SQRTREM_DIVCONQUER_CUTOFF)
        return _fmpz_poly_sqrtrem_classical(res, r, poly, len);

    /* the degree must be even */
    if (len % 2 == 0)
        return 0;

    n = (len + 1)/2;

    /* check whether a square root exists modulo 2 */
    n2 = (n + 1)/2;

    /* only check coeffs that won't be checked recursively */
    for (i = ((n - 1) | 1); i < len - n2; i += 2)
        if (!fmpz_is_even(poly + i))
            return 0;

    if (r != poly)
        _fmpz_vec_set(r, poly, len);
        
    if (!_fmpz_poly_sqrtrem_divconquer(res + n - n2, r + len - 2*n2 + 1, r + len - 2*n2 + 1, 2*n2 - 1, temp))
        return 0;

    _fmpz_vec_scalar_mul_ui(temp, res + n - n2, n2, 2);

   _fmpz_vec_set(temp + n, r + n2, 2*n - 2*n2 - 1);

    if (!_fmpz_poly_divrem(res, r + n2, temp + n, 2*n - 2*n2 - 1, temp + 2*n2 - n, n - n2, 1))
       return 0;

    _fmpz_poly_mul(temp + 2*n2 - n, res, n - n2, res, n - n2);

    _fmpz_vec_sub(r, r, temp + 2*n2 - n, 2*n - 2*n2 - 1);

    if (2*n2 > n)
        _fmpz_vec_scalar_submul_fmpz(r + n - n2, res, n2 - 1, temp);

    return 1;
}

int
fmpz_poly_sqrtrem_divconquer(fmpz_poly_t b, fmpz_poly_t r, const fmpz_poly_t a)
{
    slong blen, len = a->length;
    fmpz * temp;
    int result;

    if (len % 2 == 0)
    {
        fmpz_poly_zero(b);
        fmpz_poly_zero(r);
        return len == 0;
    }
    
    if (b == a)
    {
        fmpz_poly_t tmp;
        fmpz_poly_init(tmp);
        result = fmpz_poly_sqrtrem_divconquer(tmp, r, a);
        fmpz_poly_swap(b, tmp);
        fmpz_poly_clear(tmp);
        return result;
    }

    blen = len / 2 + 1;
    fmpz_poly_fit_length(r, len);
    fmpz_poly_fit_length(b, blen);
    _fmpz_poly_set_length(b, blen);
    temp = _fmpz_vec_init(len);
    result = _fmpz_poly_sqrtrem_divconquer(b->coeffs, r->coeffs, a->coeffs, len, temp);
    if (!result)
        _fmpz_poly_set_length(b, 0);
    else
    {
       _fmpz_poly_set_length(r, blen - 1);
       _fmpz_poly_normalise(r);
    }

    _fmpz_vec_clear(temp, len);

    return result;
}
