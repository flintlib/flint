/*
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2018 William Hart
    
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

int
_fmpz_poly_sqrtrem_classical(fmpz * res, fmpz * r, const fmpz * poly, slong len)
{
    slong i, m;
    int result;

    /* the degree must be even */
    if (len % 2 == 0)
        return 0;

    /* check whether a square root exists modulo 2 */
    m = (len + 1) / 2;
    
    for (i = ((m - 1) | 1); i < len; i += 2)
        if (!fmpz_is_even(poly + i))
            return 0;

    /* check leading coefficient */
    if (!fmpz_is_square(poly + len - 1))
        return 0;

    /* square root of leading coefficient */
    fmpz_sqrt(res + m - 1, poly + len - 1);
    result = 1;

    /* do divison style 'square root with remainder' from top to bottom */
    if (len > 1)
    {
        fmpz_t t, u;

        fmpz_init(t);
        fmpz_init(u);
        
        if (r != poly)
           _fmpz_vec_set(r, poly, len);
        
        fmpz_mul_ui(u, res + m - 1, 2);

        for (i = 1; i < m; i++)
        {
            fmpz_fdiv_qr(res + m - i - 1, t, r + len - i - 1, u);
            if (!fmpz_is_zero(t))
            {
                result = 0;
                break;
            }

            fmpz_mul_si(t, res + m - i - 1, -2);
            _fmpz_vec_scalar_addmul_fmpz(r + len - 2*i, res + m - i, i - 1, t);
            fmpz_submul(r + len - 2*i - 1, res + m - i - 1, res + m - i - 1);
        }

        fmpz_clear(t);
        fmpz_clear(u);
    }

    return result;
}

int
fmpz_poly_sqrtrem_classical(fmpz_poly_t b, fmpz_poly_t r, const fmpz_poly_t a)
{
    slong blen, len = a->length;
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
        result = fmpz_poly_sqrtrem_classical(tmp, r, a);
        fmpz_poly_swap(b, tmp);
        fmpz_poly_clear(tmp);
        return result;
    }

    blen = len / 2 + 1;
    fmpz_poly_fit_length(r, len);
    fmpz_poly_fit_length(b, blen);
    _fmpz_poly_set_length(b, blen);
    result = _fmpz_poly_sqrtrem_classical(b->coeffs, r->coeffs, a->coeffs, len);
    if (!result)
        _fmpz_poly_set_length(b, 0);
    else
    {
       _fmpz_poly_set_length(r, blen - 1);
       _fmpz_poly_normalise(r);
    }
    return result;
}
