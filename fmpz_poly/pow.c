/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_pow(fmpz * res, const fmpz * poly, long len, ulong e)
{
    if (e < 5UL)
        _fmpz_poly_pow_small(res, poly, len, e);
    else if (len == 2)
        _fmpz_poly_pow_binomial(res, poly, e);
    else
    {
        ulong limbs = (ulong) _fmpz_vec_max_limbs(poly, len);
        
        if (limbs < ((3UL * e) / 2UL + 150UL) / (ulong) len)
            _fmpz_poly_pow_multinomial(res, poly, len, e);
        else
            _fmpz_poly_pow_binexp(res, poly, len, e);
    }
}

void
fmpz_poly_pow(fmpz_poly_t res, const fmpz_poly_t poly, ulong e)
{
    const long len = poly->length;
    long rlen;

    if ((len < 2) | (e < 3UL))
    {
        if (e == 0UL)
            fmpz_poly_set_ui(res, 1);
        else if (len == 0)
            fmpz_poly_zero(res);
        else if (len == 1)
        {
            fmpz_poly_fit_length(res, 1);
            fmpz_pow_ui(res->coeffs, poly->coeffs, e);
            _fmpz_poly_set_length(res, 1);
        }
        else if (e == 1UL)
            fmpz_poly_set(res, poly);
        else  /* e == 2UL */
            fmpz_poly_sqr(res, poly);
        return;
    }

    rlen = (long) e * (len - 1) + 1;

    if (res != poly)
    {
        fmpz_poly_fit_length(res, rlen);
        _fmpz_poly_set_length(res, rlen);
        _fmpz_poly_pow(res->coeffs, poly->coeffs, len, e);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, rlen);
        _fmpz_poly_set_length(t, rlen);
        _fmpz_poly_pow(t->coeffs, poly->coeffs, len, e);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
}
