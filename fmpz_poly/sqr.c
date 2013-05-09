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

    Copyright (C) 2008, 2009 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void _fmpz_poly_sqr(fmpz * res, const fmpz * poly, len_t len)
{
    mp_size_t limbs;

    if (len < 7)
    {
        _fmpz_poly_sqr_classical(res, poly, len);
        return;
    }

    limbs = _fmpz_vec_max_limbs(poly, len);

    if (len < 16 && limbs > 12)
        _fmpz_poly_sqr_karatsuba(res, poly, len);
    else if (limbs <= 4)
        _fmpz_poly_sqr_KS(res, poly, len);
    else if (limbs/2048 > len)
        _fmpz_poly_sqr_KS(res, poly, len);
    else if (limbs*FLINT_BITS*4 < len)
       _fmpz_poly_sqr_KS(res, poly, len);
    else
       _fmpz_poly_mul_SS(res, poly, len, poly, len);
}

void fmpz_poly_sqr(fmpz_poly_t res, const fmpz_poly_t poly)
{
    len_t len = poly->length;
    len_t rlen;

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

