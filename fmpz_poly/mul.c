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

void
_fmpz_poly_mul(fmpz * res, const fmpz * poly1, long len1,
               const fmpz * poly2, long len2)
{
    mp_size_t limbs1, limbs2;

    if (len2 < 7)
    {
        _fmpz_poly_mul_classical(res, poly1, len1, poly2, len2);
        return;
    }

    limbs1 = _fmpz_vec_max_limbs(poly1, len1);
    limbs2 = _fmpz_vec_max_limbs(poly2, len2);

    if (len1 < 16 && (limbs1 > 12 || limbs2 > 12))
        _fmpz_poly_mul_karatsuba(res, poly1, len1, poly2, len2);
    else if (limbs1 + limbs2 <= 8)
        _fmpz_poly_mul_KS(res, poly1, len1, poly2, len2);
    else if ((limbs1+limbs2)/2048 > len1 + len2)
        _fmpz_poly_mul_KS(res, poly1, len1, poly2, len2);
    else if ((limbs1 + limbs2)*FLINT_BITS*4 < len1 + len2)
       _fmpz_poly_mul_KS(res, poly1, len1, poly2, len2);
    else
       _fmpz_poly_mul_SS(res, poly1, len1, poly2, len2);
}

void
fmpz_poly_mul(fmpz_poly_t res,
              const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    long len1 = poly1->length;
    long len2 = poly2->length;
    long rlen;

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
        if (len1 >= len2)
            _fmpz_poly_mul(t->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2);
        else
            _fmpz_poly_mul(t->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
    else
    {
        fmpz_poly_fit_length(res, rlen);
        if (len1 >= len2)
            _fmpz_poly_mul(res->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2);
        else
            _fmpz_poly_mul(res->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1);
    }

    _fmpz_poly_set_length(res, rlen);
}
