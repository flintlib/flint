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

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_gcd(fmpz * res, const fmpz * poly1, long len1,
               const fmpz * poly2, long len2)
{
    /* if (!_fmpz_poly_gcd_heuristic(res, poly1, len1, poly2, len2)) */
    _fmpz_poly_gcd_subresultant(res, poly1, len1, poly2, len2);
}

void
fmpz_poly_gcd(fmpz_poly_t res, const fmpz_poly_t poly1,
              const fmpz_poly_t poly2)
{
    const long len1 = poly1->length;
    const long len2 = poly2->length;
    long rlen;
    
    if (len1 == 0)
    {
        if (len2 == 0)
            fmpz_poly_zero(res);
        else
        {
            if (fmpz_sgn(poly2->coeffs + (len2 - 1)) > 0)
                fmpz_poly_set(res, poly2);
            else
                fmpz_poly_neg(res, poly2);
        }
        return;
    }
    else
    {
        if (len2 == 0)
        {
            if (fmpz_sgn(poly1->coeffs + (len1 - 1)) > 0)
                fmpz_poly_set(res, poly1);
            else
                fmpz_poly_neg(res, poly1);
            return;
        }
    }

    rlen = FLINT_MIN(len1, len2);

    if (res == poly1 || res == poly2)
    {
       fmpz_poly_t temp;
       fmpz_poly_init2(temp, rlen);
       if (len1 >= len2)
          _fmpz_poly_gcd(temp->coeffs, poly1->coeffs, len1,
                                    poly2->coeffs, len2);
       else
          _fmpz_poly_gcd(temp->coeffs, poly2->coeffs, len2,
                                    poly1->coeffs, len1);
       fmpz_poly_swap(temp, res);
       fmpz_poly_clear(temp);
    }
    else
    {
       fmpz_poly_fit_length(res, rlen);
       if (len1 >= len2)
          _fmpz_poly_gcd(res->coeffs, poly1->coeffs, len1,
                                    poly2->coeffs, len2);
       else
          _fmpz_poly_gcd(res->coeffs, poly2->coeffs, len2,
                                    poly1->coeffs, len1);
    }
    
    _fmpz_poly_set_length(res, rlen);
    _fmpz_poly_normalise(res);
    
    return;
}

