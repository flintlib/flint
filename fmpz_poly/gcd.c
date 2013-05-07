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
_fmpz_poly_gcd(fmpz * res, const fmpz * poly1, long len1,
               const fmpz * poly2, long len2)
{

    if (len1 < 6)
    {
        _fmpz_poly_gcd_subresultant(res, poly1, len1, poly2, len2);
    }
    else
    {
        long b1, b2;

        b1 = _fmpz_vec_max_bits(poly1, len1);
        b2 = _fmpz_vec_max_bits(poly2, len2);
        b1 = FLINT_ABS(b1);
        b2 = FLINT_ABS(b2);

        if (b1 + b2 < 2 * FLINT_BITS)
        {
            if (_fmpz_poly_gcd_heuristic(res, poly1, len1, poly2, len2))
                return;
        }

        _fmpz_poly_gcd_modular(res, poly1, len1, poly2, len2);
    }
}

void
fmpz_poly_gcd(fmpz_poly_t res, const fmpz_poly_t poly1,
              const fmpz_poly_t poly2)
{
    if (poly1->length < poly2->length)
    {
        fmpz_poly_gcd(res, poly2, poly1);
    }
    else /* len1 >= len2 >= 0 */
    {
        const long len1 = poly1->length;
        const long len2 = poly2->length;
        
        if (len1 == 0) /* len1 = len2 = 0 */
        {
            fmpz_poly_zero(res);
        } 
        else if (len2 == 0) /* len1 > len2 = 0 */
        {
            if (fmpz_sgn(poly1->coeffs + (len1 - 1)) > 0)
                fmpz_poly_set(res, poly1);
            else
                fmpz_poly_neg(res, poly1);
        }
        else /* len1 >= len2 >= 1 */
        {
            /* all current gcd functions automatically handle aliasing */
           
            fmpz_poly_fit_length(res, len2);
                
            _fmpz_poly_gcd(res->coeffs, poly1->coeffs, len1,
                                    poly2->coeffs, len2);
    
            _fmpz_poly_set_length(res, len2);
            _fmpz_poly_normalise(res);
        }
    }
}

