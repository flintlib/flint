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
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_gcd_subresultant(fmpz * res, const fmpz * poly1, len_t len1, 
                                        const fmpz * poly2, len_t len2)
{
    if (len2 == 1)
    {
        fmpz_t c;
        fmpz_init(c);
        _fmpz_poly_content(c, poly1, len1);
        fmpz_gcd(res, c, poly2);
        fmpz_clear(c);
    }
    else
    {
        fmpz_t a, b, d, g, h;
        fmpz *A, *B, *W;
        len_t lenA, lenB;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(d);
        fmpz_init(g);
        fmpz_init(h);

        A = W = _fmpz_vec_init(len1 + len2);
        B = W + len1;

        lenA = len1;
        lenB = len2;
        _fmpz_poly_content(a, poly1, lenA);
        _fmpz_poly_content(b, poly2, lenB);
        _fmpz_vec_scalar_divexact_fmpz(A, poly1, lenA, a);
        _fmpz_vec_scalar_divexact_fmpz(B, poly2, lenB, b);

        fmpz_gcd(d, a, b);
        fmpz_one(g);
        fmpz_one(h);

        while (1)
        {
            const len_t delta = lenA - lenB;

            _fmpz_poly_pseudo_rem_cohen(A, A, lenA, B, lenB);

            FMPZ_VEC_NORM(A, lenA);

            if (lenA <= 1)
                break;

            {                       /* Swap A and B */
                fmpz *T;
                len_t len;
                T = A, A = B, B = T, len = lenA, lenA = lenB, lenB = len;
            }

            if (delta == 1)
            {
                fmpz_mul(b, g, h);
                _fmpz_vec_scalar_divexact_fmpz(B, B, lenB, b);
                fmpz_set(g, A + (lenA - 1));
                fmpz_set(h, g);
            }
            else
            {
                fmpz_pow_ui(a, h, delta);
                fmpz_mul(b, g, a);
                _fmpz_vec_scalar_divexact_fmpz(B, B, lenB, b);
                fmpz_pow_ui(b, A + (lenA - 1), delta);
                fmpz_mul(g, h, b);
                fmpz_divexact(h, g, a);
                fmpz_set(g, A + (lenA - 1));
            }
        }

        if (lenA == 1)
        {
            fmpz_set(res, d);
            _fmpz_vec_zero(res + 1, len2 - 1);
        }
        else
        {
            _fmpz_poly_content(b, B, lenB);
            _fmpz_vec_scalar_divexact_fmpz(B, B, lenB, b);
            if (fmpz_sgn(B + (lenB - 1)) < 0)
                fmpz_neg(d, d);
            _fmpz_vec_scalar_mul_fmpz(res, B, lenB, d);
            if (len2 >= lenB)
                _fmpz_vec_zero(res + lenB, len2 - lenB);
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(d);
        fmpz_clear(g);
        fmpz_clear(h);

        _fmpz_vec_clear(W, len1 + len2);
    }
}

void
fmpz_poly_gcd_subresultant(fmpz_poly_t res, const fmpz_poly_t poly1,
              const fmpz_poly_t poly2)
{
    if (poly1->length < poly2->length)
    {
        fmpz_poly_gcd_subresultant(res, poly2, poly1);
    }
    else /* len1 >= len2 >= 0 */
    {
        const len_t len1 = poly1->length;
        const len_t len2 = poly2->length;
        
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
            /* underscore code automatically handles aliasing */
           
            fmpz_poly_fit_length(res, len2);
                
            _fmpz_poly_gcd_subresultant(res->coeffs, poly1->coeffs, len1,
                                    poly2->coeffs, len2);
    
            _fmpz_poly_set_length(res, len2);
            _fmpz_poly_normalise(res);
        }
    }
}
