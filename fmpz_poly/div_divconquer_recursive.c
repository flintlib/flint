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

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

#define FLINT_DIV_DIVCONQUER_CUTOFF  16

void
_fmpz_poly_div_divconquer_recursive(fmpz * Q, fmpz * temp, 
                                    const fmpz * A, const fmpz * B, long lenB)
{
    if (lenB <= FLINT_DIV_DIVCONQUER_CUTOFF)
    {
        _fmpz_poly_div_basecase(Q, temp, A, 2 * lenB - 1, B, lenB);
    }
    else
    {
        const long n2 = lenB / 2;
        const long n1 = lenB - n2;

        fmpz * q0 = Q;
        fmpz * q1 = Q + n2;

        /*
           t is a vector of length lenB - 1, h points to the top n2 coeffs 
           of t;  r1 is vector of length lenB >= 2 n1 - 1
         */

        fmpz * t  = temp;
        fmpz * h  = temp + (n1 - 1);
        fmpz * r1 = temp + (lenB - 1);

        /*
           Set {q1, n1}, {r1, 2 n1 - 1} to the quotient and remainder of 
           {A + 2 n2, 2 n1 - 1} divided by {B + n2, n1}
         */

        _fmpz_poly_divremlow_divconquer_recursive(q1, r1, A + 2 * n2, B + n2, n1);
        _fmpz_vec_sub(r1, A + 2 * n2, r1, n1 - 1);

        /*
           Set the top n2 coeffs of t to the top n2 coeffs of the product of 
           {q1, n1} and {B, n2}; the bottom n1 - 1 coeffs may be arbitrary

           For sufficiently large polynomials, computing the full product 
           using Kronecker segmentation is faster than computing the opposite 
           short product via Karatsuba
         */

        _fmpz_poly_mul_KS(t, q1, n1, B, n2);

        /*
           If lenB is odd, set {h, n2} to {r1, n2} - {h, n2}, otherwise, to 

               {A + lenB - 1, 1} + {x * r1, n2} - {h, n2}
         */

        if (lenB & 1L)
        {
            _fmpz_vec_sub(h, r1, h, n2);
        }
        else
        {
            _fmpz_vec_sub(h + 1, r1, h + 1, n2 - 1);
            fmpz_neg(h, h);
            fmpz_add(h, h, A + lenB - 1);
        }

        /*
           Set t to h shifted to the right by n2 - 1, and set q0 to the 
           quotient of {t, 2 n2 - 1} and {B + n1, n2}
           
           Note the bottom n2 - 1 coefficients of t are irrelevant
         */

        t += (lenB & 1L);
        
        _fmpz_poly_div_divconquer_recursive(q0, temp + lenB, t, B + n1, n2);
    }
}

