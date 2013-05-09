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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void _fmpq_poly_gcd(fmpz *G, fmpz_t denG, 
                    const fmpz *A, len_t lenA, const fmpz *B, len_t lenB)
{
    if (lenA == 1)  /* lenA == lenB == 1 */
    {
        fmpz_one(G);
        fmpz_one(denG);
    }
    else
    {
        fmpz *primA, *primB;
        fmpz_t s, t;
        len_t lenG;

        fmpz_init(s);
        fmpz_init(t);

        _fmpz_vec_content(s, A, lenA);
        _fmpz_vec_content(t, B, lenB);

        /* Set primA, primB to the primitive multiples of A, B */
        if (fmpz_is_one(s))
        {
            if (fmpz_is_one(t))
            {
                primA = (fmpz *) A;
                primB = (fmpz *) B;
            }
            else
            {
                primA = (fmpz *) A;
                primB = _fmpz_vec_init(lenB);
                _fmpz_vec_scalar_divexact_fmpz(primB, B, lenB, t);
            }
        }
        else
        {
            if (fmpz_is_one(s))
            {
                primA = _fmpz_vec_init(lenA);
                primB = (fmpz *) B;
                _fmpz_vec_scalar_divexact_fmpz(primA, A, lenA, s);
            }
            else
            {
                primA = _fmpz_vec_init(lenA + lenB);
                primB = primA + lenA;
                _fmpz_vec_scalar_divexact_fmpz(primA, A, lenA, s);
                _fmpz_vec_scalar_divexact_fmpz(primB, B, lenB, t);
            }
        }

        _fmpz_poly_gcd(G, primA, lenA, primB, lenB);

        for (lenG = lenB - 1; !G[lenG]; lenG--) ;
        lenG++;

        fmpz_set(denG, G + (lenG - 1));

        if (A != primA)
            _fmpz_vec_clear(primA, lenA + (B != primB) * lenB);
        else if (B != primB)
            _fmpz_vec_clear(primB, lenB);

        fmpz_clear(s);
        fmpz_clear(t);
    }
}

void fmpq_poly_gcd(fmpq_poly_t G, const fmpq_poly_t A, const fmpq_poly_t B)
{
    if (A->length < B->length)
    {
        fmpq_poly_gcd(G, B, A);
    }
    else
    {
        len_t lenA = A->length, lenB = B->length;

        if (lenA == 0) /* lenA = lenB = 0 */
        {
            fmpq_poly_zero(G);
        }
        else if (lenB == 0) /* lenA > lenB = 0 */
        {
            fmpq_poly_make_monic(G, A);
        }
        else /* lenA >= lenB >= 1 */
        {
            if (G == A || G == B)
            {
                fmpq_poly_t t;
                fmpq_poly_init2(t, lenB);
                
                _fmpq_poly_gcd(t->coeffs, t->den, A->coeffs, A->length, 
                                              B->coeffs, B->length);
            
                fmpq_poly_swap(t, G);
                fmpq_poly_clear(t);
            }
            else
            {
                fmpq_poly_fit_length(G, lenB);
                 
                _fmpq_poly_gcd(G->coeffs, G->den, A->coeffs, A->length, 
                                              B->coeffs, B->length);
            }

            _fmpq_poly_set_length(G, lenB);
            _fmpq_poly_normalise(G);
        }
    }
}

