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

    Copyright (C) 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

long _fmpz_mod_poly_gcd_euclidean(fmpz *G, const fmpz *A, long lenA, 
                                           const fmpz *B, long lenB, 
                                           const fmpz_t invB, const fmpz_t p)
{
    if (lenB == 1)
    {
        fmpz_one(G);
        return 1;
    }
    else  /* lenA >= lenB > 1 */
    {
        fmpz_t invR3;
        fmpz *Q, *R1, *R2, *R3, *T, *W;
        long lenR1, lenR2, lenR3, lenW;

        lenW = FLINT_MAX(lenA - lenB + 1, lenB) + lenA + 2 * lenB;
        W    = _fmpz_vec_init(lenW);
        Q    = W;
        R1   = W + FLINT_MAX(lenA - lenB + 1, lenB);

        _fmpz_mod_poly_divrem(Q, R1, A, lenA, B, lenB, invB, p);

        lenR1 = lenB - 1;
        FMPZ_VEC_NORM(R1, lenR1);

        if (lenR1 == 0)
        {
            _fmpz_vec_set(G, B, lenB);
            _fmpz_vec_clear(W, lenW);
            return lenB;
        }

        fmpz_init(invR3);

        R2 = R1 + lenA;
        R3 = R2 + lenB;
        T  = R3;
        R3 = R1;
        R1 = T;
        _fmpz_vec_set(R2, B, lenB);
        lenR2 = lenB;
        lenR3 = lenR1;
        lenR1 = lenB;

        do
        {
            fmpz_invmod(invR3, R3 + (lenR3 - 1), p);

            _fmpz_mod_poly_divrem(Q, R1, R2, lenR2, R3, lenR3, invR3, p);
            lenR1 = lenR3 - 1;
            FMPZ_VEC_NORM(R1, lenR1);

            T  = R2;
            R2 = R3;
            R3 = R1;
            R1 = T;

            lenR2 = lenR3;
            lenR3 = lenR1;
        } 
        while (lenR3 > 0);

        _fmpz_vec_set(G, R2, lenR2);

        _fmpz_vec_clear(W, lenW);
        fmpz_clear(invR3);

        return lenR2;
    }
}

void fmpz_mod_poly_gcd_euclidean(fmpz_mod_poly_t G, 
                                 const fmpz_mod_poly_t A,
                                 const fmpz_mod_poly_t B)
{
    const long lenA = A->length;
    const long lenB = B->length;
    long lenG;
    fmpz *g;
    
    if (lenA == 0)
    {
        fmpz_mod_poly_make_monic(G, B);
        return;
    } 
    if (lenB == 0)
    {
        fmpz_mod_poly_make_monic(G, A);
        return;
    }

    if (lenA == 1 || lenB == 1)
    {
        fmpz_mod_poly_fit_length(G, 1);
        fmpz_one(G->coeffs);
        _fmpz_mod_poly_set_length(G, 1);
        return;
    }

    if (G == A || G == B)
    {
        g = _fmpz_vec_init(FLINT_MIN(lenA, lenB));
    }
    else
    {
        fmpz_mod_poly_fit_length(G, FLINT_MIN(lenA, lenB));
        g = G->coeffs;
    }

    if (lenA >= lenB)
    {
        fmpz_t invB;

        fmpz_init(invB);
        fmpz_invmod(invB, fmpz_mod_poly_lead(B), &(B->p));
        lenG = _fmpz_mod_poly_gcd_euclidean(g, A->coeffs, lenA,
                                               B->coeffs, lenB, invB, &(B->p));
        fmpz_clear(invB);
    }
    else
    {
        fmpz_t invA;

        fmpz_init(invA);
        fmpz_invmod(invA, fmpz_mod_poly_lead(A), &(A->p));
        lenG = _fmpz_mod_poly_gcd_euclidean(g, B->coeffs, lenB,
                                               A->coeffs, lenA, invA, &(A->p));
        fmpz_clear(invA);
    }

    if (G == A || G == B)
    {
        _fmpz_vec_clear(G->coeffs, G->alloc);
        G->coeffs = g;
        G->alloc  = FLINT_MIN(lenA, lenB);
        G->length = FLINT_MIN(lenA, lenB);
    }
    _fmpz_mod_poly_set_length(G, lenG);
    
    if (lenG == 1)
        fmpz_one(G->coeffs);
    else
        fmpz_mod_poly_make_monic(G, G);
}

