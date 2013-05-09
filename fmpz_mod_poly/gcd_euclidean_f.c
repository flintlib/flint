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
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

len_t _fmpz_mod_poly_gcd_euclidean_f(fmpz_t f, fmpz *G, 
                                    const fmpz *A, len_t lenA, 
                                    const fmpz *B, len_t lenB, const fmpz_t p)
{
    len_t lenG = 0;

    if (lenB == 1)
    {
        fmpz_t invB;
        fmpz_init(invB);
        fmpz_gcdinv(f, invB, B, p);
        if (fmpz_is_one(f))
        {
            fmpz_one(G);
            lenG = 1;
        }
        fmpz_clear(invB);
    }
    else  /* lenA >= lenB > 1 */
    {
        const len_t lenW = FLINT_MAX(lenA - lenB + 1, lenB) + lenA + 2 * lenB;
        fmpz *Q, *R1, *R2, *R3, *T, *W;
        len_t lenR2, lenR3;

        W  = _fmpz_vec_init(lenW);
        Q  = W;
        R1 = W + FLINT_MAX(lenA - lenB + 1, lenB);
        R2 = R1 + lenA;
        R3 = R2 + lenB;

        _fmpz_mod_poly_divrem_f(f, Q, R1, A, lenA, B, lenB, p);
        if (!fmpz_is_one(f))
            goto exit;

        lenR3 = lenB - 1;
        FMPZ_VEC_NORM(R1, lenR3);

        if (lenR3 == 0)
        {
            _fmpz_vec_set(G, B, lenB);
            lenG = lenB;
        }
        else
        {
            T  = R3;
            R3 = R1;
            R1 = T;
            _fmpz_vec_set(R2, B, lenB);
            lenR2 = lenB;

            do
            {
                _fmpz_mod_poly_divrem_f(f, Q, R1, R2, lenR2, R3, lenR3, p);
                if (!fmpz_is_one(f))
                    goto exit;

                lenR2 = lenR3--;
                FMPZ_VEC_NORM(R1, lenR3);
                T = R2; R2 = R3; R3 = R1; R1 = T;
            } 
            while (lenR3 > 0);

            _fmpz_vec_set(G, R2, lenR2);
            lenG = lenR2;
        }

      exit:
        _fmpz_vec_clear(W, lenW);
    }
    return lenG;
}

void fmpz_mod_poly_gcd_euclidean_f(fmpz_t f, fmpz_mod_poly_t G, 
                                   const fmpz_mod_poly_t A,
                                   const fmpz_mod_poly_t B)
{
    if (A->length < B->length)
    {
        fmpz_mod_poly_gcd_euclidean_f(f, G, B, A);
    }
    else /* lenA >= lenB >= 0 */
    {
        const len_t lenA = A->length, lenB = B->length;
        len_t lenG;
        fmpz *g;

        if (lenA == 0) /* lenA = lenB = 0 */
        {
            fmpz_mod_poly_zero(G);
            fmpz_one(f);
        }
        else if (lenB == 0) /* lenA > lenB = 0 */
        {
            fmpz_t invA;
            fmpz_init(invA);
            fmpz_gcdinv(f, invA, A->coeffs + lenA - 1, &B->p);
            if (fmpz_is_one(f))
                fmpz_mod_poly_scalar_mul_fmpz(G, A, invA);
            else
                fmpz_mod_poly_zero(G);
            fmpz_clear(invA);
        }
        else /* lenA >= lenB >= 1 */
        {
            if (G == A || G == B)
            {
                g = _fmpz_vec_init(FLINT_MIN(lenA, lenB));
            }
            else
            {
                fmpz_mod_poly_fit_length(G, FLINT_MIN(lenA, lenB));
                g = G->coeffs;
            }

            lenG = _fmpz_mod_poly_gcd_euclidean_f(f, g, A->coeffs, lenA,
                                                    B->coeffs, lenB, &(B->p));

            if (fmpz_is_one(f))
            {
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
            else  /* Factor found, ensure G is normalised */
            {
                if (G == A || G == B)
                    _fmpz_vec_clear(g, FLINT_MIN(lenA, lenB));
                else
                {
                    _fmpz_vec_zero(G->coeffs, FLINT_MIN(lenA, lenB));
                    _fmpz_mod_poly_set_length(G, 0);
                }
            }
        }
    }
}

