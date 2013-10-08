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
    Copyright (C) 2012 Andres Goens

******************************************************************************/

#include <stdlib.h>
#include "fq_poly.h"

long _fq_poly_gcd_euclidean(fq_struct* G,const fq_struct* A, long lenA, 
                            const fq_struct* B, long lenB, const fq_ctx_t ctx)
{
    long lenR1 = lenA, lenR2 = lenB, lenT; /* assumes lenA \geq lenB */
    
    fq_struct *R1, *R2, *T;
    
    fq_t inv;

    fq_init(inv);
    T = _fq_poly_init(lenA);
    R1 = _fq_poly_init(lenA);
    R2 = _fq_poly_init(lenB);
    _fq_poly_set(R1,A,lenA);
    _fq_poly_set(R2,B,lenB);

    while(lenR2 > 0)
      {
          if(!fq_is_zero(&R2[lenR2-1])) fq_inv(inv,&R2[lenR2-1],ctx);
             else printf("something wrong with length");
	_fq_poly_rem(R1, R1, lenR1, R2, lenR2, inv, ctx);
	_fq_poly_normalise2(R1,&lenR1);
	
	_fq_poly_set(T,R1,lenR1);
	_fq_poly_set(R1,R2,lenR2);
	_fq_poly_set(R2,T,lenR1);
	lenT = lenR1; 
	lenR1 = lenR2;
	lenR2 = lenT;
      }

    _fq_poly_set(G,R1,lenR1);
    return lenR1;
}
/*

    fq_struct *F, *R1, *R2, *R3 = G, *T;
    fq_t inv;
    
    fq_init(inv);
    fq_inv(inv,&B[lenB-1],ctx);

    if (lenB == 1)
    {
        fq_set(&G[0],&B[0]);
        return 1;
    }

    F  = _fq_poly_init(2*lenB);
    R1 = F;
    R2 = R1 + lenB - 1;

    _fq_poly_rem(R1, A, lenA, B, lenB, inv, ctx);
    lenR1 = lenB - 1;
    _fq_poly_normalise2(R1,&lenR1);

    if (lenR1 > 1)
    {
        fq_inv(inv,R1+(lenR1-1),ctx);
        _fq_poly_rem(R2, B, lenB, R1, lenR1, inv, ctx);
        lenR2 = lenR1 - 1;
        _fq_poly_normalise2(R2, &lenR2);
    }
    else
    {
        if (lenR1 == 0)
        {
            _fq_poly_set(G, B, lenB);
            _fq_poly_clear(F,2*lenB);
            return lenB;
        }
        else
        {
            fq_set(&G[0],&R1[0]);
            _fq_poly_clear(F,2*lenB);
            return 1;
        }
    }

    for (steps = 2; lenR2 > 1; steps++)
    {

        fq_inv(inv,R2+(lenR2 -1),ctx);
        _fq_poly_normalise2(R2,&lenR2);
        _fq_poly_rem(R3, R1, lenR1, R2, lenR2, inv, ctx);
        lenR1 = lenR2--;
        _fq_poly_normalise2(R3, &lenR2);
        T = R1; R1 = R2; R2 = R3; R3 = T;
    }

    if (lenR2 == 1)
    {
        lenG = 1;
        if (steps % 3) 
            fq_set(&G[0],&R2[0]);
    }
    else
    {
        lenG = lenR1;
        if (steps % 3 != 1)
            _fq_poly_set(G, R1, lenR1);
    }

    _fq_poly_clear(F,2*lenB);
    fq_clear(inv);
    return lenG;
}

*/

void fq_poly_gcd_euclidean(fq_poly_t G, 
                           const fq_poly_t A, const fq_poly_t B, const fq_ctx_t ctx)
{
    if (A->length < B->length)
    {
        fq_poly_gcd_euclidean(G, B, A, ctx);
    }
    else /* lenA >= lenB >= 0 */
    {
        long lenA = A->length, lenB = B->length, lenG;
        fq_poly_t tG;
        fq_struct *g;

        if (lenA == 0) /* lenA = lenB = 0 */
        {
            fq_poly_zero(G);
        } 
        else if (lenB == 0) /* lenA > lenB = 0 */
        {
            fq_poly_make_monic(G, A, ctx);
        }
        else /* lenA >= lenB >= 1 */
        {
            if (G == A || G == B)
            {
                fq_poly_init2(tG, FLINT_MIN(lenA, lenB));
                g = tG->coeffs;
            }
            else
            {
                fq_poly_fit_length(G, FLINT_MIN(lenA, lenB));
                g = G->coeffs;
            }

            lenG = _fq_poly_gcd_euclidean(g, A->coeffs, lenA,
                                               B->coeffs, lenB, ctx);

            if (G == A || G == B)
            {
                fq_poly_swap(tG, G);
                fq_poly_clear(tG);
            }
            G->length = lenG;

            if (G->length == 1)
                fq_one(&(G->coeffs[0]));
            else
                fq_poly_make_monic(G, G,ctx);
        }
    }
}
