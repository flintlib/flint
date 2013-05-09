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

******************************************************************************/

#include <stdlib.h>
#include "nmod_poly.h"
#include "mpn_extras.h"

len_t _nmod_poly_gcd_euclidean(mp_ptr G, mp_srcptr A, len_t lenA, 
                                        mp_srcptr B, len_t lenB, nmod_t mod)
{
    len_t steps;
    len_t lenR1, lenR2 = 0, lenG = 0;

    mp_ptr F, R1, R2, R3 = G, T;
    
    if (lenB == 1)
    {
        G[0] = B[0];
        return 1;
    }

    F  = _nmod_vec_init(2*lenB - 3);
    R1 = F;
    R2 = R1 + lenB - 1;

    _nmod_poly_rem(R1, A, lenA, B, lenB, mod);
    lenR1 = lenB - 1;
    MPN_NORM(R1, lenR1);

    if (lenR1 > 1)
    {
        _nmod_poly_rem(R2, B, lenB, R1, lenR1, mod);
        lenR2 = lenR1 - 1;
        MPN_NORM(R2, lenR2);
    }
    else
    {
        if (lenR1 == 0)
        {
            flint_mpn_copyi(G, B, lenB);
            _nmod_vec_clear(F);
            return lenB;
        }
        else
        {
            G[0] = R1[0];
            _nmod_vec_clear(F);
            return 1;
        }
    }

    for (steps = 2; lenR2 > 1; steps++)
    {
        _nmod_poly_rem(R3, R1, lenR1, R2, lenR2, mod);
        lenR1 = lenR2--;
        MPN_NORM(R3, lenR2);
        T = R1; R1 = R2; R2 = R3; R3 = T;
    }

    if (lenR2 == 1)
    {
        lenG = 1;
        if (steps % 3) 
            G[0] = R2[0];
    }
    else
    {
        lenG = lenR1;
        if (steps % 3 != 1)
            flint_mpn_copyi(G, R1, lenR1);
    }

    _nmod_vec_clear(F);
    return lenG;
}

void nmod_poly_gcd_euclidean(nmod_poly_t G, 
                             const nmod_poly_t A, const nmod_poly_t B)
{
    if (A->length < B->length)
    {
        nmod_poly_gcd_euclidean(G, B, A);
    }
    else /* lenA >= lenB >= 0 */
    {
        len_t lenA = A->length, lenB = B->length, lenG;
        nmod_poly_t tG;
        mp_ptr g;

        if (lenA == 0) /* lenA = lenB = 0 */
        {
            nmod_poly_zero(G);
        } 
        else if (lenB == 0) /* lenA > lenB = 0 */
        {
            nmod_poly_make_monic(G, A);
        }
        else /* lenA >= lenB >= 1 */
        {
            if (G == A || G == B)
            {
                nmod_poly_init2(tG, A->mod.n, FLINT_MIN(lenA, lenB));
                g = tG->coeffs;
            }
            else
            {
                nmod_poly_fit_length(G, FLINT_MIN(lenA, lenB));
                g = G->coeffs;
            }

            lenG = _nmod_poly_gcd_euclidean(g, A->coeffs, lenA,
                                               B->coeffs, lenB, A->mod);

            if (G == A || G == B)
            {
                nmod_poly_swap(tG, G);
                nmod_poly_clear(tG);
            }
            G->length = lenG;

            if (G->length == 1)
                G->coeffs[0] = 1;
            else
                nmod_poly_make_monic(G, G);
        }
    }
}
