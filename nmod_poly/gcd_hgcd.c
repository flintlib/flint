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
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"
#include "mpn_extras.h"

#define __set(B, lenB, A, lenA)      \
do {                                 \
    _nmod_vec_set((B), (A), (lenA)); \
    (lenB) = (lenA);                 \
} while (0)

#define __rem(R, lenR, A, lenA, B, lenB)                    \
do {                                                        \
    if ((lenA) >= (lenB))                                   \
    {                                                       \
        _nmod_poly_rem((R), (A), (lenA), (B), (lenB), mod); \
        (lenR) = (lenB) - 1;                                \
        MPN_NORM((R), (lenR));                              \
    }                                                       \
    else                                                    \
    {                                                       \
        _nmod_vec_set((R), (A), (lenA));                    \
        (lenR) = (lenA);                                    \
    }                                                       \
} while (0)

/*
    XXX: Incidentally, this implementation currently supports aliasing.  
    But since this may change in the future, no function other than 
    nmod_poly_gcd_hgcd() should rely on this.
 */

long _nmod_poly_gcd_hgcd(mp_ptr G, mp_srcptr A, long lenA, 
                                   mp_srcptr B, long lenB, nmod_t mod)
{
    const long cutoff = FLINT_BIT_COUNT(mod.n) <= 8 ? 
                        NMOD_POLY_SMALL_GCD_CUTOFF : NMOD_POLY_GCD_CUTOFF;

    mp_ptr J = _nmod_vec_init(2 * lenB);
    mp_ptr R = J + lenB;

    long lenG, lenJ, lenR;

    __rem(R, lenR, A, lenA, B, lenB);

    if (lenR == 0)
    {
        __set(G, lenG, B, lenB);
    }
    else
    {
        _nmod_poly_hgcd(NULL, NULL, G, &(lenG), J, &(lenJ), B, lenB, R, lenR, mod);

        while (lenJ != 0)
        {
            __rem(R, lenR, G, lenG, J, lenJ);

            if (lenR == 0)
            {
                __set(G, lenG, J, lenJ);
                break;
            }
            if (lenJ < cutoff)
            {
                lenG = _nmod_poly_gcd_euclidean(G, J, lenJ, R, lenR, mod);
                break;
            }

            _nmod_poly_hgcd(NULL, NULL, G, &(lenG), J, &(lenJ), J, lenJ, R, lenR, mod);
        }
    }
    _nmod_vec_clear(J);

    return lenG;
}

void nmod_poly_gcd_hgcd(nmod_poly_t G, const nmod_poly_t A, const nmod_poly_t B)
{
    const long lenA = A->length, lenB = B->length;

    if (lenA == 0)
    {
        if (lenB == 0) 
            nmod_poly_zero(G);
        else 
            nmod_poly_make_monic(G, B);
    }
    else if (lenB == 0)
    {
        nmod_poly_make_monic(G, A);
    }
    else
    {
        nmod_poly_fit_length(G, FLINT_MIN(lenA, lenB));

        if (lenA >= lenB)
        {
            G->length = _nmod_poly_gcd_hgcd(G->coeffs, A->coeffs, A->length, 
                                                       B->coeffs, B->length, A->mod);
        }
        else
        {
            G->length = _nmod_poly_gcd_hgcd(G->coeffs, B->coeffs, B->length, 
                                                       A->coeffs, A->length, A->mod);
        }

        _nmod_poly_make_monic(G->coeffs, G->coeffs, G->length, G->mod);
    }
}

#undef __set
#undef __rem

