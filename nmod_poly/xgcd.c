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
    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include "nmod_poly.h"
#include "mpn_extras.h"

len_t _nmod_poly_xgcd(mp_ptr G, mp_ptr S, mp_ptr T, 
                     mp_srcptr A, len_t lenA, mp_srcptr B, len_t lenB, nmod_t mod)
{
    const len_t cutoff = FLINT_BIT_COUNT(mod.n) <= 8 ? 
                        NMOD_POLY_SMALL_GCD_CUTOFF : NMOD_POLY_GCD_CUTOFF;

    if (lenA < cutoff)
        return _nmod_poly_xgcd_euclidean(G, S, T, A, lenA, B, lenB, mod);
    else
        return _nmod_poly_xgcd_hgcd(G, S, T, A, lenA, B, lenB, mod);
}

void
nmod_poly_xgcd(nmod_poly_t G, nmod_poly_t S, nmod_poly_t T,
                         const nmod_poly_t A, const nmod_poly_t B)
{
    if (A->length < B->length)
    {
        nmod_poly_xgcd(G, T, S, B, A);
    }
    else  /* lenA >= lenB >= 0 */
    {
        const len_t lenA = A->length, lenB = B->length;
        mp_limb_t inv;

        if (lenA == 0)  /* lenA = lenB = 0 */
        {
            nmod_poly_zero(G);
            nmod_poly_zero(S);
            nmod_poly_zero(T);
        }
        else if (lenB == 0)  /* lenA > lenB = 0 */
        {
            inv = n_invmod(A->coeffs[lenA - 1], A->mod.n);
            nmod_poly_scalar_mul_nmod(G, A, inv);
            nmod_poly_zero(T);
            nmod_poly_set_coeff_ui(S, 0, inv);
            S->length = 1;
        }
        else if (lenB == 1)  /* lenA >= lenB = 1 */
        {
            nmod_poly_fit_length(T, 1);
            T->length = 1;
            T->coeffs[0] = n_invmod(B->coeffs[0], A->mod.n);
            nmod_poly_one(G);
            nmod_poly_zero(S);
        }
        else  /* lenA >= lenB >= 2 */
        {
            mp_ptr g, s, t;
            len_t lenG;

            if (G == A || G == B)
            {
                g = _nmod_vec_init(FLINT_MIN(lenA, lenB));
            }
            else
            {
                nmod_poly_fit_length(G, FLINT_MIN(lenA, lenB));
                g = G->coeffs;
            }
            if (S == A || S == B)
            {
                s = _nmod_vec_init(lenB - 1);
            }
            else
            {
                nmod_poly_fit_length(S, lenB - 1);
                s = S->coeffs;
            }
            if (T == A || T == B)
            {
                t = _nmod_vec_init(lenA - 1);
            }
            else
            {
                nmod_poly_fit_length(T, lenA - 1);
                t = T->coeffs;
            }

            if (lenA >= lenB)
                lenG = _nmod_poly_xgcd(g, s, t, A->coeffs, lenA,
                                                          B->coeffs, lenB, A->mod);
            else
                lenG = _nmod_poly_xgcd(g, t, s, B->coeffs, lenB,
                                                          A->coeffs, lenA, A->mod);

            if (G == A || G == B)
            {
                free(G->coeffs);
                G->coeffs = g;
                G->alloc  = FLINT_MIN(lenA, lenB);
            }
            if (S == A || S == B)
            {
                free(S->coeffs);
                S->coeffs = s;
                S->alloc  = lenB - 1;
            }
            if (T == A || T == B)
            {
                free(T->coeffs);
                T->coeffs = t;
                T->alloc  = lenA - 1;
            }

            G->length = lenG;
            S->length = FLINT_MAX(lenB - lenG, 1);
            T->length = FLINT_MAX(lenA - lenG, 1);
            MPN_NORM(S->coeffs, S->length);
            MPN_NORM(T->coeffs, T->length);

            if (G->coeffs[lenG - 1] != 1)
            {
                inv = n_invmod(G->coeffs[lenG - 1], A->mod.n);
                nmod_poly_scalar_mul_nmod(G, G, inv);
                nmod_poly_scalar_mul_nmod(S, S, inv);
                nmod_poly_scalar_mul_nmod(T, T, inv);
            }
        }
    }
}

