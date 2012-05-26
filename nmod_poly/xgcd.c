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

long _nmod_poly_xgcd(mp_ptr G, mp_ptr S, mp_ptr T, 
                     mp_srcptr A, long lenA, mp_srcptr B, long lenB, nmod_t mod)
{
    const long cutoff = FLINT_BIT_COUNT(mod.n) <= 8 ? 
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
    const long lenA = A->length, lenB = B->length;
    mp_limb_t inv;

    if (lenA == 0)
    {
        if (lenB == 0) 
        {
            nmod_poly_zero(G);
            nmod_poly_zero(S);
            nmod_poly_zero(T);
        }
        else 
        {
            inv = n_invmod(B->coeffs[lenB - 1], B->mod.n);
            nmod_poly_scalar_mul_nmod(G, B, inv);
            nmod_poly_zero(S);
            nmod_poly_set_coeff_ui(T, 0, inv);
            T->length = 1;
        }
    } 
    else if (lenB == 0)
    {
        inv = n_invmod(A->coeffs[lenA - 1], A->mod.n);
        nmod_poly_scalar_mul_nmod(G, A, inv);
        nmod_poly_zero(T);
        nmod_poly_set_coeff_ui(S, 0, inv);
        S->length = 1;
    }
    else
    {
        nmod_poly_t tG, tS, tT;
        mp_ptr g, s, t;
        long lenG;

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
        if (S == A || S == B)
        {
            nmod_poly_init2(tS, A->mod.n, lenB);
            s = tS->coeffs;
        }
        else
        {
            nmod_poly_fit_length(S, lenB);
            s = S->coeffs;
        }
        if (T == A || T == B)
        {
            nmod_poly_init2(tT, A->mod.n, lenA);
            t = tT->coeffs;
        }
        else
        {
            nmod_poly_fit_length(T, lenA);
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
            nmod_poly_swap(tG, G);
            nmod_poly_clear(tG);
        }
        if (S == A || S == B)
        {
            nmod_poly_swap(tS, S);
            nmod_poly_clear(tS);
        }
        if (T == A || T == B)
        {
            nmod_poly_swap(tT, T);
            nmod_poly_clear(tT);
        }
        
        G->length = lenG;
        S->length = lenB - lenG;
        T->length = lenA - lenG;
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

