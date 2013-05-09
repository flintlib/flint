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

len_t _nmod_poly_xgcd_euclidean(mp_ptr G, mp_ptr S, mp_ptr T, 
                               mp_srcptr A, len_t lenA, 
                               mp_srcptr B, len_t lenB, nmod_t mod)
{
    flint_mpn_zero(G, lenB);
    flint_mpn_zero(S, lenB - 1);
    flint_mpn_zero(T, lenA - 1);

    if (lenB == 1)
    {
        G[0] = B[0];
        T[0] = 1;
        return 1;
    }
    else
    {
        mp_ptr Q, R;
        len_t lenQ, lenR, lenG;

        Q = _nmod_vec_init(2 * lenA);
        R = Q + lenA;

        _nmod_poly_divrem(Q, R, A, lenA, B, lenB, mod);
        lenR = lenB - 1;
        MPN_NORM(R, lenR);

        if (lenR == 0)
        {
            _nmod_vec_set(G, B, lenB);
            T[0] = 1;
            lenG = lenB;
        }
        else
        {
            mp_ptr D, U, V1, V3, W;
            len_t lenD, lenU, lenV1, lenV3, lenW;

            W  = _nmod_vec_init(FLINT_MAX(5 * lenB, lenA + lenB));
            D  = W  + lenB;
            U  = D  + lenB;
            V1 = U  + lenB;
            V3 = V1 + lenB;

            lenU = 0;
            _nmod_vec_set(D, B, lenB);
            lenD = lenB;
            V1[0] = 1;
            lenV1 = 1;
            lenV3 = 0;
            MPN_SWAP(V3, lenV3, R, lenR);

            do {
                _nmod_poly_divrem(Q, R, D, lenD, V3, lenV3, mod);
                lenQ = lenD - lenV3 + 1;
                lenR = lenV3 - 1;
                MPN_NORM(R, lenR);

                if (lenV1 >= lenQ)
                    _nmod_poly_mul(W, V1, lenV1, Q, lenQ, mod);
                else
                    _nmod_poly_mul(W, Q, lenQ, V1, lenV1, mod);
                lenW = lenQ + lenV1 - 1;

                _nmod_poly_sub(U, U, lenU, W, lenW, mod);
                lenU = FLINT_MAX(lenU, lenW);
                MPN_NORM(U, lenU);

                MPN_SWAP(U, lenU, V1, lenV1);
                {
                    mp_ptr __t;
                    len_t __tn;

                    __t = D;
                    D   = V3;
                    V3  = R;
                    R   = __t;
                    __tn  = lenD;
                    lenD  = lenV3;
                    lenV3 = lenR;
                    lenR  = __tn;
                }

            } while (lenV3 != 0);

            _nmod_vec_set(G, D, lenD);
            _nmod_vec_set(S, U, lenU);

            {
                lenQ = lenA + lenU - 1;

                _nmod_poly_mul(Q, A, lenA, S, lenU, mod);
                _nmod_vec_neg(Q, Q, lenQ, mod);
                _nmod_poly_add(Q, G, lenD, Q, lenQ, mod);

                _nmod_poly_divrem(T, W, Q, lenQ, B, lenB, mod);
            }

            _nmod_vec_clear(W);
            lenG = lenD;
        }

        _nmod_vec_clear(Q);
        return lenG;
    }
}

void
nmod_poly_xgcd_euclidean(nmod_poly_t G, nmod_poly_t S, nmod_poly_t T,
                         const nmod_poly_t A, const nmod_poly_t B)
{
    if (A->length < B->length)
    {
        nmod_poly_xgcd_euclidean(G, T, S, B, A);
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
                lenG = _nmod_poly_xgcd_euclidean(g, s, t, A->coeffs, lenA,
                                                          B->coeffs, lenB, A->mod);
            else
                lenG = _nmod_poly_xgcd_euclidean(g, t, s, B->coeffs, lenB,
                                                          A->coeffs, lenA, A->mod);

            if (G == A || G == B)
            {
                flint_free(G->coeffs);
                G->coeffs = g;
                G->alloc  = FLINT_MIN(lenA, lenB);
            }
            if (S == A || S == B)
            {
                flint_free(S->coeffs);
                S->coeffs = s;
                S->alloc  = lenB - 1;
            }
            if (T == A || T == B)
            {
                flint_free(T->coeffs);
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

