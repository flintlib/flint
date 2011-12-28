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
#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "mpn_extras.h"

long _nmod_poly_xgcd_euclidean(mp_ptr res, mp_ptr s, mp_ptr t, 
              mp_srcptr poly1, long len1, mp_srcptr poly2, long len2, nmod_t mod)
{
    mp_ptr Q, R, u1, u2, v1, v2, prod, A, B;
    long u1_len, u2_len, v1_len, v2_len, Q_len, R_len, lenA, lenB, len;
    int steps;

    /* clear s, t and res */
    mpn_zero(res, len2);
    mpn_zero(s, len2-1);
    mpn_zero(t, len1-1);

    if (len2 == 1)
    {
        res[0] = poly2[0];
        t[0] = 1;
        return 1;
    }

    /* initialise arrays to store intermediate info */
    Q = _nmod_vec_init(len1);
    R = _nmod_vec_init(len2);
    u1 = _nmod_vec_init(len1);
    u2 = _nmod_vec_init(len1);
    v1 = _nmod_vec_init(len1);
    v2 = _nmod_vec_init(len1);
    prod = _nmod_vec_init(len1);
    mpn_zero(u1, len1);
    mpn_zero(u2, len1);
    mpn_zero(v1, len1);
    mpn_zero(v2, len1);
    mpn_zero(prod, len1);

    u1[0] = 1;
    u1_len = 1;
    u2_len = 0;
    v2[0] = 1;
    v2_len = 1;
    v1_len = 0;    

    steps = 0;

    A = (mp_ptr) poly1;
    lenA = len1;
    B = (mp_ptr) poly2;
    lenB = len2;

    while (lenB > 1)
    {
        _nmod_poly_divrem(Q, R, A, lenA, B, lenB, mod);
        Q_len = lenA - lenB + 1;
        R_len = lenB - 1;
        MPN_NORM(Q, Q_len);
        MPN_NORM(R, R_len);

        if (Q_len && u2_len) /* u1 = u1 - Q*u2 */
        {
            if (Q_len >= u2_len)
                _nmod_poly_mul(prod, Q, Q_len, u2, u2_len, mod);
            else
                _nmod_poly_mul(prod, u2, u2_len, Q, Q_len, mod);
            _nmod_poly_sub(u1, u1, u1_len, prod, u2_len + Q_len - 1, mod);
            u1_len = FLINT_MAX(u1_len, u2_len + Q_len - 1);
            MPN_NORM(u1, u1_len);
        }

        MPN_SWAP(u1, u1_len, u2, u2_len);

        if (Q_len && v2_len) /* v1 = v1 - Q*v2 */
        {
            if (Q_len >= v2_len)
                _nmod_poly_mul(prod, Q, Q_len, v2, v2_len, mod);
            else
                _nmod_poly_mul(prod, v2, v2_len, Q, Q_len, mod);
            _nmod_poly_sub(v1, v1, v1_len, prod, v2_len + Q_len - 1, mod);
            v1_len = FLINT_MAX(v1_len, v2_len + Q_len - 1);
            MPN_NORM(v1, v1_len);
        }

        MPN_SWAP(v1, v1_len, v2, v2_len);

        MPN_SWAP(A, lenA, B, lenB);
        MPN_SWAP(B, lenB, R, R_len);

        if (steps < 2)
        {
            R = _nmod_vec_init(lenB); /* initialise to original R_len */ 
            steps++;
        }
    }

    if (lenB == 1) 
    {      
        MPN_SWAP(u1, u1_len, u2, u2_len);
        MPN_SWAP(v1, v1_len, v2, v2_len);

        mpn_copyi(res, B, lenB);
        len = lenB;
    }
    else 
    {
        mpn_copyi(res, A, lenA);
        len = lenA;
    }

    mpn_copyi(s, u1, u1_len);
    mpn_copyi(t, v1, v1_len);

    if (steps == 2) 
        _nmod_vec_clear(A);

    _nmod_vec_clear(u1);
    _nmod_vec_clear(u2);
    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);
    _nmod_vec_clear(prod);
    _nmod_vec_clear(B);
    _nmod_vec_clear(R);
    _nmod_vec_clear(Q);

    return len;
}

void
nmod_poly_xgcd_euclidean(nmod_poly_t G, nmod_poly_t S, nmod_poly_t T,
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
            nmod_poly_init2(tS, A->mod.n, lenB - 1);
            s = tS->coeffs;
        }
        else
        {
            nmod_poly_fit_length(S, lenB);
            s = S->coeffs;
        }
        if (T == A || T == B)
        {
            nmod_poly_init2(tT, A->mod.n, lenA - 1);
            t = tT->coeffs;
        }
        else
        {
            nmod_poly_fit_length(T, lenA);
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

