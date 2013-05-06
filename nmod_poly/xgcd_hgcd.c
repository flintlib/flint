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
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"
#include "mpn_extras.h"

/*
    We define a whole bunch of macros here which essentially provide 
    the nmod_poly functionality as far as the setting of coefficient 
    data and lengths is concerned, but which do not do any separate 
    memory allocation.  None of these macros support aliasing.
 */

#define __set(B, lenB, A, lenA)      \
do {                                 \
    _nmod_vec_set((B), (A), (lenA)); \
    (lenB) = (lenA);                 \
} while (0)

#define __add(C, lenC, A, lenA, B, lenB)                \
do {                                                    \
    _nmod_poly_add((C), (A), (lenA), (B), (lenB), mod); \
    (lenC) = FLINT_MAX((lenA), (lenB));                 \
    MPN_NORM((C), (lenC));                              \
} while (0)

#define __sub(C, lenC, A, lenA, B, lenB)                \
do {                                                    \
    _nmod_poly_sub((C), (A), (lenA), (B), (lenB), mod); \
    (lenC) = FLINT_MAX((lenA), (lenB));                 \
    MPN_NORM((C), (lenC));                              \
} while (0)

#define __mul(C, lenC, A, lenA, B, lenB)                        \
do {                                                            \
    if ((lenA) != 0 && (lenB) != 0)                             \
    {                                                           \
        if ((lenA) >= (lenB))                                   \
            _nmod_poly_mul((C), (A), (lenA), (B), (lenB), mod); \
        else                                                    \
            _nmod_poly_mul((C), (B), (lenB), (A), (lenA), mod); \
        (lenC) = (lenA) + (lenB) - 1;                           \
    }                                                           \
    else                                                        \
    {                                                           \
        (lenC) = 0;                                             \
    }                                                           \
} while (0)

#define __divrem(Q, lenQ, R, lenR, A, lenA, B, lenB)                \
do {                                                                \
    if ((lenA) >= (lenB))                                           \
    {                                                               \
        _nmod_poly_divrem((Q), (R), (A), (lenA), (B), (lenB), mod); \
        (lenQ) = (lenA) - (lenB) + 1;                               \
        (lenR) = (lenB) - 1;                                        \
        MPN_NORM((R), (lenR));                                      \
    }                                                               \
    else                                                            \
    {                                                               \
        _nmod_vec_set((R), (A), (lenA));                            \
        (lenQ) = 0;                                                 \
        (lenR) = (lenA);                                            \
    }                                                               \
} while (0)

#define __div(Q, lenQ, A, lenA, B, lenB)                    \
do {                                                        \
    if ((lenA) >= (lenB))                                   \
    {                                                       \
        _nmod_poly_div((Q), (A), (lenA), (B), (lenB), mod); \
        (lenQ) = (lenA) - (lenB) + 1;                       \
    }                                                       \
    else                                                    \
    {                                                       \
        (lenQ) = 0;                                         \
    }                                                       \
} while (0)

long _nmod_poly_xgcd_hgcd(mp_ptr G, mp_ptr S, mp_ptr T, 
                          mp_srcptr A, long lenA, mp_srcptr B, long lenB, 
                          nmod_t mod)
{
	const long cutoff = FLINT_BIT_COUNT(mod.n) <= 8 ? 
                        NMOD_POLY_SMALL_GCD_CUTOFF : NMOD_POLY_GCD_CUTOFF;

    long lenG, lenS, lenT;

    if (lenB == 1)
    {
        G[0] = B[0];
        T[0] = 1;
        lenG = 1;
        lenS = 0;
        lenT = 1;
    }
    else
    {
        mp_ptr q = _nmod_vec_init(lenA + lenB);
        mp_ptr r = q + lenA;

        long lenq, lenr;

        __divrem(q, lenq, r, lenr, A, lenA, B, lenB);

        if (lenr == 0)
        {
            __set(G, lenG, B, lenB);
            T[0] = 1;
            lenS = 0;
            lenT = 1;
        }
        else
        {
            mp_ptr h, j, v, w, R[4], X;
            long lenh, lenj, lenv, lenw, lenR[4];
            int sgnR;

            lenh = lenj = lenB;
            lenv = lenw = lenA + lenB - 2;
            lenR[0] = lenR[1] = lenR[2] = lenR[3] = (lenB + 1) / 2;

            X = _nmod_vec_init(2 * lenh + 2 * lenv + 4 * lenR[0]);
            h = X;
            j = h + lenh;
            v = j + lenj;
            w = v + lenv;
            R[0] = w + lenw;
            R[1] = R[0] + lenR[0];
            R[2] = R[1] + lenR[1];
            R[3] = R[2] + lenR[2];

            sgnR = _nmod_poly_hgcd(R, lenR, h, &lenh, j, &lenj, B, lenB, r, lenr, mod);

            if (sgnR > 0)
            {
                _nmod_vec_neg(S, R[1], lenR[1], mod);
                _nmod_vec_set(T, R[0], lenR[0]);
            }
            else
            {
                _nmod_vec_set(S, R[1], lenR[1]);
                _nmod_vec_neg(T, R[0], lenR[0], mod);
            }
            lenS = lenR[1];
            lenT = lenR[0];

            while (lenj != 0)
            {
                __divrem(q, lenq, r, lenr, h, lenh, j, lenj);
                __mul(v, lenv, q, lenq, T, lenT);
                {
                    long l;
                    _nmod_vec_swap(S, T, FLINT_MAX(lenS, lenT));
                    l = lenS;lenS = lenT; lenT = l;
                }
                __sub(T, lenT, T, lenT, v, lenv);

                if (lenr == 0)
                {
                    __set(G, lenG, j, lenj);

                    goto cofactor;
                }
                if (lenj < cutoff)
                {
                    mp_ptr u0 = R[0], u1 = R[1];
                    long lenu0 = lenr - 1, lenu1 = lenj - 1;

                    lenG = _nmod_poly_xgcd_euclidean(G, u0, u1, j, lenj, r, lenr, mod);
                    MPN_NORM(u0, lenu0);
                    MPN_NORM(u1, lenu1);

                    __mul(v, lenv, S, lenS, u0, lenu0);
                    __mul(w, lenw, T, lenT, u1, lenu1);
                    __add(S, lenS, v, lenv, w, lenw);

                    goto cofactor;
                }

                sgnR = _nmod_poly_hgcd(R, lenR, h, &lenh, j, &lenj, j,lenj, r, lenr, mod);

                __mul(v, lenv, R[1], lenR[1], T, lenT);
                __mul(w, lenw, R[2], lenR[2], S, lenS);

                __mul(q, lenq, S, lenS, R[3], lenR[3]);
                if (sgnR > 0)
                    __sub(S, lenS, q, lenq, v, lenv);
                else
                    __sub(S, lenS, v, lenv, q, lenq);

                __mul(q, lenq, T, lenT, R[0], lenR[0]);
                if (sgnR > 0L)
                    __sub(T, lenT, q, lenq, w, lenw);
                else
                    __sub(T, lenT, w, lenw, q, lenq);
            }
            __set(G, lenG, h, lenh);

            cofactor:

            __mul(v, lenv, S, lenS, A, lenA);
            __sub(w, lenw, G, lenG, v, lenv);
            __div(T, lenT, w, lenw, B, lenB);

            _nmod_vec_clear(X);
        }
        _nmod_vec_clear(q);
    }
    flint_mpn_zero(S + lenS, lenB - 1 - lenS);
    flint_mpn_zero(T + lenT, lenA - 1 - lenT);

    return lenG;
}

void
nmod_poly_xgcd_hgcd(nmod_poly_t G, nmod_poly_t S, nmod_poly_t T,
                         const nmod_poly_t A, const nmod_poly_t B)
{
    if (A->length < B->length)
    {
        nmod_poly_xgcd_hgcd(G, T, S, B, A);
    }
    else  /* lenA >= lenB >= 0 */
    {
        const long lenA = A->length, lenB = B->length;
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
            long lenG;

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
                lenG = _nmod_poly_xgcd_hgcd(g, s, t, A->coeffs, lenA,
                                                          B->coeffs, lenB, A->mod);
            else
                lenG = _nmod_poly_xgcd_hgcd(g, t, s, B->coeffs, lenB,
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

#undef __set
#undef __add
#undef __sub
#undef __mul
#undef __divrem
#undef __div

