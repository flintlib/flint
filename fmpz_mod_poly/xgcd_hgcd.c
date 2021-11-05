/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"
#include "mpn_extras.h"

/*
    We define a whole bunch of macros here which essentially provide 
    the nmod_poly functionality as far as the setting of coefficient 
    data and lengths is concerned, but which do not do any separate 
    memory allocation.  None of these macros support aliasing.
 */

#define __set(B, lenB, A, lenA)          \
do {                                     \
    _fmpz_vec_set((B), (A), (lenA));     \
    (lenB) = (lenA);                     \
} while (0)

#define __add(C, lenC, A, lenA, B, lenB)                    \
do {                                                        \
    _fmpz_mod_poly_add((C), (A), (lenA), (B), (lenB), mod); \
    (lenC) = FLINT_MAX((lenA), (lenB));                     \
    FMPZ_VEC_NORM((C), (lenC));                             \
} while (0)

#define __sub(C, lenC, A, lenA, B, lenB)                    \
do {                                                        \
    _fmpz_mod_poly_sub((C), (A), (lenA), (B), (lenB), mod); \
    (lenC) = FLINT_MAX((lenA), (lenB));                     \
    FMPZ_VEC_NORM((C), (lenC));                             \
} while (0)

#define __mul(C, lenC, A, lenA, B, lenB)                            \
do {                                                                \
    if ((lenA) != 0 && (lenB) != 0)                                 \
    {                                                               \
        if ((lenA) >= (lenB))                                       \
            _fmpz_mod_poly_mul((C), (A), (lenA), (B), (lenB), mod); \
        else                                                        \
            _fmpz_mod_poly_mul((C), (B), (lenB), (A), (lenA), mod); \
        (lenC) = (lenA) + (lenB) - 1;                               \
    }                                                               \
    else                                                            \
    {                                                               \
        (lenC) = 0;                                                 \
    }                                                               \
} while (0)

#define __divrem(Q, lenQ, R, lenR, A, lenA, B, lenB)                          \
do {                                                                          \
    if ((lenA) >= (lenB))                                                     \
    {                                                                         \
        fmpz_invmod(invB, B + lenB - 1, mod);                                 \
        _fmpz_mod_poly_divrem((Q), (R), (A), (lenA), (B), (lenB), invB, mod); \
        (lenQ) = (lenA) - (lenB) + 1;                                         \
        (lenR) = (lenB) - 1;                                                  \
        FMPZ_VEC_NORM((R), (lenR));                                           \
    }                                                                         \
    else                                                                      \
    {                                                                         \
        _fmpz_vec_set((R), (A), (lenA));                                      \
        (lenQ) = 0;                                                           \
        (lenR) = (lenA);                                                      \
    }                                                                         \
} while (0)

#define __div(Q, lenQ, A, lenA, B, lenB)                                      \
do {                                                                          \
    if ((lenA) >= (lenB))                                                     \
    {                                                                         \
        fmpz * __t = _fmpz_vec_init(lenB - 1);                                \
        fmpz_invmod(invB, B + lenB - 1, mod);                                 \
        _fmpz_mod_poly_divrem((Q), __t, (A), (lenA), (B), (lenB), invB, mod); \
        _fmpz_vec_clear(__t, lenB - 1);                                       \
        (lenQ) = (lenA) - (lenB) + 1;                                         \
    }                                                                         \
    else                                                                      \
    {                                                                         \
        (lenQ) = 0;                                                           \
    }                                                                         \
} while (0)

slong _fmpz_mod_poly_xgcd_hgcd(fmpz *G, fmpz *S, fmpz *T, 
                          const fmpz *A, slong lenA, const fmpz *B, slong lenB, 
                          const fmpz_t mod)
{
	 slong lenG, lenS, lenT;

    if (lenB == 1)
    {
        fmpz_set(G + 0, B + 0);
        fmpz_set_ui(T, 1);
        lenG = 1;
        lenS = 0;
        lenT = 1;
    }
    else
    {
        slong lenq, lenr, len1 = lenA + lenB;

        fmpz *q = _fmpz_vec_init(len1);
        fmpz *r = q + lenA;
        fmpz_t invB;

        fmpz_init(invB);

        __divrem(q, lenq, r, lenr, A, lenA, B, lenB);
        
        if (lenr == 0)
        {
            __set(G, lenG, B, lenB);
            fmpz_set_ui(T + 0, 1);
            lenS = 0;
            lenT = 1;
        }
        else
        {
            fmpz *h, *j, *v, *w, *R[4], *X;
            slong lenh, lenj, lenv, lenw, lenR[4], len2;
            int sgnR;

            lenh = lenj = lenB;
            lenv = lenw = lenA + lenB - 2;
            lenR[0] = lenR[1] = lenR[2] = lenR[3] = (lenB + 1) / 2;

            len2 = 2 * lenh + 2 * lenv + 4 * lenR[0];
            X = _fmpz_vec_init(len2);
            h = X;
            j = h + lenh;
            v = j + lenj;
            w = v + lenv;
            R[0] = w + lenw;
            R[1] = R[0] + lenR[0];
            R[2] = R[1] + lenR[1];
            R[3] = R[2] + lenR[2];

            sgnR = _fmpz_mod_poly_hgcd(R, lenR, h, &lenh, j, &lenj, B, lenB, r, lenr, mod);

            if (sgnR > 0)
            {
                _fmpz_mod_poly_neg(S, R[1], lenR[1], mod);
                _fmpz_vec_set(T, R[0], lenR[0]);
            }
            else
            {
                _fmpz_vec_set(S, R[1], lenR[1]);
                _fmpz_mod_poly_neg(T, R[0], lenR[0], mod);
            }
            lenS = lenR[1];
            lenT = lenR[0];

            while (lenj != 0)
            {
                __divrem(q, lenq, r, lenr, h, lenh, j, lenj);
                
                __mul(v, lenv, q, lenq, T, lenT);
                {
                    slong l;
                    _fmpz_vec_swap(S, T, FLINT_MAX(lenS, lenT));
                    l = lenS; lenS = lenT; lenT = l;
                }
                if (lenr != 0) /* prevent overflow of T on last iteration */
                   __sub(T, lenT, T, lenT, v, lenv);
                else /* lenr == 0 */
                {
                    __set(G, lenG, j, lenj);

                    goto cofactor;
                }
                if (lenj < FMPZ_MOD_POLY_GCD_CUTOFF)
                {
                    fmpz *u0 = R[0], *u1 = R[1];
                    slong lenu0 = lenr - 1, lenu1 = lenj - 1;

                    fmpz_invmod(invB, r + lenr - 1, mod);
                    lenG = _fmpz_mod_poly_xgcd_euclidean(G, u0, u1, j, lenj, r, lenr, invB, mod);
                    FMPZ_VEC_NORM(u0, lenu0);
                    FMPZ_VEC_NORM(u1, lenu1);

                    __mul(v, lenv, S, lenS, u0, lenu0);
                    __mul(w, lenw, T, lenT, u1, lenu1);
                    __add(S, lenS, v, lenv, w, lenw);

                    goto cofactor;
                }

                sgnR = _fmpz_mod_poly_hgcd(R, lenR, h, &lenh, j, &lenj, j,lenj, r, lenr, mod);

                __mul(v, lenv, R[1], lenR[1], T, lenT);
                __mul(w, lenw, R[2], lenR[2], S, lenS);

                __mul(q, lenq, S, lenS, R[3], lenR[3]);
                if (sgnR > 0)
                    __sub(S, lenS, q, lenq, v, lenv);
                else
                    __sub(S, lenS, v, lenv, q, lenq);

                __mul(q, lenq, T, lenT, R[0], lenR[0]);
                if (sgnR > WORD(0))
                    __sub(T, lenT, q, lenq, w, lenw);
                else
                    __sub(T, lenT, w, lenw, q, lenq);
            }
            __set(G, lenG, h, lenh);

            cofactor:

            __mul(v, lenv, S, lenS, A, lenA);
            __sub(w, lenw, G, lenG, v, lenv);
            __div(T, lenT, w, lenw, B, lenB);

            _fmpz_vec_clear(X, len2);
        }

        _fmpz_vec_clear(q, len1);

        fmpz_clear(invB);
    }

    _fmpz_vec_zero(S + lenS, lenB - 1 - lenS);
    _fmpz_vec_zero(T + lenT, lenA - 1 - lenT);

    return lenG;
}

void fmpz_mod_poly_xgcd_hgcd(fmpz_mod_poly_t G, fmpz_mod_poly_t S,
                             fmpz_mod_poly_t T, const fmpz_mod_poly_t A,
                             const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
{
    if (A->length < B->length)
    {
        fmpz_mod_poly_xgcd_hgcd(G, T, S, B, A, ctx);
    }
    else  /* lenA >= lenB >= 0 */
    {
        const fmpz * p = fmpz_mod_ctx_modulus(ctx);
        const slong lenA = A->length, lenB = B->length;
        slong lenS, lenT;
        fmpz_t inv;

        fmpz_init(inv);

        if (lenA == 0)  /* lenA = lenB = 0 */
        {
            fmpz_mod_poly_zero(G, ctx);
            fmpz_mod_poly_zero(S, ctx);
            fmpz_mod_poly_zero(T, ctx);
        }
        else if (lenB == 0)  /* lenA > lenB = 0 */
        {
            fmpz_invmod(inv, A->coeffs + lenA - 1, p);
            fmpz_mod_poly_scalar_mul_fmpz(G, A, inv, ctx);
            fmpz_mod_poly_zero(T, ctx);
            fmpz_mod_poly_set_coeff_fmpz(S, 0, inv, ctx);
            _fmpz_mod_poly_set_length(S, 1);
        }
        else if (lenB == 1)  /* lenA >= lenB = 1 */
        {
            fmpz_mod_poly_fit_length(T, 1, ctx);
            _fmpz_mod_poly_set_length(T, 1);
            fmpz_invmod(inv, B->coeffs + 0, p);
            fmpz_set(T->coeffs + 0, inv);
            fmpz_mod_poly_set_coeff_ui(G, 0, 1, ctx);
            _fmpz_mod_poly_set_length(G, 1);
            fmpz_mod_poly_zero(S, ctx);
        }
        else  /* lenA >= lenB >= 2 */
        {
            fmpz *g, *s, *t;
            slong lenG;

            if (G == A || G == B)
            {
                g = _fmpz_vec_init(FLINT_MIN(lenA, lenB));
            }
            else
            {
                fmpz_mod_poly_fit_length(G, FLINT_MIN(lenA, lenB), ctx);
                g = G->coeffs;
            }
            if (S == A || S == B)
            {
                s = _fmpz_vec_init(FLINT_MAX(lenB - 1, 2));
            }
            else
            {
                fmpz_mod_poly_fit_length(S, FLINT_MAX(lenB - 1, 2), ctx);
                s = S->coeffs;
            }
            if (T == A || T == B)
            {
                t = _fmpz_vec_init(FLINT_MAX(lenA - 1, 2));
            }
            else
            {
                fmpz_mod_poly_fit_length(T, FLINT_MAX(lenA - 1, 2), ctx);
                t = T->coeffs;
            }

            if (lenA >= lenB)
                lenG = _fmpz_mod_poly_xgcd_hgcd(g, s, t, A->coeffs, lenA,
                                                           B->coeffs, lenB, p);
            else
                lenG = _fmpz_mod_poly_xgcd_hgcd(g, t, s, B->coeffs, lenB,
                                                           A->coeffs, lenA, p);

            if (G == A || G == B)
            {
                _fmpz_vec_clear(G->coeffs, FLINT_MIN(lenA, lenB));
                G->coeffs = g;
                G->alloc  = FLINT_MIN(lenA, lenB);
            }
            if (S == A || S == B)
            {
                _fmpz_vec_clear(S->coeffs, FLINT_MAX(lenB - 1, 2));
                S->coeffs = s;
                S->alloc  = FLINT_MAX(lenB - 1, 2);
            }
            if (T == A || T == B)
            {
                _fmpz_vec_clear(T->coeffs, FLINT_MAX(lenA - 1, 2));
                T->coeffs = t;
                T->alloc  = FLINT_MAX(lenA - 1, 2);
            }

            _fmpz_mod_poly_set_length(G, lenG);
            lenS = FLINT_MAX(lenB - lenG, 1);
            lenT = FLINT_MAX(lenA - lenG, 1);
            FMPZ_VEC_NORM(S->coeffs, lenS);
            FMPZ_VEC_NORM(T->coeffs, lenT);

            _fmpz_mod_poly_set_length(S, lenS);
            _fmpz_mod_poly_set_length(T, lenT);
            
            if (!fmpz_is_one(G->coeffs + lenG - 1))
            {
                fmpz_invmod(inv, G->coeffs + lenG - 1, p);
                fmpz_mod_poly_scalar_mul_fmpz(G, G, inv, ctx);
                fmpz_mod_poly_scalar_mul_fmpz(S, S, inv, ctx);
                fmpz_mod_poly_scalar_mul_fmpz(T, T, inv, ctx);
            }
        }

        fmpz_clear(inv);
    }
}

#undef __set
#undef __add
#undef __sub
#undef __mul
#undef __divrem
#undef __div

