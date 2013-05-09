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

#define __attach_shift(B, lenB, A, lenA, m)      \
do {                                             \
    (B) = (A) + (m);                             \
    (lenB) = ((lenA) >= (m)) ? (lenA) - (m) : 0; \
} while (0)

#define __attach_truncate(B, lenB, A, lenA, m) \
do {                                           \
    (B) = (A);                                 \
    (lenB) = ((lenA) < (m)) ? (lenA) : (m);    \
} while (0)

#define __set(B, lenB, A, lenA)      \
do {                                 \
    _nmod_vec_set((B), (A), (lenA)); \
    (lenB) = (lenA);                 \
} while (0)

#define __swap  MPN_SWAP

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

static __inline__ void __mat_one(mp_ptr *M, len_t *lenM)
{
    M[0][0] = 1L;
    M[3][0] = 1L;
    lenM[0] = 1;
    lenM[1] = 0;
    lenM[2] = 0;
    lenM[3] = 1;
}

/*
    Computes the matrix product C of the two 2x2 matrices A and B, 
    using classical multiplication.

    Does not support aliasing.

    Expects T to be temporary space sufficient for any of the 
    polynomial products involved.
 */

static void __mat_mul_classical(mp_ptr *C, len_t *lenC, 
    mp_ptr *A, len_t *lenA, mp_ptr *B, len_t *lenB, mp_ptr T, nmod_t mod)
{
    len_t lenT;

    __mul(C[0], lenC[0], A[0], lenA[0], B[0], lenB[0]);
    __mul(T, lenT, A[1], lenA[1], B[2], lenB[2]);
    __add(C[0], lenC[0], C[0], lenC[0], T, lenT);

    __mul(C[1], lenC[1], A[0], lenA[0], B[1], lenB[1]);
    __mul(T, lenT, A[1], lenA[1], B[3], lenB[3]);
    __add(C[1], lenC[1], C[1], lenC[1], T, lenT);

    __mul(C[2], lenC[2], A[2], lenA[2], B[0], lenB[0]);
    __mul(T, lenT, A[3], lenA[3], B[2], lenB[2]);
    __add(C[2], lenC[2], C[2], lenC[2], T, lenT);

    __mul(C[3], lenC[3], A[2], lenA[2], B[1], lenB[1]);
    __mul(T, lenT, A[3], lenA[3], B[3], lenB[3]);
    __add(C[3], lenC[3], C[3], lenC[3], T, lenT);
}

/*
    Computes the matrix product C of the two 2x2 matrices A and B, 
    using Strassen multiplication.

    Does not support aliasing.

    Expects T0, T1 to be temporary space sufficient for any of the 
    polynomial products involved.
 */

static void __mat_mul_strassen(mp_ptr *C, len_t *lenC, 
    mp_ptr *A, len_t *lenA, mp_ptr *B, len_t *lenB, mp_ptr T0, mp_ptr T1, 
    nmod_t mod)
{
    len_t lenT0, lenT1;

    __sub(T0, lenT0, A[0], lenA[0], A[2], lenA[2]);
    __sub(T1, lenT1, B[3], lenB[3], B[1], lenB[1]);
    __mul(C[2], lenC[2], T0, lenT0, T1, lenT1);

    __add(T0, lenT0, A[2], lenA[2], A[3], lenA[3]);
    __sub(T1, lenT1, B[1], lenB[1], B[0], lenB[0]);
    __mul(C[3], lenC[3], T0, lenT0, T1, lenT1);

    __sub(T0, lenT0, T0, lenT0, A[0], lenA[0]);
    __sub(T1, lenT1, B[3], lenB[3], T1, lenT1);
    __mul(C[1], lenC[1], T0, lenT0, T1, lenT1);

    __sub(T0, lenT0, A[1], lenA[1], T0, lenT0);
    __mul(C[0], lenC[0], T0, lenT0, B[3], lenB[3]);

    __mul(T0, lenT0, A[0], lenA[0], B[0], lenB[0]);

    __add(C[1], lenC[1], T0, lenT0, C[1], lenC[1]);
    __add(C[2], lenC[2], C[1], lenC[1], C[2], lenC[2]);
    __add(C[1], lenC[1], C[1], lenC[1], C[3], lenC[3]);
    __add(C[3], lenC[3], C[2], lenC[2], C[3], lenC[3]);
    __add(C[1], lenC[1], C[1], lenC[1], C[0], lenC[0]);
    __sub(T1, lenT1, T1, lenT1, B[2], lenB[2]);
    __mul(C[0], lenC[0], A[3], lenA[3], T1, lenT1);

    __sub(C[2], lenC[2], C[2], lenC[2], C[0], lenC[0]);
    __mul(C[0], lenC[0], A[1], lenA[1], B[2], lenB[2]);

    __add(C[0], lenC[0], C[0], lenC[0], T0, lenT0);
}

/*
    Computs the matrix product C of the two 2x2 matrices A and B, 
    using either classical or Strassen multiplication depending 
    on the degrees of the input polynomials.

    Does not support aliasing.

    Expects T0, T1 to be temporary space sufficient for any of the 
    polynomial products involved.
 */

static void __mat_mul(mp_ptr *C, len_t *lenC, 
    mp_ptr *A, len_t *lenA, mp_ptr *B, len_t *lenB, mp_ptr T0, mp_ptr T1, 
    nmod_t mod)
{
    len_t min = lenA[0];

    min = FLINT_MIN(min, lenA[1]);
    min = FLINT_MIN(min, lenA[2]);
    min = FLINT_MIN(min, lenA[3]);
    min = FLINT_MIN(min, lenB[0]);
    min = FLINT_MIN(min, lenB[1]);
    min = FLINT_MIN(min, lenB[2]);
    min = FLINT_MIN(min, lenB[3]);

    if (min < 20)
    {
        __mat_mul_classical(C, lenC, A, lenA, B, lenB, T0, mod);
    }
    else
    {
        __mat_mul_strassen(C, lenC, A, lenA, B, lenB, T0, T1, mod);
    }
}

/*
    HGCD Iterative step.

    Only supports aliasing in {*A,a} and {*B,b}.

    Assumes that lena > lenb > 0.

    Assumes that the pointers {*A, *B, *T} as well as 
    {M + 0, M + 1, M + 2, M + 3, t} may be swapped. 
    With the underlying HGCD implementation in mind, 
    this is to say that the blocks of memory implicitly 
    reserved for these pointers probably should have 
    the same size.

    Expects {*A, *B, *T} to be of size at least lena, 
    {M + 0, M + 1, M + 2, M + 3, *t} and Q of size at 
    least (lena + 1)/2.
 */

len_t _nmod_poly_hgcd_recursive_iter(mp_ptr *M, len_t *lenM, 
    mp_ptr *A, len_t *lenA, mp_ptr *B, len_t *lenB, 
    mp_srcptr a, len_t lena, mp_srcptr b, len_t lenb, 
    mp_ptr Q, mp_ptr *T, mp_ptr *t, nmod_t mod)
{
    const len_t m = lena / 2;
    len_t sgn = 1;

    __mat_one(M, lenM);
    __set(*A, *lenA, a, lena);
    __set(*B, *lenB, b, lenb);

    while (*lenB >= m + 1)
    {
        len_t lenQ, lenT, lent;

        __divrem(Q, lenQ, *T, lenT, *A, *lenA, *B, *lenB);
        __swap(*B, *lenB, *T, lenT);
        __swap(*A, *lenA, *T, lenT);

        __mul(*T, lenT, Q, lenQ, M[2], lenM[2]);
        __add(*t, lent, M[3], lenM[3], *T, lenT);
        __swap(M[3], lenM[3], M[2], lenM[2]);
        __swap(M[2], lenM[2], *t, lent);

        __mul(*T, lenT, Q, lenQ, M[0], lenM[0]);
        __add(*t, lent, M[1], lenM[1], *T, lenT);
        __swap(M[1], lenM[1], M[0], lenM[0]);
        __swap(M[0], lenM[0], *t, lent);

        sgn = -sgn;
    }

    return sgn;
}

/* 
    Assumes that lena > lenb > 0.

    The current implementation requires P to point to a memory pool 
    of size at least 6 lena + 10 (lena + 1)/2 just in this iteration.

    Supports aliasing only between {*A, a} and {*B, b}.

    Only computes the matrix {M, lenM} if flag is non-zero, in 
    which case these arrays are supposed to be sufficiently allocated. 
    Does not permute the pointers in {M, lenM}.  When flag is zero, 
    the first two arguments are allowed to be NULL.
 */

len_t _nmod_poly_hgcd_recursive(mp_ptr *M, len_t *lenM, 
    mp_ptr A, len_t *lenA, mp_ptr B, len_t *lenB, 
    mp_srcptr a, len_t lena, mp_srcptr b, len_t lenb, 
    mp_ptr P, nmod_t mod, int flag)
{
    const len_t m = lena / 2;

    if (lenb < m + 1)
    {
        if (flag)
        {
            __mat_one(M, lenM);
        }
        __set(A, *lenA, a, lena);
        __set(B, *lenB, b, lenb);
        return 1;
    }
    else
    {
        /* Readonly pointers */
        mp_ptr a0, b0, s, t, a4, b4, c0, d0;
        len_t lena0, lenb0, lens, lent, lena4, lenb4, lenc0, lend0;

        /* Pointers to independently allocated memory */
        mp_ptr a2, b2, a3, b3, q, d, T0, T1;
        len_t lena2, lenb2, lena3, lenb3, lenq, lend, lenT0;

        mp_ptr R[4], S[4];
        len_t lenR[4], lenS[4];
        len_t sgnR, sgnS;

        a2 = P;
        b2 = a2 + lena;
        a3 = b2 + lena;
        b3 = a3 + lena;
        q  = b3 + lena;
        d  = q  + (lena + 1)/2;
        T0 = d  + lena;
        T1 = T0 + lena;

        R[0] = T1   + (lena + 1)/2;
        R[1] = R[0] + (lena + 1)/2;
        R[2] = R[1] + (lena + 1)/2;
        R[3] = R[2] + (lena + 1)/2;
        S[0] = R[3] + (lena + 1)/2;
        S[1] = S[0] + (lena + 1)/2;
        S[2] = S[1] + (lena + 1)/2;
        S[3] = S[2] + (lena + 1)/2;

        P += 6 * lena + 10 * (lena + 1)/2;

        __attach_shift(a0, lena0, (mp_ptr) a, lena, m);
        __attach_shift(b0, lenb0, (mp_ptr) b, lenb, m);

        if (lena0 < NMOD_POLY_HGCD_CUTOFF)
            sgnR = _nmod_poly_hgcd_recursive_iter(R, lenR, &a3, &lena3, &b3, &lenb3, 
                                            a0, lena0, b0, lenb0, 
                                            q, &T0, &T1, mod);
        else 
            sgnR = _nmod_poly_hgcd_recursive(R, lenR, a3, &lena3, b3, &lenb3, 
                                       a0, lena0, b0, lenb0, P, mod, 1);

        __attach_truncate(s, lens, (mp_ptr) a, lena, m);
        __attach_truncate(t, lent, (mp_ptr) b, lenb, m);

        __mul(b2, lenb2, R[2], lenR[2], s, lens);
        __mul(T0, lenT0, R[0], lenR[0], t, lent);

        if (sgnR < 0)
            __sub(b2, lenb2, b2, lenb2, T0, lenT0);
        else
            __sub(b2, lenb2, T0, lenT0, b2, lenb2);

        flint_mpn_zero(b2 + lenb2, m + lenb3 - lenb2);

        __attach_shift(b4, lenb4, b2, lenb2, m);
        __add(b4, lenb4, b4, lenb4, b3, lenb3);
        lenb2 = FLINT_MAX(m + lenb3, lenb2);
        MPN_NORM(b2, lenb2);

        __mul(a2, lena2, R[3], lenR[3], s, lens);
        __mul(T0, lenT0, R[1], lenR[1], t, lent);

        if (sgnR < 0)
            __sub(a2, lena2, T0, lenT0, a2, lena2);
        else
            __sub(a2, lena2, a2, lena2, T0, lenT0);

        flint_mpn_zero(a2 + lena2, m + lena3 - lena2);
        __attach_shift(a4, lena4, a2, lena2, m);
        __add(a4, lena4, a4, lena4, a3, lena3);
        lena2 = FLINT_MAX(m + lena3, lena2);
        MPN_NORM(a2, lena2);

        if (lenb2 < m + 1)
        {
            __set(A, *lenA, a2, lena2);
            __set(B, *lenB, b2, lenb2);

            if (flag)
            {
                __set(M[0], lenM[0], R[0], lenR[0]);
                __set(M[1], lenM[1], R[1], lenR[1]);
                __set(M[2], lenM[2], R[2], lenR[2]);
                __set(M[3], lenM[3], R[3], lenR[3]);
            }

            return sgnR;
        }
        else
        {
            len_t k = 2 * m - lenb2 + 1;

            __divrem(q, lenq, d, lend, a2, lena2, b2, lenb2);

            __attach_shift(c0, lenc0, b2, lenb2, k);
            __attach_shift(d0, lend0, d, lend, k);

            if (lenc0 < NMOD_POLY_HGCD_CUTOFF)
                sgnS = _nmod_poly_hgcd_recursive_iter(S, lenS, &a3, &lena3, &b3, &lenb3, 
                                                c0, lenc0, d0, lend0, 
                                                a2, &T0, &T1, mod); /* a2 as temp */
            else 
                sgnS = _nmod_poly_hgcd_recursive(S, lenS, a3, &lena3, b3, &lenb3, 
                                           c0, lenc0, d0, lend0, P, mod, 1);

            __attach_truncate(s, lens, b2, lenb2, k);
            __attach_truncate(t, lent, d, lend, k);

            __mul(B, *lenB, S[2], lenS[2], s, lens);
            __mul(T0, lenT0, S[0], lenS[0], t, lent);

            if (sgnS < 0)
                __sub(B, *lenB, B, *lenB, T0, lenT0);
            else
                __sub(B, *lenB, T0, lenT0, B, *lenB);

            flint_mpn_zero(B + *lenB, k + lenb3 - *lenB);
            __attach_shift(b4, lenb4, B, *lenB, k);
            __add(b4, lenb4, b4, lenb4, b3, lenb3);
            *lenB = FLINT_MAX(k + lenb3, *lenB);
            MPN_NORM(B, *lenB);

            __mul(A, *lenA, S[3], lenS[3], s, lens);
            __mul(T0, lenT0, S[1], lenS[1], t, lent);

            if (sgnS < 0)
                __sub(A, *lenA, T0, lenT0, A, *lenA);
            else
                __sub(A, *lenA, A, *lenA, T0, lenT0);

            flint_mpn_zero(A + *lenA, k + lena3 - *lenA);
            __attach_shift(a4, lena4, A, *lenA, k);
            __add(a4, lena4, a4, lena4, a3, lena3);
            *lenA = FLINT_MAX(k + lena3, *lenA);
            MPN_NORM(A, *lenA);

            if (flag)
            {
                __swap(S[0], lenS[0], S[2], lenS[2]);
                __swap(S[1], lenS[1], S[3], lenS[3]);
                __mul(T0, lenT0, S[2], lenS[2], q, lenq);
                __add(S[0], lenS[0], S[0], lenS[0], T0, lenT0);
                __mul(T0, lenT0, S[3], lenS[3], q, lenq);
                __add(S[1], lenS[1], S[1], lenS[1], T0, lenT0);

                __mat_mul(M, lenM, R, lenR, S, lenS, a2, b2, mod);
            }

            return - (sgnR * sgnS);
        }
    }
}

/*
    XXX: Currently supports aliasing between {A,a} and {B,b}.
 */

len_t _nmod_poly_hgcd(mp_ptr *M, len_t *lenM, 
                     mp_ptr A, len_t *lenA, mp_ptr B, len_t *lenB, 
                     mp_srcptr a, len_t lena, mp_srcptr b, len_t lenb, 
                     nmod_t mod)
{
    const len_t lenW = 22 * lena + 16 * (FLINT_CLOG2(lena) + 1);
    len_t sgnM;
    mp_ptr W;

    W = _nmod_vec_init(lenW);

    if (M == NULL)
    {
        sgnM = _nmod_poly_hgcd_recursive(NULL, NULL, 
                                         A, lenA, B, lenB, 
                                         a, lena, b, lenb, W, mod, 0);
    }
    else
    {
        sgnM = _nmod_poly_hgcd_recursive(M, lenM, 
                                         A, lenA, B, lenB, 
                                         a, lena, b, lenb, W, mod, 1);
    }
    _nmod_vec_clear(W);

    return sgnM;
}

