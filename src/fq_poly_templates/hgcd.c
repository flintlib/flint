/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

/*
    We define a whole bunch of macros here which essentially provide 
    the TEMPLATE(T, poly) functionality as far as the setting of coefficient 
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

#define __set(B, lenB, A, lenA)                     \
do {                                                \
    _TEMPLATE(T, vec_set)((B), (A), (lenA), ctx);   \
    (lenB) = (lenA);                                \
} while (0)

#define __swap(u, l, v, m) \
  do {								\
    { TEMPLATE(T, struct)* _; _ = (u), (u) = (v), (v) = _;}     \
    { slong _; _ = (l), (l) = (m), (m) = _;}			\
  } while (0)

#define __add(C, lenC, A, lenA, B, lenB)                \
do {                                                    \
    _TEMPLATE(T, poly_add)((C), (A), (lenA), (B), (lenB), ctx); \
    (lenC) = FLINT_MAX((lenA), (lenB));                 \
    TEMPLATE(CAP_T, VEC_NORM)((C), (lenC), ctx);        \
} while (0)

#define __sub(C, lenC, A, lenA, B, lenB)               \
do {                                                   \
    _TEMPLATE(T, poly_sub)((C), (A), (lenA), (B), (lenB), ctx); \
    (lenC) = FLINT_MAX((lenA), (lenB));                 \
    TEMPLATE(CAP_T, VEC_NORM)((C), (lenC), ctx);        \
} while (0)

#define __mul(C, lenC, A, lenA, B, lenB)                        \
do {                                                            \
    if ((lenA) != 0 && (lenB) != 0)                             \
    {                                                           \
        if ((lenA) >= (lenB))                                   \
            _TEMPLATE(T, poly_mul)((C), (A), (lenA), (B), (lenB), ctx); \
        else                                                    \
            _TEMPLATE(T, poly_mul)((C), (B), (lenB), (A), (lenA), ctx); \
        (lenC) = (lenA) + (lenB) - 1;                           \
    }                                                           \
    else                                                        \
    {                                                           \
        (lenC) = 0;                                             \
    }                                                           \
} while (0)

#define __divrem(Q, lenQ, R, lenR, A, lenA, B, lenB, invB)               \
do {                                                                \
    if ((lenA) >= (lenB))                                           \
    {                                                               \
        _TEMPLATE(T, poly_divrem)((Q), (R), (A), (lenA), (B), (lenB), (invB), ctx); \
        (lenQ) = (lenA) - (lenB) + 1;                               \
        (lenR) = (lenB) - 1;                                        \
        TEMPLATE(CAP_T, VEC_NORM)((R), (lenR), ctx);                    \
    }                                                               \
    else                                                            \
    {                                                               \
        _TEMPLATE(T, vec_set)((R), (A), (lenA), ctx);                   \
        (lenQ) = 0;                                                 \
        (lenR) = (lenA);                                            \
    }                                                               \
} while (0)

static __inline__ void
__mat_one(TEMPLATE(T, struct) ** M, slong * lenM, const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, one) (M[0], ctx);
    TEMPLATE(T, one) (M[3], ctx);
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

static void
__mat_mul_classical(TEMPLATE(T, struct) ** C, slong * lenC,
                    TEMPLATE(T, struct) ** A, slong * lenA,
                    TEMPLATE(T, struct) ** B, slong * lenB,
                    TEMPLATE(T, struct) * T, const TEMPLATE(T, ctx_t) ctx)
{
    slong lenT;

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

static void
__mat_mul_strassen(TEMPLATE(T, struct) ** C, slong * lenC,
                   TEMPLATE(T, struct) ** A, slong * lenA,
                   TEMPLATE(T, struct) ** B, slong * lenB,
                   TEMPLATE(T, struct) * T0,
                   TEMPLATE(T, struct) * T1, const TEMPLATE(T, ctx_t) ctx)
{
    slong lenT0, lenT1;

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

static void
__mat_mul(TEMPLATE(T, struct) ** C, slong * lenC,
          TEMPLATE(T, struct) ** A, slong * lenA,
          TEMPLATE(T, struct) ** B, slong * lenB,
          TEMPLATE(T, struct) * T0, TEMPLATE(T, struct) * T1,
          const TEMPLATE(T, ctx_t) ctx)
{
    slong min = lenA[0];

    min = FLINT_MIN(min, lenA[1]);
    min = FLINT_MIN(min, lenA[2]);
    min = FLINT_MIN(min, lenA[3]);
    min = FLINT_MIN(min, lenB[0]);
    min = FLINT_MIN(min, lenB[1]);
    min = FLINT_MIN(min, lenB[2]);
    min = FLINT_MIN(min, lenB[3]);

    if (min < 20)
    {
        __mat_mul_classical(C, lenC, A, lenA, B, lenB, T0, ctx);
    }
    else
    {
        __mat_mul_strassen(C, lenC, A, lenA, B, lenB, T0, T1, ctx);
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


slong
_TEMPLATE(T, poly_hgcd_recursive_iter) (
    TEMPLATE(T, struct) ** M, slong * lenM,
    TEMPLATE(T, struct) ** A, slong * lenA,
    TEMPLATE(T, struct) ** B, slong * lenB,
    const TEMPLATE(T, struct) * a, slong lena,
    const TEMPLATE(T, struct) * b, slong lenb,
    TEMPLATE(T, struct) * Q, TEMPLATE(T, struct) ** T,
    TEMPLATE(T, struct) **t,
    const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, t) invB;
    const slong m = lena / 2;
    slong sgn = 1;

    TEMPLATE(CAP_T, VEC_NORM) (b, lenb, ctx);

    __mat_one(M, lenM, ctx);
    __set(*A, *lenA, a, lena);
    __set(*B, *lenB, b, lenb);

    TEMPLATE(T, init) (invB, ctx);

    while (*lenB >= m + 1)
    {
        slong lenQ, lenT, lent;

        TEMPLATE(T, inv) (invB, *B + *lenB - 1, ctx);
        __divrem(Q, lenQ, *T, lenT, *A, *lenA, *B, *lenB, invB);
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

    TEMPLATE(T, clear) (invB, ctx);

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

slong _TEMPLATE(T, poly_hgcd_recursive) (
    TEMPLATE(T, struct)**M, slong * lenM,
    TEMPLATE(T, struct) * A, slong * lenA,
    TEMPLATE(T, struct) * B, slong * lenB,
    const TEMPLATE(T, struct) * a, slong lena,
    const TEMPLATE(T, struct) * b, slong lenb,
    TEMPLATE(T, struct) * P,
    const TEMPLATE(T, ctx_t) ctx, int flag)
{
    const slong m = lena / 2;

    if (lenb < m + 1)
    {
        if (flag)
        {
            __mat_one(M, lenM, ctx);
        }
        __set(A, *lenA, a, lena);
        __set(B, *lenB, b, lenb);
        return 1;
    }
    else
    {
        /* Readonly pointers */
        TEMPLATE(T, struct) * a0, *b0, *s, *t, *a4, *b4, *c0, *d0;
        slong lena0, lenb0, lens, lent, lena4, lenb4, lenc0, lend0;

        /* Pointers to independently allocated memory */
        TEMPLATE(T, struct) * a2, *b2, *a3, *b3, *q, *d, *T0, *T1;
        slong lena2, lenb2, lena3, lenb3, lenq, lend, lenT0;

        TEMPLATE(T, struct) * R[4], *S[4];
        slong lenR[4], lenS[4];
        slong sgnR, sgnS;

        a2 = P;
        b2 = a2 + lena;
        a3 = b2 + lena;
        b3 = a3 + lena;
        q = b3 + lena;
        d = q + (lena + 1) / 2;
        T0 = d + lena;
        T1 = T0 + lena;

        R[0] = T1 + (lena + 1) / 2;
        R[1] = R[0] + (lena + 1) / 2;
        R[2] = R[1] + (lena + 1) / 2;
        R[3] = R[2] + (lena + 1) / 2;
        S[0] = R[3] + (lena + 1) / 2;
        S[1] = S[0] + (lena + 1) / 2;
        S[2] = S[1] + (lena + 1) / 2;
        S[3] = S[2] + (lena + 1) / 2;

        P += 6 * lena + 10 * (lena + 1) / 2;

        __attach_shift(a0, lena0, (TEMPLATE(T, struct) *) a, lena, m);
        __attach_shift(b0, lenb0, (TEMPLATE(T, struct) *) b, lenb, m);

        if (lena0 < TEMPLATE(CAP_T, POLY_HGCD_CUTOFF))
            sgnR =
                _TEMPLATE(T, poly_hgcd_recursive_iter) (R, lenR, &a3, &lena3,
                                                        &b3, &lenb3, a0, lena0,
                                                        b0, lenb0, q, &T0, &T1,
                                                        ctx);
        else
            sgnR =
                _TEMPLATE(T, poly_hgcd_recursive) (R, lenR, a3, &lena3, b3,
                                                   &lenb3, a0, lena0, b0,
                                                   lenb0, P, ctx, 1);

        __attach_truncate(s, lens, (TEMPLATE(T, struct) *) a, lena, m);
        __attach_truncate(t, lent, (TEMPLATE(T, struct) *) b, lenb, m);

        __mul(b2, lenb2, R[2], lenR[2], s, lens);
        __mul(T0, lenT0, R[0], lenR[0], t, lent);

        if (sgnR < 0)
            __sub(b2, lenb2, b2, lenb2, T0, lenT0);
        else
            __sub(b2, lenb2, T0, lenT0, b2, lenb2);

        _TEMPLATE(T, vec_zero) (b2 + lenb2, m + lenb3 - lenb2, ctx);

        __attach_shift(b4, lenb4, b2, lenb2, m);
        __add(b4, lenb4, b4, lenb4, b3, lenb3);
        lenb2 = FLINT_MAX(m + lenb3, lenb2);
        TEMPLATE(CAP_T, VEC_NORM) (b2, lenb2, ctx);

        __mul(a2, lena2, R[3], lenR[3], s, lens);
        __mul(T0, lenT0, R[1], lenR[1], t, lent);

        if (sgnR < 0)
            __sub(a2, lena2, T0, lenT0, a2, lena2);
        else
            __sub(a2, lena2, a2, lena2, T0, lenT0);

        _TEMPLATE(T, vec_zero) (a2 + lena2, m + lena3 - lena2, ctx);
        __attach_shift(a4, lena4, a2, lena2, m);
        __add(a4, lena4, a4, lena4, a3, lena3);
        lena2 = FLINT_MAX(m + lena3, lena2);
        TEMPLATE(CAP_T, VEC_NORM) (a2, lena2, ctx);

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
            TEMPLATE(T, t) invB;
            slong k = 2 * m - lenb2 + 1;

            TEMPLATE(T, init) (invB, ctx);
            TEMPLATE(T, inv) (invB, b2 + lenb2 - 1, ctx);
            __divrem(q, lenq, d, lend, a2, lena2, b2, lenb2, invB);
            TEMPLATE(T, clear) (invB, ctx);
            __attach_shift(c0, lenc0, b2, lenb2, k);
            __attach_shift(d0, lend0, d, lend, k);

            if (lenc0 < TEMPLATE(CAP_T, POLY_HGCD_CUTOFF))
                sgnS = _TEMPLATE(T, poly_hgcd_recursive_iter) (
                    S, lenS, &a3, &lena3, &b3, &lenb3, c0, lenc0, d0, lend0, a2,
                    &T0, &T1, ctx); /* a2 as temp */
            else
                sgnS =
                    _TEMPLATE(T, poly_hgcd_recursive) (S, lenS, a3, &lena3, b3,
                                                       &lenb3, c0, lenc0, d0,
                                                       lend0, P, ctx, 1);

            __attach_truncate(s, lens, b2, lenb2, k);
            __attach_truncate(t, lent, d, lend, k);

            __mul(B, *lenB, S[2], lenS[2], s, lens);
            __mul(T0, lenT0, S[0], lenS[0], t, lent);

            if (sgnS < 0)
                __sub(B, *lenB, B, *lenB, T0, lenT0);
            else
                __sub(B, *lenB, T0, lenT0, B, *lenB);

            _TEMPLATE(T, vec_zero) (B + *lenB, k + lenb3 - *lenB, ctx);
            __attach_shift(b4, lenb4, B, *lenB, k);
            __add(b4, lenb4, b4, lenb4, b3, lenb3);
            *lenB = FLINT_MAX(k + lenb3, *lenB);
            TEMPLATE(CAP_T, VEC_NORM) (B, *lenB, ctx);

            __mul(A, *lenA, S[3], lenS[3], s, lens);
            __mul(T0, lenT0, S[1], lenS[1], t, lent);

            if (sgnS < 0)
                __sub(A, *lenA, T0, lenT0, A, *lenA);
            else
                __sub(A, *lenA, A, *lenA, T0, lenT0);

            _TEMPLATE(T, vec_zero) (A + *lenA, k + lena3 - *lenA, ctx);
            __attach_shift(a4, lena4, A, *lenA, k);
            __add(a4, lena4, a4, lena4, a3, lena3);
            *lenA = FLINT_MAX(k + lena3, *lenA);
            TEMPLATE(CAP_T, VEC_NORM) (A, *lenA, ctx);

            if (flag)
            {
                __swap(S[0], lenS[0], S[2], lenS[2]);
                __swap(S[1], lenS[1], S[3], lenS[3]);
                __mul(T0, lenT0, S[2], lenS[2], q, lenq);
                __add(S[0], lenS[0], S[0], lenS[0], T0, lenT0);
                __mul(T0, lenT0, S[3], lenS[3], q, lenq);
                __add(S[1], lenS[1], S[1], lenS[1], T0, lenT0);

                __mat_mul(M, lenM, R, lenR, S, lenS, a2, b2, ctx);
            }

            return -(sgnR * sgnS);
        }
    }
}

/*
    XXX: Currently supports aliasing between {A,a} and {B,b}.
 */

slong _TEMPLATE(T, poly_hgcd) (TEMPLATE(T, struct)**M, slong * lenM,
                               TEMPLATE(T, struct) * A, slong * lenA,
                               TEMPLATE(T, struct) * B, slong * lenB,
                               const TEMPLATE(T, struct) * a, slong lena,
                               const TEMPLATE(T, struct) * b, slong lenb,
                               const TEMPLATE(T, ctx_t) ctx)
{
    const slong lenW = 22 * lena + 16 * (FLINT_CLOG2(lena) + 1);
    slong sgnM;
    TEMPLATE(T, struct) * W;

    W = _TEMPLATE(T, vec_init) (lenW, ctx);

    if (M == NULL)
    {
        sgnM = _TEMPLATE(T, poly_hgcd_recursive) (NULL, NULL,
                                                  A, lenA, B, lenB,
                                                  a, lena, b, lenb, W, ctx, 0);
    }
    else
    {
        sgnM = _TEMPLATE(T, poly_hgcd_recursive) (M, lenM,
                                                  A, lenA, B, lenB,
                                                  a, lena, b, lenb, W, ctx, 1);
    }
    _TEMPLATE(T, vec_clear) (W, lenW, ctx);

    return sgnM;
}



#endif
