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

#include "gr_poly.h"
#include "gr_vec.h"

typedef struct
{
   gr_ptr res;
   gr_ptr lc;
   slong len0;
   slong len1;
   slong off;
} gr_poly_res_struct;

typedef gr_poly_res_struct gr_poly_res_t[1];


#define OFFSET(x, i) GR_ENTRY((x), (i), sz)

#define GR_VEC_NORM(status, R, lenR, sz, ctx) \
    do { \
        (void) sz; \
        (status) |= _gr_vec_normalise(&(lenR), (R), (lenR), (ctx)); \
    } while (0)
/*
    We define a whole bunch of macros here which essentially provide
    the gr_poly functionality as far as the setting of coefficient
    data and lengths is concerned, but which do not do any separate
    memory allocation.  None of these macros support aliasing.
 */

#define __attach_shift(B, lenB, A, lenA, m)      \
do {                                             \
    (B) = OFFSET(A, m);                          \
    (lenB) = ((lenA) >= (m)) ? (lenA) - (m) : 0; \
} while (0)

#define __attach_truncate(B, lenB, A, lenA, m) \
do {                                           \
    (B) = (A);                                 \
    (lenB) = ((lenA) < (m)) ? (lenA) : (m);    \
} while (0)

#define __set(B, lenB, A, lenA)                     \
do {                                                \
    status |= _gr_vec_set((B), (A), (lenA), ctx);   \
    (lenB) = (lenA);                                \
} while (0)

#define __swap(u, l, v, m) \
  do {                                                  \
    { gr_ptr _; _ = (u), (u) = (v), (v) = _;}           \
    { slong _; _ = (l), (l) = (m), (m) = _;}			\
  } while (0)

#define __add(C, lenC, A, lenA, B, lenB)                         \
do {                                                             \
    status |= _gr_poly_add((C), (A), (lenA), (B), (lenB), ctx);  \
    (lenC) = FLINT_MAX((lenA), (lenB));                          \
    GR_VEC_NORM(status, C, lenC, sz, ctx);                       \
} while (0)

#define __sub(C, lenC, A, lenA, B, lenB)                         \
do {                                                             \
    status |= _gr_poly_sub((C), (A), (lenA), (B), (lenB), ctx);  \
    (lenC) = FLINT_MAX((lenA), (lenB));                          \
    GR_VEC_NORM(status, C, lenC, sz, ctx);                       \
} while (0)

#define __mul(C, lenC, A, lenA, B, lenB)                        \
do {                                                            \
    if ((lenA) != 0 && (lenB) != 0)                             \
    {                                                           \
        if ((lenA) >= (lenB))                                   \
            status |= _gr_poly_mul((C), (A), (lenA), (B), (lenB), ctx); \
        else                                                    \
            status |= _gr_poly_mul((C), (B), (lenB), (A), (lenA), ctx); \
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
        status |= _gr_poly_divrem((Q), (R), (A), (lenA), (B), (lenB), ctx); \
        (lenQ) = (lenA) - (lenB) + 1;                               \
        (lenR) = (lenB) - 1;                                        \
        GR_VEC_NORM(status, R, lenR, sz, ctx);                      \
    }                                                               \
    else                                                            \
    {                                                               \
        status |= _gr_vec_set((R), (A), (lenA), ctx);               \
        (lenQ) = 0;                                                 \
        (lenR) = (lenA);                                            \
    }                                                               \
} while (0)

static inline int
__mat_one(gr_ptr * M, slong * lenM, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    status |= gr_one(M[0], ctx);
    status |= gr_one(M[3], ctx);
    lenM[0] = 1;
    lenM[1] = 0;
    lenM[2] = 0;
    lenM[3] = 1;
    return status;
}

/*
    Computes the matrix product C of the two 2x2 matrices A and B,
    using classical multiplication.

    Does not support aliasing.

    Expects T to be temporary space sufficient for any of the
    polynomial products involved.
 */

static int
__mat_mul_classical(gr_ptr * C, slong * lenC,
                    gr_ptr * A, slong * lenA,
                    gr_ptr * B, slong * lenB,
                    gr_ptr T, gr_ctx_t ctx)
{
    slong lenT;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

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

    return status;
}

/*
    Computes the matrix product C of the two 2x2 matrices A and B,
    using Strassen multiplication.

    Does not support aliasing.

    Expects T0, T1 to be temporary space sufficient for any of the
    polynomial products involved.
 */

static int
__mat_mul_strassen(gr_ptr * C, slong * lenC,
                   gr_ptr * A, slong * lenA,
                   gr_ptr * B, slong * lenB,
                   gr_ptr T0,
                   gr_ptr T1, gr_ctx_t ctx)
{
    slong lenT0, lenT1;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

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

    return status;
}

/*
    Computs the matrix product C of the two 2x2 matrices A and B,
    using either classical or Strassen multiplication depending
    on the degrees of the input polynomials.

    Does not support aliasing.

    Expects T0, T1 to be temporary space sufficient for any of the
    polynomial products involved.
 */

static int
__mat_mul(gr_ptr * C, slong * lenC,
          gr_ptr * A, slong * lenA,
          gr_ptr * B, slong * lenB,
          gr_ptr T0, gr_ptr T1,
          gr_ctx_t ctx)
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
        return __mat_mul_classical(C, lenC, A, lenA, B, lenB, T0, ctx);
    }
    else
    {
        return __mat_mul_strassen(C, lenC, A, lenA, B, lenB, T0, T1, ctx);
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


int
_gr_poly_hgcd_recursive_iter(
    slong * ans_sgn,
    gr_ptr * M, slong * lenM,
    gr_ptr * A, slong * lenA,
    gr_ptr * B, slong * lenB,
    gr_srcptr a, slong lena,
    gr_srcptr b, slong lenb,
    gr_ptr Q, gr_ptr * T,
    gr_ptr * t,
    gr_ctx_t ctx,
    gr_poly_res_t res)
{
    const slong m = lena / 2;
    slong sgn = 1;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    GR_VEC_NORM(status, b, lenb, sz, ctx);

    __mat_one(M, lenM, ctx);
    __set(*A, *lenA, a, lena);
    __set(*B, *lenB, b, lenb);

    while (*lenB >= m + 1)
    {
        slong lenQ, lenT, lent;

        if (res != NULL)
           status |= gr_set(res->lc, OFFSET(*B, *lenB - 1), ctx);

        __divrem(Q, lenQ, *T, lenT, *A, *lenA, *B, *lenB);

        if (res != NULL)
        {
            if (lenT >= m + 1)
            {
                if (lenT >= 1)
                {
                    status |= gr_pow_ui(res->lc, res->lc, *lenA - lenT, ctx);
                    status |= gr_mul(res->res, res->res, res->lc, ctx);

                    if ((((*lenA + res->off) | (*lenB + res->off)) & 1) == 0)
                        status |= gr_neg(res->res, res->res, ctx);
                }
                else
                {
                    if (*lenB == 1)
                    {
                        status |= gr_pow_ui(res->lc, res->lc, *lenA - 1, ctx);
                        status |= gr_mul(res->res, res->res, res->lc, ctx);
                    }
                    else
                        status |= gr_zero(res->res, ctx);
                }
            }
            else
            {
                res->len0 = *lenA;
                res->len1 = *lenB;
            }
        }

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

    *ans_sgn = sgn;

    return status;
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

int _gr_poly_hgcd_recursive(
    slong * ans_sgn,
    gr_ptr * M, slong * lenM,
    gr_ptr A, slong * lenA,
    gr_ptr B, slong * lenB,
    gr_srcptr a, slong lena,
    gr_srcptr b, slong lenb,
    gr_ptr P,
    gr_ctx_t ctx, slong cutoff, int flag, gr_poly_res_t res)
{
    const slong m = lena / 2;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (lenb < m + 1)
    {
        if (flag)
        {
            __mat_one(M, lenM, ctx);
        }
        __set(A, *lenA, a, lena);
        __set(B, *lenB, b, lenb);

        *ans_sgn = 1;
        return status;
    }
    else
    {
        /* Readonly pointers */
        gr_ptr a0, b0, s, t, a4, b4, c0, d0;
        slong lena0, lenb0, lens, lent, lena4, lenb4, lenc0, lend0;

        /* Pointers to independently allocated memory */
        gr_ptr a2, b2, a3, b3, q, d, T0, T1;
        slong lena2, lenb2, lena3, lenb3, lenq, lend, lenT0;

        gr_ptr R[4], S[4];
        slong lenR[4], lenS[4];
        slong sgnR, sgnS;

        a2 = P;
        b2 = OFFSET(a2, lena);
        a3 = OFFSET(b2, lena);
        b3 = OFFSET(a3, lena);
        q = OFFSET(b3, lena);
        d = OFFSET(q, (lena + 1) / 2);
        T0 = OFFSET(d, lena);
        T1 = OFFSET(T0, lena);

        R[0] = OFFSET(T1, (lena + 1) / 2);
        R[1] = OFFSET(R[0], (lena + 1) / 2);
        R[2] = OFFSET(R[1], (lena + 1) / 2);
        R[3] = OFFSET(R[2], (lena + 1) / 2);
        S[0] = OFFSET(R[3], (lena + 1) / 2);
        S[1] = OFFSET(S[0], (lena + 1) / 2);
        S[2] = OFFSET(S[1], (lena + 1) / 2);
        S[3] = OFFSET(S[2], (lena + 1) / 2);

        P = OFFSET(P, 6 * lena + 10 * (lena + 1) / 2);

        __attach_shift(a0, lena0, (gr_ptr) a, lena, m);
        __attach_shift(b0, lenb0, (gr_ptr) b, lenb, m);

        if (res != NULL)
        {
            status |= gr_set(res->lc, OFFSET(b, lenb - 1), ctx);
            res->len0 -= m;
            res->len1 -= m;
            res->off += m;
        }

        if (lena0 < cutoff)
            status |= _gr_poly_hgcd_recursive_iter(&sgnR, R, lenR, &a3, &lena3,
                                                        &b3, &lenb3, a0, lena0,
                                                        b0, lenb0, q, &T0, &T1,
                                                        ctx, res);
        else
            status |= _gr_poly_hgcd_recursive(&sgnR, R, lenR, a3, &lena3, b3,
                                                   &lenb3, a0, lena0, b0,
                                                   lenb0, P, ctx, cutoff, 1, res);

        if (res != NULL)
        {
            res->off -= m;
            res->len0 += m;
            res->len1 += m;
        }

        __attach_truncate(s, lens, (gr_ptr) a, lena, m);
        __attach_truncate(t, lent, (gr_ptr) b, lenb, m);

        __mul(b2, lenb2, R[2], lenR[2], s, lens);
        __mul(T0, lenT0, R[0], lenR[0], t, lent);

        if (sgnR < 0)
            __sub(b2, lenb2, b2, lenb2, T0, lenT0);
        else
            __sub(b2, lenb2, T0, lenT0, b2, lenb2);

        status |= _gr_vec_zero(OFFSET(b2, lenb2), m + lenb3 - lenb2, ctx);

        __attach_shift(b4, lenb4, b2, lenb2, m);
        __add(b4, lenb4, b4, lenb4, b3, lenb3);
        lenb2 = FLINT_MAX(m + lenb3, lenb2);

        GR_VEC_NORM(status, b2, lenb2, sz, ctx);

        __mul(a2, lena2, R[3], lenR[3], s, lens);
        __mul(T0, lenT0, R[1], lenR[1], t, lent);

        if (sgnR < 0)
            __sub(a2, lena2, T0, lenT0, a2, lena2);
        else
            __sub(a2, lena2, a2, lena2, T0, lenT0);

        status |= _gr_vec_zero(OFFSET(a2, lena2), m + lena3 - lena2, ctx);
        __attach_shift(a4, lena4, a2, lena2, m);
        __add(a4, lena4, a4, lena4, a3, lena3);
        lena2 = FLINT_MAX(m + lena3, lena2);

        GR_VEC_NORM(status, a2, lena2, sz, ctx);

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

            *ans_sgn = sgnR;
            return status;
        }
        else
        {
            slong k = 2 * m - lenb2 + 1;

            if (res != NULL)
            {
                if (lenb2 < lenb) /* ensure something happened */
                {
                    if (lenb2 >= 1)
                    {
                        status |= gr_pow_ui(res->lc, res->lc, res->len0 - lenb2, ctx);
                        status |= gr_mul(res->res, res->res, res->lc, ctx);

                        if ((((res->len0 + res->off) | (res->len1 + res->off)) & 1) == 0)
                        {
                            status |= gr_neg(res->res, res->res, ctx);
                        }
                    }
                    else
                    {
                        if (res->len1 == 1)
                        {
                            status |= gr_pow_ui(res->lc, res->lc, res->len0 - 1, ctx);
                            status |= gr_mul(res->res, res->res, res->lc, ctx);
                        }
                        else
                        {
                            status |= gr_zero(res->res, ctx);
                        }
                    }
                }

                status |= gr_set(res->lc, OFFSET(b2, lenb2 - 1), ctx);

                res->len0 = lena2;
                res->len1 = lenb2;
            }

            __divrem(q, lenq, d, lend, a2, lena2, b2, lenb2);
            __attach_shift(c0, lenc0, b2, lenb2, k);
            __attach_shift(d0, lend0, d, lend, k);

            if (res != NULL)
            {
                if (lend >= m + 1)
                {
                    if (lend >= 1)
                    {
                        status |= gr_pow_ui(res->lc, res->lc, lena2 - lend, ctx);
                        status |= gr_mul(res->res, res->res, res->lc, ctx);

                        if ((((lena2 + res->off) | (lenb2 + res->off)) & 1) == 0)
                            status |= gr_neg(res->res, res->res, ctx);
                    }
                    else
                    {
                        if (lenb2 == 1)
                        {
                            status |= gr_pow_ui(res->lc, res->lc, lena2 - 1, ctx);
                            status |= gr_mul(res->res, res->res, res->lc, ctx);
                        }
                        else
                            status |= gr_zero(res->res, ctx);
                    }

                    res->len0 = lenb2;
                    res->len1 = lend;
               }

               res->len0 -= k;
               res->len1 -= k;
               res->off += k;
            }

            if (lenc0 < cutoff)
                status |= _gr_poly_hgcd_recursive_iter(&sgnS,
                    S, lenS, &a3, &lena3, &b3, &lenb3, c0, lenc0, d0, lend0, a2,
                    &T0, &T1, ctx, res); /* a2 as temp */
            else
                status |= _gr_poly_hgcd_recursive(&sgnS, S, lenS, a3, &lena3, b3,
                                                       &lenb3, c0, lenc0, d0,
                                                       lend0, P, ctx, cutoff, 1, res);

            if (res != NULL)
            {
                res->len0 += k;
                res->len1 += k;
                res->off -= k;
            }

            __attach_truncate(s, lens, b2, lenb2, k);
            __attach_truncate(t, lent, d, lend, k);

            __mul(B, *lenB, S[2], lenS[2], s, lens);
            __mul(T0, lenT0, S[0], lenS[0], t, lent);

            if (sgnS < 0)
                __sub(B, *lenB, B, *lenB, T0, lenT0);
            else
                __sub(B, *lenB, T0, lenT0, B, *lenB);

            status |= _gr_vec_zero(OFFSET(B, *lenB), k + lenb3 - *lenB, ctx);
            __attach_shift(b4, lenb4, B, *lenB, k);
            __add(b4, lenb4, b4, lenb4, b3, lenb3);
            *lenB = FLINT_MAX(k + lenb3, *lenB);
            GR_VEC_NORM(status, B, *lenB, sz, ctx);

            __mul(A, *lenA, S[3], lenS[3], s, lens);
            __mul(T0, lenT0, S[1], lenS[1], t, lent);

            if (sgnS < 0)
                __sub(A, *lenA, T0, lenT0, A, *lenA);
            else
                __sub(A, *lenA, A, *lenA, T0, lenT0);

            status |= _gr_vec_zero(OFFSET(A, *lenA), k + lena3 - *lenA, ctx);
            __attach_shift(a4, lena4, A, *lenA, k);
            __add(a4, lena4, a4, lena4, a3, lena3);
            *lenA = FLINT_MAX(k + lena3, *lenA);
            GR_VEC_NORM(status, A, *lenA, sz, ctx);

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

            *ans_sgn = -(sgnR * sgnS);
            return status;
        }
    }
}

/*
    XXX: Currently supports aliasing between {A,a} and {B,b}.
 */

int _gr_poly_hgcd(gr_ptr r, slong * sgn, gr_ptr * M, slong * lenM,
                               gr_ptr A, slong * lenA,
                               gr_ptr B, slong * lenB,
                               gr_srcptr a, slong lena,
                               gr_srcptr b, slong lenb,
                               slong cutoff,
                               gr_ctx_t ctx)
{
    slong lenW = 22 * lena + 16 * (FLINT_CLOG2(lena) + 1);
    slong sgnM;
    gr_ptr W;
    int status = GR_SUCCESS;
    gr_poly_res_t res;
    slong sz = ctx->sizeof_elem;

    if (lena == 0 || lenb == 0)
    {
        if (sgn != NULL)
            *sgn = 0;
        *lenA = *lenB = 0;
        if (lenM != NULL)
            lenM[0] = lenM[1] = lenM[2] = lenM[3] = 0;
        return GR_DOMAIN;
    }

    if (r != NULL)
    {
        GR_TMP_INIT2(res->res, res->lc, ctx);
        status |= gr_set(res->res, r, ctx);
        status |= gr_set(res->lc, GR_ENTRY(b, lenb - 1, sz), ctx);
        res->len0 = lena;
        res->len1 = lenb;
        res->off = 0;
    }

    if (lenb < lena / 2 + 1)
        lenW = 0;

    GR_TMP_INIT_VEC(W, lenW, ctx);

    if (M == NULL)
    {
        status = _gr_poly_hgcd_recursive(&sgnM, NULL, NULL, A, lenA, B, lenB, a, lena, b, lenb, W, ctx, cutoff, 0, r == NULL ? NULL : res);
    }
    else
    {
        status = _gr_poly_hgcd_recursive(&sgnM, M, lenM, A, lenA, B, lenB, a, lena, b, lenb, W, ctx, cutoff, 1, r == NULL ? NULL : res);
    }

    if (r != NULL)
    {
        if (*lenB < lenb) /* make sure something happened */
        {
            if (*lenB >= 1)
            {
                status |= gr_pow_ui(res->lc, res->lc, res->len0 - *lenB, ctx);
                status |= gr_mul(res->res, res->res, res->lc, ctx);

                if (((res->len0 | res->len1) & 1) == 0)
                    status |= gr_neg(res->res, res->res, ctx);
            }
            else
            {
                if (res->len1 == 1)
                {
                    status |= gr_pow_ui(res->lc, res->lc, res->len0 - 1, ctx);
                    status |= gr_mul(res->res, res->res, res->lc, ctx);
                }
                else
                    status |= gr_zero(res->res, ctx);
            }
        }

        status |= gr_set(r, res->res, ctx);

        GR_TMP_CLEAR2(res->res, res->lc, ctx);
    }

    GR_TMP_CLEAR_VEC(W, lenW, ctx);

    if (sgn != NULL)
        *sgn = sgnM;

    return status;
}
