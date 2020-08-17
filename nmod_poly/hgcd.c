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

static __inline__ void __mat_one(mp_ptr *M, slong *lenM)
{
    M[0][0] = WORD(1);
    M[3][0] = WORD(1);
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

static void __mat_mul_classical(mp_ptr *C, slong *lenC, 
    mp_ptr *A, slong *lenA, mp_ptr *B, slong *lenB, mp_ptr T, nmod_t mod)
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

static void __mat_mul_strassen(mp_ptr *C, slong *lenC, 
    mp_ptr *A, slong *lenA, mp_ptr *B, slong *lenB, mp_ptr T0, mp_ptr T1, 
    nmod_t mod)
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

static void __mat_mul(mp_ptr *C, slong *lenC, 
    mp_ptr *A, slong *lenA, mp_ptr *B, slong *lenB, mp_ptr T0, mp_ptr T1, 
    nmod_t mod)
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

slong _nmod_poly_hgcd_recursive_iter(mp_ptr *M, slong *lenM, 
    mp_ptr *A, slong *lenA, mp_ptr *B, slong *lenB, 
    mp_srcptr a, slong lena, mp_srcptr b, slong lenb, 
    mp_ptr Q, mp_ptr *T, mp_ptr *t, nmod_t mod, nmod_poly_res_t res)
{
    const slong m = lena / 2;
    slong sgn = 1;

    __mat_one(M, lenM);
    __set(*A, *lenA, a, lena);
    __set(*B, *lenB, b, lenb);

    while (*lenB >= m + 1)
    {
        slong lenQ, lenT, lent;
        
        if (res)
           res->lc = (*B)[*lenB - 1];

        __divrem(Q, lenQ, *T, lenT, *A, *lenA, *B, *lenB);

        if (res)
        {
           if (lenT >= m + 1)
           {
              if (lenT >= 1)
              {
                 res->lc  = n_powmod2_preinv(res->lc, *lenA - lenT, mod.n, mod.ninv);
                 res->res = n_mulmod2_preinv(res->res, res->lc, mod.n, mod.ninv);
              
                 if ((((*lenA + res->off) | (*lenB + res->off)) & 1) == 0)
                    res->res = nmod_neg(res->res, mod);
              } else
              {
                 if (*lenB == 1) 
                 {
                    res->lc  = n_powmod2_preinv(res->lc, *lenA - 1, mod.n, mod.ninv);
                    res->res = n_mulmod2_preinv(res->res, res->lc, mod.n, mod.ninv);
                 } else
                    res->res = 0;
              }
           } else
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

    The res struct, if not NULL, passes information required to compute
    the sign changes and powers of leading terms used to compute the
    resultant.
 */

slong _nmod_poly_hgcd_recursive(mp_ptr *M, slong *lenM, 
    mp_ptr A, slong *lenA, mp_ptr B, slong *lenB, 
    mp_srcptr a, slong lena, mp_srcptr b, slong lenb, 
    mp_ptr P, nmod_t mod, int flag, nmod_poly_res_t res)
{
    const slong m = lena / 2;

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
        slong lena0, lenb0, lens, lent, lena4, lenb4, lenc0, lend0;

        /* Pointers to independently allocated memory */
        mp_ptr a2, b2, a3, b3, q, d, T0, T1;
        slong lena2, lenb2, lena3, lenb3, lenq, lend, lenT0;

        mp_ptr R[4], S[4];
        slong lenR[4], lenS[4];
        slong sgnR, sgnS;

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

        if (res)
        {
           res->lc = b[lenb - 1];

           res->len0 -= m;
           res->len1 -= m;
           res->off += m;
        }

        if (lena0 < NMOD_POLY_HGCD_CUTOFF)
            sgnR = _nmod_poly_hgcd_recursive_iter(R, lenR, &a3, &lena3, &b3, &lenb3, 
                                            a0, lena0, b0, lenb0, 
                                            q, &T0, &T1, mod, res);
        else 
            sgnR = _nmod_poly_hgcd_recursive(R, lenR, a3, &lena3, b3, &lenb3, 
                                       a0, lena0, b0, lenb0, P, mod, 1, res);

        if (res)
        {
           res->off -= m;
           res->len0 += m;
           res->len1 += m;
        }

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
            slong k = 2 * m - lenb2 + 1;

            if (res) 
            {
               if (lenb2 < lenb) /* ensure something happened */
               {
                  if (lenb2 >= 1)
                  {
                     res->lc  = n_powmod2_preinv(res->lc, res->len0 - lenb2, mod.n, mod.ninv);
                     res->res = n_mulmod2_preinv(res->res, res->lc, mod.n, mod.ninv);

                     if ((((res->len0 + res->off) | (res->len1 + res->off)) & 1) == 0)
                        res->res = nmod_neg(res->res, mod);
                  } else
                  {
                     if (res->len1 == 1) 
                     {
                        res->lc  = n_powmod2_preinv(res->lc, res->len0 - 1, mod.n, mod.ninv);
                        res->res = n_mulmod2_preinv(res->res, res->lc, mod.n, mod.ninv);
                     } else
                        res->res = 0;
                  }
               }

               res->lc = b2[lenb2 - 1];
            
               res->len0 = lena2;
               res->len1 = lenb2;
            }

            __divrem(q, lenq, d, lend, a2, lena2, b2, lenb2);

            __attach_shift(c0, lenc0, b2, lenb2, k);
            __attach_shift(d0, lend0, d, lend, k);
            
            if (res)
            {
               if (lend >= m + 1)
               {
                  if (lend >= 1)
                  {
                     res->lc  = n_powmod2_preinv(res->lc, lena2 - lend, mod.n, mod.ninv);
                     res->res = n_mulmod2_preinv(res->res, res->lc, mod.n, mod.ninv);

                     if ((((lena2 + res->off) | (lenb2 + res->off)) & 1) == 0)
                        res->res = nmod_neg(res->res, mod);
                  } else
                  {
                     if (lenb2 == 1) 
                     {
                        res->lc  = n_powmod2_preinv(res->lc, lena2 - 1, mod.n, mod.ninv);
                        res->res = n_mulmod2_preinv(res->res, res->lc, mod.n, mod.ninv);
                     } else
                        res->res = 0;
                  }
                  
                  res->len0 = lenb2;
                  res->len1 = lend;
               }

               res->len0 -= k;
               res->len1 -= k;
               res->off += k;
            } 
            
            if (lenc0 < NMOD_POLY_HGCD_CUTOFF)
                sgnS = _nmod_poly_hgcd_recursive_iter(S, lenS, &a3, &lena3, &b3, &lenb3, 
                                                c0, lenc0, d0, lend0, 
                                                a2, &T0, &T1, mod, res); /* a2 as temp */
            else 
                sgnS = _nmod_poly_hgcd_recursive(S, lenS, a3, &lena3, b3, &lenb3, 
                                           c0, lenc0, d0, lend0, P, mod, 1, res);

            if (res)
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

slong _nmod_poly_hgcd(mp_ptr *M, slong *lenM, 
                     mp_ptr A, slong *lenA, mp_ptr B, slong *lenB, 
                     mp_srcptr a, slong lena, mp_srcptr b, slong lenb, 
                     nmod_t mod)
{
    const slong lenW = 22 * lena + 16 * (FLINT_CLOG2(lena) + 1);
    slong sgnM;
    mp_ptr W;
    
    W = _nmod_vec_init(lenW);

    if (M == NULL)
    {
        sgnM = _nmod_poly_hgcd_recursive(NULL, NULL, 
                                         A, lenA, B, lenB, 
                                         a, lena, b, lenb, W, mod, 0, NULL);
    }
    else
    {
        sgnM = _nmod_poly_hgcd_recursive(M, lenM, 
                                         A, lenA, B, lenB, 
                                         a, lena, b, lenb, W, mod, 1, NULL);
    }
    _nmod_vec_clear(W);

    return sgnM;
}


/*
    Assuming deg(a) > deg(b) (in particular b could be zero, but a must not be)
    compute a matrix M = [[m11, m12] [m21 m22]] and A, B so that

    [ a ] = [ m11  m12 ][ A ]
    [ b ]   [ m21  m22 ][ B ]

    with

    (1) A and B are consecutive remainders in the euclidean remainder
        sequence for a, b satsifying 2*deg(A) >= deg(a) > 2*deg(B)

    (2) M is a product of [[qi 1][1 0]] where the qi are the quotients
        obtained in (1)

    No aliasing. The return is det(M), which is +-1. All of the moduli of the
    arguments should be the same and prime. Here is a reference implementation
    in case something is wrong with _nmod_poly_hgcd (there doesn't seem to be).
*/
slong nmod_poly_hgcd_ref(
    nmod_poly_t m11, nmod_poly_t m12,
    nmod_poly_t m21, nmod_poly_t m22,
    nmod_poly_t A, nmod_poly_t B,
    const nmod_poly_t a, const nmod_poly_t b)
{
    slong sgnM;
    slong dega = nmod_poly_degree(a);
    nmod_poly_t q, r, t;

    if (nmod_poly_degree(a) <= nmod_poly_degree(b))
    {
        flint_throw(FLINT_ERROR, "Exception in nmod_poly_hgcd_ref:"
                                              " Input degrees are invalid.\n");
    }

    nmod_poly_init_mod(q, a->mod);
    nmod_poly_init_mod(r, a->mod);
    nmod_poly_init_mod(t, a->mod);

    nmod_poly_one(m11);
    nmod_poly_zero(m12);
    nmod_poly_zero(m21);
    nmod_poly_one(m22);
    nmod_poly_set(A, a);
    nmod_poly_set(B, b);

    sgnM = 1;
    while (dega <= 2*nmod_poly_degree(B))
    {
        nmod_poly_divrem(q, r, A, B);
        nmod_poly_swap(A, B);
        nmod_poly_swap(B, r);

        /* multipliy M by [[q 1] [1 0]] on the right */
        nmod_poly_mul(t, q, m11);
        nmod_poly_add(r, m12, t);
        nmod_poly_swap(m11, m12);
        nmod_poly_swap(m11, r);

        nmod_poly_mul(t, q, m21);
        nmod_poly_add(r, m22, t);
        nmod_poly_swap(m21, m22);
        nmod_poly_swap(m21, r);

        sgnM = -sgnM;
    }

    nmod_poly_clear(q);
    nmod_poly_clear(r);
    nmod_poly_clear(t);

    return sgnM;
}


slong nmod_poly_hgcd(
    nmod_poly_t m11, nmod_poly_t m12,
    nmod_poly_t m21, nmod_poly_t m22,
    nmod_poly_t A, nmod_poly_t B,
    const nmod_poly_t a, const nmod_poly_t b)
{
    mp_limb_t * M[4];
    slong lenM[4];
    slong sgnM;

    if (nmod_poly_degree(a) <= nmod_poly_degree(b))
    {
        flint_throw(FLINT_ERROR, "Exception in nmod_poly_hgcd:"
                                              " Input degrees are invalid.\n");
    }

    if (nmod_poly_length(b) == 0)
    {
        nmod_poly_one(m11);
        nmod_poly_zero(m12);
        nmod_poly_zero(m21);
        nmod_poly_one(m22);
        nmod_poly_set(A, a);
        nmod_poly_set(B, b);
        return 1;
    }

    /* now nmod_poly_length(a) > nmod_poly_length(b) > 0 */

    nmod_poly_fit_length(m11, nmod_poly_length(a));
    nmod_poly_fit_length(m12, nmod_poly_length(a));
    nmod_poly_fit_length(m21, nmod_poly_length(a));
    nmod_poly_fit_length(m22, nmod_poly_length(a));
    nmod_poly_fit_length(A, nmod_poly_length(a));
    nmod_poly_fit_length(B, nmod_poly_length(a));

    /*
        Looks like _nmod_poly_hgcd produces an M with det(M) = sgnM = +-1 and
            [ a ] = [ m11  m12 ] [ A ]
            [ b ]   [ m21  m22 ] [ B ]
    */

    M[0] = m11->coeffs;
    M[1] = m12->coeffs;
    M[2] = m21->coeffs;
    M[3] = m22->coeffs;
    sgnM = _nmod_poly_hgcd(M, lenM, A->coeffs, &A->length, B->coeffs, &B->length,
                          a->coeffs, a->length, b->coeffs, b->length, A->mod);
    m11->length = lenM[0];
    m12->length = lenM[1];
    m21->length = lenM[2];
    m22->length = lenM[3];

    FLINT_ASSERT(2*nmod_poly_degree(A) >= nmod_poly_degree(a));
    FLINT_ASSERT(nmod_poly_degree(a) > 2*nmod_poly_degree(B));

    return sgnM;
}
