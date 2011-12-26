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
#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"
#include "mpn_extras.h"

#define NMOD_POLY_HGCD_CUTOFF 60        /* HGCD: Basecase -> Recursion */
#define NMOD_POLY_GCD_CUTOFF 184        /* GCD:  Euclidean -> HGCD */
#define NMOD_POLY_SMALL_GCD_CUTOFF 174  /* GCD (small n): Euclidean -> HGCD */

static int FLAG = 0;

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

#define __mat_one(M, lenM)\
do {\
    (M)[0][0] = 1L;\
    (M)[1][0] = 0L;\
    (M)[2][0] = 0L;\
    (M)[3][0] = 1L;\
    (lenM)[0] = 1;\
    (lenM)[1] = 0;\
    (lenM)[2] = 0;\
    (lenM)[3] = 1;\
} while (0)

static __inline__ void __mat_mul(mp_ptr *C, long *lenC, 
    mp_ptr *A, long *lenA, mp_ptr *B, long *lenB, const nmod_t mod)
{
    mp_ptr T;
    long lenT = lenA[0];

    {
        lenT = FLINT_MAX(lenT, lenA[1]);
        lenT = FLINT_MAX(lenT, lenA[2]);
        lenT = FLINT_MAX(lenT, lenA[3]);
        lenT = FLINT_MAX(lenT, lenB[0]);
        lenT = FLINT_MAX(lenT, lenB[1]);
        lenT = FLINT_MAX(lenT, lenB[2]);
        lenT = FLINT_MAX(lenT, lenB[3]);
    }

    T = _nmod_vec_init(2 * lenT - 1);

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

    _nmod_vec_clear(T);
}

/*
    Only supports aliasing in {A,a} and {B,b}.
 */

long _nmod_poly_half_gcd_iter(mp_ptr *M, long *lenM, 
    mp_ptr *A, long *lenA, mp_ptr *B, long *lenB, 
    mp_srcptr a, long lena, mp_srcptr b, long lenb, 
    mp_ptr *Q, mp_ptr *T, const nmod_t mod)
{
    const long m = lena / 2;
    long sgn = 1;

    __mat_one(M, lenM);
    __set(*A, *lenA, a, lena);
    __set(*B, *lenB, b, lenb);

    while (*lenB >= m + 1)
    {
        long lenQ, lenT;

        __divrem(*Q, lenQ, *T, lenT, *A, *lenA, *B, *lenB);
        __swap(*B, *lenB, *T, lenT);
        __swap(*A, *lenA, *T, lenT);

        __mul(*T, lenT, *Q, lenQ, M[2], lenM[2]);
        __add(*T, lenT, M[3], lenM[3], *T, lenT);
        __swap(M[3], lenM[3], M[2], lenM[2]);
        __swap(M[2], lenM[2], *T, lenT);

        __mul(*T, lenT, *Q, lenQ, M[0], lenM[0]);
        __add(*T, lenT, M[1], lenM[1], *T, lenT);
        __swap(M[1], lenM[1], M[0], lenM[0]);
        __swap(M[0], lenM[0], *T, lenT);

        sgn = -sgn;
    }

    return sgn;
}

long nmod_poly_half_gcd_iter(nmod_poly_mat_t res, 
    nmod_poly_t a2, nmod_poly_t b2, const nmod_poly_t a, const nmod_poly_t b)
{
    const long m = a->length / 2;

    nmod_poly_mat_one(res);
    nmod_poly_set(a2, a);
    nmod_poly_set(b2, b);

    if (b->length < m + 1) 
    {
        return 1L;
    }
    else
    {
        nmod_poly_t Q, temp;
        long sign = 1L;

        nmod_poly_init_preinv(Q, a->mod.n, a->mod.ninv);
        nmod_poly_init_preinv(temp, a->mod.n, a->mod.ninv);

        while (b2->length >= m + 1)
        {
            nmod_poly_divrem(Q, a2, a2, b2);
            nmod_poly_swap(a2, b2);

            nmod_poly_mul(temp, Q, nmod_poly_mat_entry(res,1,0));
            nmod_poly_add(temp, nmod_poly_mat_entry(res,1,1), temp);
            nmod_poly_swap(nmod_poly_mat_entry(res,1,1), nmod_poly_mat_entry(res,1,0));
            nmod_poly_swap(nmod_poly_mat_entry(res,1,0), temp);

            nmod_poly_mul(temp, Q, nmod_poly_mat_entry(res,0,0));
            nmod_poly_add(temp, nmod_poly_mat_entry(res,0,1), temp);
            nmod_poly_swap(nmod_poly_mat_entry(res,0,1), nmod_poly_mat_entry(res,0,0));
            nmod_poly_swap(nmod_poly_mat_entry(res,0,0), temp);
            sign = -sign;
        }

        nmod_poly_clear(temp);
        nmod_poly_clear(Q);

        return sign;
    }
}


static __inline__ 
void _nmod_poly_attach(nmod_poly_t output, nmod_poly_t input)
{
   output->length = input->length;
   output->coeffs = input->coeffs;
   output->mod = input->mod;
}

static __inline__ 
void nmod_poly_attach(nmod_poly_t output, nmod_poly_t input)
{
   _nmod_poly_attach(output, input);
}

/*
   Attach input shifted right by n to output
*/

static __inline__ 
void _nmod_poly_attach_shift(nmod_poly_t output, 
                   nmod_poly_t input, long n)
{
   if (input->length >= n) output->length = input->length - n;
   else output->length = 0;
   output->coeffs = input->coeffs + n;
   output->mod = input->mod;
}

static __inline__ 
void nmod_poly_attach_shift(nmod_poly_t output, 
                   nmod_poly_t input, long n)
{
   _nmod_poly_attach_shift(output, input, n);
}

/*
   Attach input to first n coefficients of input
*/

static __inline__ 
void _nmod_poly_attach_truncate(nmod_poly_t output, 
                     nmod_poly_t input, long n)
{
   if (input->length < n) output->length = input->length;
   else output->length = n;
   output->coeffs = input->coeffs;
   output->mod = input->mod;
}

static __inline__ 
void nmod_poly_attach_truncate(nmod_poly_t output, 
                     nmod_poly_t input, long n)
{
   _nmod_poly_attach_truncate(output, input, n);
}


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

/*
    Assumes that lena >= lenb.
 */

#define trace_vec(str, vec, len, mod)\
do {\
if (FLAG) {\
printf(str), _nmod_vec_print(vec, len, mod), printf("\n");\
fflush(stdout);\
}\
} while (0)

#define trace_poly(str, poly)\
do{\
if (FLAG){\
printf(str), nmod_poly_print(poly), printf("\n");\
fflush(stdout);}\
}while(0)

long _nmod_poly_half_gcd(mp_ptr *M, long *lenM, 
    mp_ptr *A, long *lenA, mp_ptr *B, long *lenB, 
    mp_srcptr a, long lena, mp_srcptr b, long lenb, 
    const nmod_t mod, int flag)
{
    const long m = lena / 2;

if (FLAG) {
printf("---{%ld}{%ld}-----------------------------------------\n", lena, lenb); fflush(stdout);
}
trace_vec("a = ", a, lena, mod);
trace_vec("b = ", b, lenb, mod);

    if (lenb < m + 1)
    {
        __mat_one(M, lenM);
        __set(*A, *lenA, a, lena);
        __set(*B, *lenB, b, lenb);
        return 1;
    }
    else
    {
        long sgnR, sgnS;

        const long lenX = 2 * lena;
        mp_ptr X = _nmod_vec_init(16 * lenX);

        mp_ptr a0, b0, s, t, a4, b4, c0, d0;
        long lena0, lenb0, lens, lent, lena4, lenb4, lenc0, lend0;

        mp_ptr a2, b2, a3, b3, q, d, temp, temp1;
        long lena2, lenb2, lena3, lenb3, lenq, lend, lentemp, lentemp1;

        mp_ptr R[4], S[4];
        long lenR[4], lenS[4];

        a2 = X;
        b2 = a2 + lenX;
        a3 = b2 + lenX;
        b3 = a3 + lenX;
        q  = b3 + lenX;
        d  = q  + lenX;
        temp = d  + lenX;
        temp1 = temp + lenX;

        R[0] = temp1 + lenX;
        R[1] = R[0] + lenX;
        R[2] = R[1] + lenX;
        R[3] = R[2] + lenX;
        S[0] = R[3] + lenX;
        S[1] = S[0] + lenX;
        S[2] = S[1] + lenX;
        S[3] = S[2] + lenX;

        __attach_shift(a0, lena0, (mp_ptr) a, lena, m);
        __attach_shift(b0, lenb0, (mp_ptr) b, lenb, m);

if (FLAG) {
printf("a0 = "), _nmod_vec_print(a0, lena0,mod), printf("\n");
printf("b0 = "), _nmod_vec_print(b0, lenb0,mod), printf("\n");
printf("lena0 < CUTOFF ? %ld\n", lena0 < NMOD_POLY_HGCD_CUTOFF);
}

        if (lena0 < NMOD_POLY_HGCD_CUTOFF)
            sgnR = _nmod_poly_half_gcd_iter(R, lenR, &a3, &lena3, &b3, &lenb3, 
                                            a0, lena0, b0, lenb0, 
                                            &temp, &temp1, mod);
        else 
            sgnR = _nmod_poly_half_gcd(R, lenR, &a3, &lena3, &b3, &lenb3, 
                                       a0, lena0, b0, lenb0, mod, 1);
if (FLAG) {
printf("sgnR = %ld\n", sgnR);
printf("R[0] = "), _nmod_vec_print(R[0], lenR[0],mod), printf("\n");
printf("R[1] = "), _nmod_vec_print(R[1], lenR[1],mod), printf("\n");
printf("R[2] = "), _nmod_vec_print(R[2], lenR[2],mod), printf("\n");
printf("R[3] = "), _nmod_vec_print(R[3], lenR[3],mod), printf("\n");
}

        __attach_truncate(s, lens, (mp_ptr) a, lena, m);
        __attach_truncate(t, lent, (mp_ptr) b, lenb, m);

        __mul(b2, lenb2, R[2], lenR[2], s, lens);
        __mul(temp, lentemp, R[0], lenR[0], t, lent);

        if (sgnR < 0)
            __sub(b2, lenb2, b2, lenb2, temp, lentemp);
        else
            __sub(b2, lenb2, temp, lentemp, b2, lenb2);

        if (m + lenb3 > lenX)
        {
            printf("ERROR. NOT ENOUGH SPACE\n"); 
            abort();
        }

        mpn_zero(b2 + lenb2, m + lenb3 - lenb2);

        __attach_shift(b4, lenb4, b2, lenb2, m);
        __add(b4, lenb4, b4, lenb4, b3, lenb3);
        lenb2 = FLINT_MAX(m + lenb3, lenb2);
        MPN_NORM(b2, lenb2);

        __mul(a2, lena2, R[3], lenR[3], s, lens);
        __mul(temp, lentemp, R[1], lenR[1], t, lent);

        if (sgnR < 0)
            __sub(a2, lena2, temp, lentemp, a2, lena2);
        else
            __sub(a2, lena2, a2, lena2, temp, lentemp);

        if (m + lena3 > lenX)
        {
            printf("CRAZY!!\n"); abort();
        }

        mpn_zero(a2 + lena2, m + lena3 - lena2);
        __attach_shift(a4, lena4, a2, lena2, m);
        __add(a4, lena4, a4, lena4, a3, lena3);
        lena2 = FLINT_MAX(m + lena3, lena2);
        MPN_NORM(a2, lena2);

        if (lenb2 < m + 1)
        {
            __set(*A, *lenA, a2, lena2);
            __set(*B, *lenB, b2, lenb2);

            __set(M[0], lenM[0], R[0], lenR[0]);
            __set(M[1], lenM[1], R[1], lenR[1]);
            __set(M[2], lenM[2], R[2], lenR[2]);
            __set(M[3], lenM[3], R[3], lenR[3]);

            return sgnR;
        }

        __divrem(q, lenq, d, lend, a2, lena2, b2, lenb2);

trace_vec("q = ", q, lenq, mod);
trace_vec("d = ", d, lend, mod);
trace_vec("a2 = ", a2, lena2, mod);
trace_vec("b2 = ", b2, lenb2, mod);

        long k = 2 * m - lenb2 + 1;

        __attach_shift(c0, lenc0, b2, lenb2, k);
        __attach_shift(d0, lend0, d, lend, k);

        if (lenc0 < NMOD_POLY_HGCD_CUTOFF)
            sgnS = _nmod_poly_half_gcd_iter(S, lenS, &a3, &lena3, &b3, &lenb3, 
                                            c0, lenc0, d0, lend0, 
                                            &temp, &temp1, mod);
        else 
            sgnS = _nmod_poly_half_gcd(S, lenS, &a3, &lena3, &b3, &lenb3, 
                                       c0, lenc0, d0, lend0, mod, 1);

        __attach_truncate(s, lens, b2, lenb2, k);
        __attach_truncate(t, lent, d, lend, k);

        __mul(*B, *lenB, S[2], lenS[2], s, lens);
        __mul(temp, lentemp, S[0], lenS[0], t, lent);

        if (sgnS < 0)
            __sub(*B, *lenB, *B, *lenB, temp, lentemp);
        else
            __sub(*B, *lenB, temp, lentemp, *B, *lenB);

        mpn_zero(*B + *lenB, k + lenb3 - *lenB);
        __attach_shift(b4, lenb4, *B, *lenB, k);
        __add(b4, lenb4, b4, lenb4, b3, lenb3);
        *lenB = FLINT_MAX(k + lenb3, *lenB);
        MPN_NORM(*B, *lenB);

        __mul(*A, *lenA, S[3], lenS[3], s, lens);
        __mul(temp, lentemp, S[1], lenS[1], t, lent);

        if (sgnS < 0)
            __sub(*A, *lenA, temp, lentemp, *A, *lenA);
        else
            __sub(*A, *lenA, *A, *lenA, temp, lentemp);

        mpn_zero(*A + *lenA, k + lena3 - *lenA);
        __attach_shift(a4, lena4, *A, *lenA, k);
        __add(a4, lena4, a4, lena4, a3, lena3);
        *lenA = FLINT_MAX(k + lena3, *lenA);
        MPN_NORM(*A, *lenA);

        if (flag)
        {
            __swap(S[0], lenS[0], S[2], lenS[2]);
            __swap(S[1], lenS[1], S[3], lenS[3]);
            __mul(temp, lentemp, S[2], lenS[2], q, lenq);
            __add(S[0], lenS[0], S[0], lenS[0], temp, lentemp);
            __mul(temp, lentemp, S[3], lenS[3], q, lenq);
            __add(S[1], lenS[1], S[1], lenS[1], temp, lentemp);

            __mat_mul(M, lenM, R, lenR, S, lenS, mod);
        }

        return - (sgnR * sgnS);
    }
}

long nmod_poly_half_gcd(nmod_poly_mat_t R, nmod_poly_t A, nmod_poly_t B, 
    const nmod_poly_t a, const nmod_poly_t b, int flag)
{
    long max = FLINT_MAX(a->length, b->length);

    int ans;

    mp_ptr M[4];
    long lenM[4] = {0,0,0,0};

    M[0] = _nmod_vec_init(max);
    M[1] = _nmod_vec_init(max);
    M[2] = _nmod_vec_init(max);
    M[3] = _nmod_vec_init(max);

    nmod_poly_fit_length(A, max);
    nmod_poly_fit_length(B, max);

    ans =  _nmod_poly_half_gcd(M, lenM, &(A->coeffs), &(A->length), &(B->coeffs), &(B->length), 
                                        a->coeffs, a->length, b->coeffs, b->length, 
                                        a->mod, flag);

    nmod_poly_fit_length(nmod_poly_mat_entry(R,0,0), max);
    nmod_poly_fit_length(nmod_poly_mat_entry(R,0,1), max);
    nmod_poly_fit_length(nmod_poly_mat_entry(R,1,0), max);
    nmod_poly_fit_length(nmod_poly_mat_entry(R,1,1), max);

if (lenM[0] > max) {
printf("ABORT...BLAH\n"); abort();
}

    _nmod_vec_set(nmod_poly_mat_entry(R,0,0)->coeffs, M[0], lenM[0]);
    _nmod_vec_set(nmod_poly_mat_entry(R,0,1)->coeffs, M[1], lenM[1]);
    _nmod_vec_set(nmod_poly_mat_entry(R,1,0)->coeffs, M[2], lenM[2]);
    _nmod_vec_set(nmod_poly_mat_entry(R,1,1)->coeffs, M[3], lenM[3]);

    nmod_poly_mat_entry(R, 0, 0)->length = lenM[0];
    nmod_poly_mat_entry(R, 0, 1)->length = lenM[1];
    nmod_poly_mat_entry(R, 1, 0)->length = lenM[2];
    nmod_poly_mat_entry(R, 1, 1)->length = lenM[3];

    return ans;
}


long nmod_poly_half_gcd2(nmod_poly_mat_t R, nmod_poly_t A, nmod_poly_t B, 
    const nmod_poly_t a, const nmod_poly_t b, int flag)
{
    const long m = a->length/2;

if (FLAG) {
printf("---{%ld}{%ld}-----------------------------------------\n", a->length, b->length); fflush(stdout);
printf("a = "), nmod_poly_print(a); printf("\n");
printf("b = "), nmod_poly_print(b); printf("\n");
}

    if (b->length < m + 1)
    {
        nmod_poly_mat_one(R);
        nmod_poly_set(A, a);
        nmod_poly_set(B, b);
        return 1L;
    }
    else
    {
        long signR;
        nmod_poly_t a0, b0;
        nmod_poly_t temp, a2, b2; 
        nmod_poly_t a3, b3, a4, b4, s, t;

        _nmod_poly_attach_shift(a0, a, m);
        _nmod_poly_attach_shift(b0, b, m);

        nmod_poly_init_preinv(temp, a->mod.n, a->mod.ninv);
        nmod_poly_init_preinv(a2, a->mod.n, a->mod.ninv);
        nmod_poly_init_preinv(b2, a->mod.n, a->mod.ninv);

        nmod_poly_init_preinv(a3, a->mod.n, a->mod.ninv);
        nmod_poly_init_preinv(b3, a->mod.n, a->mod.ninv);

if (FLAG) {
printf("a0 = "), nmod_poly_print(a0), printf("\n");
printf("b0 = "), nmod_poly_print(b0), printf("\n");
printf("lena0 < CUTOFF ? %ld\n", a0->length < NMOD_POLY_HGCD_CUTOFF);
}

        if (a0->length < NMOD_POLY_HGCD_CUTOFF) 
            signR = nmod_poly_half_gcd_iter(R, a3, b3, a0, b0);
        else 
            signR = nmod_poly_half_gcd2(R, a3, b3, a0, b0, 1);

if (FLAG) {
printf("sgnR = %ld\n", signR);
printf("R[0] = "), nmod_poly_print(nmod_poly_mat_entry(R,0,0)), printf("\n");
printf("R[1] = "), nmod_poly_print(nmod_poly_mat_entry(R,0,1)), printf("\n");
printf("R[2] = "), nmod_poly_print(nmod_poly_mat_entry(R,1,0)), printf("\n");
printf("R[3] = "), nmod_poly_print(nmod_poly_mat_entry(R,1,1)), printf("\n");
}

        nmod_poly_attach_truncate(s, a, m);
        nmod_poly_attach_truncate(t, b, m);

        nmod_poly_mul(b2, nmod_poly_mat_entry(R,1,0), s);
        nmod_poly_mul(temp, nmod_poly_mat_entry(R,0,0), t);

        if (signR < 0L) 
            nmod_poly_sub(b2, b2, temp);
        else 
            nmod_poly_sub(b2, temp, b2);

        nmod_poly_fit_length(b2, m + b3->length);

        long i;
        for (i = b2->length; i < m + b3->length; i++) 
            b2->coeffs[i] = 0L;

        nmod_poly_attach_shift(b4, b2, m);
        b4->alloc = FLINT_MAX(b4->length, b3->length);
        nmod_poly_add(b4, b4, b3);
        b2->length = FLINT_MAX(m + b3->length, b2->length);
        _nmod_poly_normalise(b2);

        nmod_poly_mul(a2, nmod_poly_mat_entry(R,1,1), s);
        nmod_poly_mul(temp, nmod_poly_mat_entry(R,0,1), t);

        if (signR < 0L) 
            nmod_poly_sub(a2, temp, a2);
        else 
            nmod_poly_sub(a2, a2, temp);

        nmod_poly_fit_length(a2, m + a3->length);
        for (i = a2->length; i < m + a3->length; i++) 
            a2->coeffs[i] = 0L;
        nmod_poly_attach_shift(a4, a2, m);
        a4->alloc = FLINT_MAX(a4->length, a3->length);
        nmod_poly_add(a4, a4, a3);
        a2->length = FLINT_MAX(m + a3->length, a2->length);
        _nmod_poly_normalise(a2);

        if (b2->length < m + 1)
        {
            nmod_poly_set(A, a2);
            nmod_poly_set(B, b2);
            nmod_poly_clear(temp);
            nmod_poly_clear(a2);
            nmod_poly_clear(b2);
            nmod_poly_clear(a3);
            nmod_poly_clear(b3);
            return signR;
        }

        nmod_poly_t q, d;
        nmod_poly_init(q, a->mod.n);
        nmod_poly_init(d, a->mod.n);

        nmod_poly_divrem(q, d, a2, b2);

trace_poly("q = ", q);
trace_poly("d = ", d);
trace_poly("a2 = ", a2);
trace_poly("b2 = ", b2);

        long k = 2*m - b2->length + 1;

        nmod_poly_t c0, d0;
        _nmod_poly_attach_shift(c0, b2, k);
        _nmod_poly_attach_shift(d0, d, k);
        nmod_poly_mat_t S;
        nmod_poly_mat_init(S, 2, 2, a->mod.n);

        long signS;

        if (c0->length < NMOD_POLY_HGCD_CUTOFF) 
            signS = nmod_poly_half_gcd_iter(S, a3, b3, c0, d0);
        else 
            signS = nmod_poly_half_gcd2(S, a3, b3, c0, d0, 1);

        nmod_poly_attach_truncate(s, b2, k);
        nmod_poly_attach_truncate(t, d, k);

        nmod_poly_mul(B, nmod_poly_mat_entry(S,1,0), s);
        nmod_poly_mul(temp, nmod_poly_mat_entry(S,0,0), t);

        if (signS < 0L) 
            nmod_poly_sub(B, B, temp);
        else 
            nmod_poly_sub(B, temp, B);

        nmod_poly_fit_length(B, k + b3->length);
        for (i = B->length; i < k + b3->length; i++) 
            B->coeffs[i] = 0L;
        nmod_poly_attach_shift(b4, B, k);
        b4->alloc = FLINT_MAX(b4->length, b3->length);
        nmod_poly_add(b4, b4, b3);
        B->length = FLINT_MAX(k + b3->length, B->length);
        _nmod_poly_normalise(B);

        nmod_poly_mul(A, nmod_poly_mat_entry(S,1,1), s);
        nmod_poly_mul(temp, nmod_poly_mat_entry(S,0,1), t);

        if (signS < 0L) 
            nmod_poly_sub(A, temp, A);
        else 
            nmod_poly_sub(A, A, temp);

        nmod_poly_fit_length(A, k + a3->length);
        for (i = A->length; i < k + a3->length; i++) 
            A->coeffs[i] = 0L;
        nmod_poly_attach_shift(a4, A, k);
        a4->alloc = FLINT_MAX(a4->length, a3->length);
        nmod_poly_add(a4, a4, a3);
        A->length = FLINT_MAX(k + a3->length, A->length);
        _nmod_poly_normalise(A);

if (flag) {

        nmod_poly_swap(nmod_poly_mat_entry(S,0,0), nmod_poly_mat_entry(S,1,0));
        nmod_poly_swap(nmod_poly_mat_entry(S,0,1), nmod_poly_mat_entry(S,1,1));
        nmod_poly_mul(temp, nmod_poly_mat_entry(S,1,0), q);
        nmod_poly_add(nmod_poly_mat_entry(S,0,0), nmod_poly_mat_entry(S,0,0), temp);
        nmod_poly_mul(temp, nmod_poly_mat_entry(S,1,1), q);
        nmod_poly_add(nmod_poly_mat_entry(S,0,1), nmod_poly_mat_entry(S,0,1), temp);

        nmod_poly_mat_mul(R, R, S);
}

        nmod_poly_mat_clear(S);
        nmod_poly_clear(temp);
        nmod_poly_clear(a3); 
        nmod_poly_clear(b3);
        nmod_poly_clear(q);
        nmod_poly_clear(d);
        nmod_poly_clear(a2);
        nmod_poly_clear(b2);

        return - (signR * signS);
    }
}

void nmod_poly_half_gcd_no_matrix(nmod_poly_t a_out, nmod_poly_t b_out, nmod_poly_t a, nmod_poly_t b)
{
    nmod_poly_mat_t mat;

    if (a->mod.n == 13 && (a->length == 824 && b->length == 823))
    {
        FLAG = 0;
    }

printf("Call #1--------------------------------------------------\n"); fflush(stdout);
    nmod_poly_mat_init(mat, 2, 2, a->mod.n);

    nmod_poly_half_gcd(mat, a_out, b_out, a, b, 0);

    nmod_poly_mat_clear(mat);

    FLAG = 0;
}

void nmod_poly_gcd_hgcd(nmod_poly_t res, nmod_poly_t f, const nmod_poly_t g)
{
    const nmod_t mod  = f->mod;
    const long CUTOFF = FLINT_BIT_COUNT(mod.n) <= 8 ? 
    NMOD_POLY_SMALL_GCD_CUTOFF : NMOD_POLY_GCD_CUTOFF;
    nmod_poly_t h, j, r;

printf("Call #0--------------------------------------------------\n"); fflush(stdout);

    if (f->length == 0)
    {
        if (g->length == 0) nmod_poly_zero(res);
        else nmod_poly_make_monic(res, g);
        return;
    }
    if (g->length == 0)
    {
        nmod_poly_make_monic(res, f);
        return;
    }

    nmod_poly_init_preinv(r, mod.n, mod.ninv);

    nmod_poly_rem(r, f, g);
    if (r->length == 0)
    {
        nmod_poly_make_monic(res, g);
        nmod_poly_clear(r);
        return;
    }

    nmod_poly_init_preinv(j, mod.n, mod.ninv);
    nmod_poly_init_preinv(h, mod.n, mod.ninv);

    nmod_poly_half_gcd_no_matrix(h, j, g, r);

    while (j->length != 0)
    {
        nmod_poly_rem(r, h, j);

        if (r->length == 0)
        {
            nmod_poly_make_monic(res, j);
            goto gcd_exit;
        }

        if (j->length < CUTOFF)
        {
            nmod_poly_gcd_euclidean(res, j, r);
            goto gcd_exit;
        }

        nmod_poly_half_gcd_no_matrix(h, j, j, r);
    }

    nmod_poly_make_monic(res, h);

  gcd_exit: 

    nmod_poly_clear(j);
    nmod_poly_clear(h);
    nmod_poly_clear(r);
}
