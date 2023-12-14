/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/* Some helper macros that eventually should be public. */

/* {s0,s1,s2} = u[0]v[n-1] + u[1]v[n-2] + ... */
/* Assumes n >= 2 */
#define NN_DOTREV_S3_1X1(s2, s1, s0, u, v, n) \
    do { \
        mp_limb_t __dt0, __dt1, __ds0, __ds1, __ds2; \
        slong __i; \
        FLINT_ASSERT((n) >= 2); \
        umul_ppmm(__ds1, __ds0, (u)[0], (v)[(n) - 1]); \
        umul_ppmm(__dt1, __dt0, (u)[1], (v)[(n) - 2]); \
        add_sssaaaaaa(__ds2, __ds1, __ds0, 0, __ds1, __ds0, 0, __dt1, __dt0); \
        for (__i = 2; __i < (n); __i++) \
        { \
            umul_ppmm(__dt1, __dt0, (u)[i], (v)[(n) - 1 - __i]); \
            add_sssaaaaaa(__ds2, __ds1, __ds0, __ds2, __ds1, __ds0, 0, __dt1, __dt0); \
        } \
        (s0) = __ds0; (s1) = __ds1; (s2) = __ds2; \
    } while (0) \

/* {r0,r1,r2} = {s0,s1,s2} + u[0]v[n-1] + u[1]v[n-2] + ... */
/* Assumes n >= 1. May have s2 != 0, but the final sum is assumed to fit in 3 limbs. */
#define NN_DOTREV_S3_A3_1X1(r2, r1, r0, s2, s1, s0, u, v, n) \
    do { \
        mp_limb_t __dt0, __dt1, __ds0, __ds1, __ds2; \
        slong __i; \
        FLINT_ASSERT((n) >= 1); \
        __ds0 = (s0); __ds1 = (s1); __ds2 = (s2); \
        for (__i = 0; __i < (n); __i++) \
        { \
            umul_ppmm(__dt1, __dt0, (u)[__i], (v)[(n) - 1 - __i]); \
            add_sssaaaaaa(__ds2, __ds1, __ds0, __ds2, __ds1, __ds0, 0, __dt1, __dt0); \
        } \
        (r0) = __ds0; (r1) = __ds1; (r2) = __ds2; \
    } while (0) \

#define NN_MUL_1X1 umul_ppmm

/* {r0,r1} = {s0,s1} + x * y, with no carry-out. */
#define NN_ADDMUL_S2_A2_1X1(r1, r0, s1, s0, x, y) \
    do { \
        mp_limb_t __dt0, __dt1; \
        umul_ppmm(__dt1, __dt0, (x), (y)); \
        add_ssaaaa(r1, r0, s1, s0, __dt1, __dt0); \
    } while (0); \


/*
Generic version:

void flint_mpn_mul_n_basecase(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
{
    mp_limb_t b, a;
    slong i;

    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);

    for (i = 2; i < n; i++)
        NN_DOTREV_S3_A3_1X1(b, a, res[i], 0, b, a, u, v, i + 1);

    for (i = n; i < 2 * n - 2; i++)
        NN_DOTREV_S3_A3_1X1(b, a, res[i], 0, b, a, u + i - n + 1, v + i - n + 1, 2 * n - i - 1);

    NN_ADDMUL_S2_A2_1X1(res[2 * n - 1], res[2 * n - 2], b, a, u[n - 1], v[n - 1]);
}

The compiler may refuse to unroll the nested loops, so we generate the code:

def mul(n):
    print("void flint_mpn_mul_%i(mp_ptr res, mp_srcptr u, mp_srcptr v)" % n)
    print("{")
    print("    mp_limb_t b, a;")
    print("    NN_MUL_1X1(a, res[0], u[0], v[0]);")
    print("    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);")
    for i in range(2, n):
        print("    NN_DOTREV_S3_A3_1X1(b, a, res[%i], 0, b, a, u, v, %i);" % (i, i + 1))
    for i in range(n, 2 * n - 2):
        print("    NN_DOTREV_S3_A3_1X1(b, a, res[%i], 0, b, a, u + %i, v + %i, %i);" % (i, i - n + 1, i - n + 1, 2 * n - i - 1))
    print("    NN_ADDMUL_S2_A2_1X1(res[%i], res[%i], b, a, u[%i], v[%i]);" % (2 * n - 1, 2 * n - 2, n - 1, n - 1))
    print("}")

*/

FLINT_FORCE_INLINE
void flint_mpn_mul_1(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    NN_MUL_1X1(res[1], res[0], u[0], v[0]);
}

FLINT_FORCE_INLINE
void flint_mpn_mul_2(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_ADDMUL_S2_A2_1X1(res[3], res[2], b, a, u[1], v[1]);
}

void flint_mpn_mul_3(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u + 1, v + 1, 2);
    NN_ADDMUL_S2_A2_1X1(res[5], res[4], b, a, u[2], v[2]);
}

void flint_mpn_mul_4(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u + 1, v + 1, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u + 2, v + 2, 2);
    NN_ADDMUL_S2_A2_1X1(res[7], res[6], b, a, u[3], v[3]);
}

void flint_mpn_mul_5(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u, v, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u + 1, v + 1, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u + 2, v + 2, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[7], 0, b, a, u + 3, v + 3, 2);
    NN_ADDMUL_S2_A2_1X1(res[9], res[8], b, a, u[4], v[4]);
}

void flint_mpn_mul_6(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u, v, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u, v, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u + 1, v + 1, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[7], 0, b, a, u + 2, v + 2, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[8], 0, b, a, u + 3, v + 3, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[9], 0, b, a, u + 4, v + 4, 2);
    NN_ADDMUL_S2_A2_1X1(res[11], res[10], b, a, u[5], v[5]);
}

void flint_mpn_mul_7(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u, v, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u, v, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u, v, 7);
    NN_DOTREV_S3_A3_1X1(b, a, res[7], 0, b, a, u + 1, v + 1, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[8], 0, b, a, u + 2, v + 2, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[9], 0, b, a, u + 3, v + 3, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[10], 0, b, a, u + 4, v + 4, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[11], 0, b, a, u + 5, v + 5, 2);
    NN_ADDMUL_S2_A2_1X1(res[13], res[12], b, a, u[6], v[6]);
}

void flint_mpn_mul_8(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u, v, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u, v, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u, v, 7);
    NN_DOTREV_S3_A3_1X1(b, a, res[7], 0, b, a, u, v, 8);
    NN_DOTREV_S3_A3_1X1(b, a, res[8], 0, b, a, u + 1, v + 1, 7);
    NN_DOTREV_S3_A3_1X1(b, a, res[9], 0, b, a, u + 2, v + 2, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[10], 0, b, a, u + 3, v + 3, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[11], 0, b, a, u + 4, v + 4, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[12], 0, b, a, u + 5, v + 5, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[13], 0, b, a, u + 6, v + 6, 2);
    NN_ADDMUL_S2_A2_1X1(res[15], res[14], b, a, u[7], v[7]);
}

void flint_mpn_mul_9(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u, v, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u, v, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u, v, 7);
    NN_DOTREV_S3_A3_1X1(b, a, res[7], 0, b, a, u, v, 8);
    NN_DOTREV_S3_A3_1X1(b, a, res[8], 0, b, a, u, v, 9);
    NN_DOTREV_S3_A3_1X1(b, a, res[9], 0, b, a, u + 1, v + 1, 8);
    NN_DOTREV_S3_A3_1X1(b, a, res[10], 0, b, a, u + 2, v + 2, 7);
    NN_DOTREV_S3_A3_1X1(b, a, res[11], 0, b, a, u + 3, v + 3, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[12], 0, b, a, u + 4, v + 4, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[13], 0, b, a, u + 5, v + 5, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[14], 0, b, a, u + 6, v + 6, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[15], 0, b, a, u + 7, v + 7, 2);
    NN_ADDMUL_S2_A2_1X1(res[17], res[16], b, a, u[8], v[8]);
}

void flint_mpn_mul_10(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u, v, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u, v, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u, v, 7);
    NN_DOTREV_S3_A3_1X1(b, a, res[7], 0, b, a, u, v, 8);
    NN_DOTREV_S3_A3_1X1(b, a, res[8], 0, b, a, u, v, 9);
    NN_DOTREV_S3_A3_1X1(b, a, res[9], 0, b, a, u, v, 10);
    NN_DOTREV_S3_A3_1X1(b, a, res[10], 0, b, a, u + 1, v + 1, 9);
    NN_DOTREV_S3_A3_1X1(b, a, res[11], 0, b, a, u + 2, v + 2, 8);
    NN_DOTREV_S3_A3_1X1(b, a, res[12], 0, b, a, u + 3, v + 3, 7);
    NN_DOTREV_S3_A3_1X1(b, a, res[13], 0, b, a, u + 4, v + 4, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[14], 0, b, a, u + 5, v + 5, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[15], 0, b, a, u + 6, v + 6, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[16], 0, b, a, u + 7, v + 7, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[17], 0, b, a, u + 8, v + 8, 2);
    NN_ADDMUL_S2_A2_1X1(res[19], res[18], b, a, u[9], v[9]);
}

void flint_mpn_mul_11(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u, v, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u, v, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u, v, 7);
    NN_DOTREV_S3_A3_1X1(b, a, res[7], 0, b, a, u, v, 8);
    NN_DOTREV_S3_A3_1X1(b, a, res[8], 0, b, a, u, v, 9);
    NN_DOTREV_S3_A3_1X1(b, a, res[9], 0, b, a, u, v, 10);
    NN_DOTREV_S3_A3_1X1(b, a, res[10], 0, b, a, u, v, 11);
    NN_DOTREV_S3_A3_1X1(b, a, res[11], 0, b, a, u + 1, v + 1, 10);
    NN_DOTREV_S3_A3_1X1(b, a, res[12], 0, b, a, u + 2, v + 2, 9);
    NN_DOTREV_S3_A3_1X1(b, a, res[13], 0, b, a, u + 3, v + 3, 8);
    NN_DOTREV_S3_A3_1X1(b, a, res[14], 0, b, a, u + 4, v + 4, 7);
    NN_DOTREV_S3_A3_1X1(b, a, res[15], 0, b, a, u + 5, v + 5, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[16], 0, b, a, u + 6, v + 6, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[17], 0, b, a, u + 7, v + 7, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[18], 0, b, a, u + 8, v + 8, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[19], 0, b, a, u + 9, v + 9, 2);
    NN_ADDMUL_S2_A2_1X1(res[21], res[20], b, a, u[10], v[10]);
}

void flint_mpn_mul_12(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u, v, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u, v, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u, v, 7);
    NN_DOTREV_S3_A3_1X1(b, a, res[7], 0, b, a, u, v, 8);
    NN_DOTREV_S3_A3_1X1(b, a, res[8], 0, b, a, u, v, 9);
    NN_DOTREV_S3_A3_1X1(b, a, res[9], 0, b, a, u, v, 10);
    NN_DOTREV_S3_A3_1X1(b, a, res[10], 0, b, a, u, v, 11);
    NN_DOTREV_S3_A3_1X1(b, a, res[11], 0, b, a, u, v, 12);
    NN_DOTREV_S3_A3_1X1(b, a, res[12], 0, b, a, u + 1, v + 1, 11);
    NN_DOTREV_S3_A3_1X1(b, a, res[13], 0, b, a, u + 2, v + 2, 10);
    NN_DOTREV_S3_A3_1X1(b, a, res[14], 0, b, a, u + 3, v + 3, 9);
    NN_DOTREV_S3_A3_1X1(b, a, res[15], 0, b, a, u + 4, v + 4, 8);
    NN_DOTREV_S3_A3_1X1(b, a, res[16], 0, b, a, u + 5, v + 5, 7);
    NN_DOTREV_S3_A3_1X1(b, a, res[17], 0, b, a, u + 6, v + 6, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[18], 0, b, a, u + 7, v + 7, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[19], 0, b, a, u + 8, v + 8, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[20], 0, b, a, u + 9, v + 9, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[21], 0, b, a, u + 10, v + 10, 2);
    NN_ADDMUL_S2_A2_1X1(res[23], res[22], b, a, u[11], v[11]);
}

void
flint_mpn_mul_n(mp_ptr z, mp_srcptr x, mp_srcptr y, mp_size_t n)
{
    FLINT_ASSERT(n >= 1);

    switch (n)
    {
        case 0: FLINT_UNREACHABLE;
        case 1: flint_mpn_mul_1(z, x, y); break;
        case 2: flint_mpn_mul_2(z, x, y); break;
        case 3: flint_mpn_mul_3(z, x, y); break;
        case 4: flint_mpn_mul_4(z, x, y); break;
        case 5: flint_mpn_mul_5(z, x, y); break;
        case 6: flint_mpn_mul_6(z, x, y); break;
        case 7: flint_mpn_mul_7(z, x, y); break;
        case 8: flint_mpn_mul_8(z, x, y); break;
        case 9: flint_mpn_mul_9(z, x, y); break;
        case 10: flint_mpn_mul_10(z, x, y); break;
        case 11: flint_mpn_mul_11(z, x, y); break;
        case 12: flint_mpn_mul_12(z, x, y); break;
        default:
            if (n < FLINT_MPN_MUL_THRESHOLD)
                mpn_mul_n(z, x, y, n);
            else
                flint_mpn_mul_large(z, x, n, y, n);
    }
}

