/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

void __gmpn_mul_basecase(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);

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
Generic version of mul_n:

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

The compiler may refuse to unroll the nested loops, so we generate the code.

Schema for a general n x m multiply (here n = 7, m = 4):

    v0   u0 u1 u2 u3 u4 u5 u6 .
    v1      u0 u1 u2 u3 u4 u5 u6 .
    v2         u0 u1 u2 u3 u4 u5 u6 .
    v3            u0 u1 u2 u3 u4 u5 u6 .

def mul1v(n):
    print("mp_limb_t flint_mpn_mul_%ix1v(mp_ptr res, mp_srcptr u, mp_limb_t v0)" % n)
    print("{")
    print("    mp_limb_t a;")
    print("    NN_MUL_1X1(a, res[0], u[0], v0);")
    for i in range(1, n-1):
        print("    NN_ADDMUL_S2_A2_1X1(a, res[%i], 0, a, u[%i], v0);" % (i, i))
    print("    NN_ADDMUL_S2_A2_1X1(a, res[%i], 0, a, u[%i], v0);" % (n - 1, n - 1))
    print("    return a;")
    print("}")

def mulnm(n, m):
    if m == 1:
        print("mp_limb_t flint_mpn_mul_%ix1(mp_ptr res, mp_srcptr u, mp_srcptr v)" % n)
        print("{")
        print("    mp_limb_t a, v0 = v[0];")
        print("    NN_MUL_1X1(a, res[0], u[0], v0);")
        for i in range(1, n-1):
            print("    NN_ADDMUL_S2_A2_1X1(a, res[%i], 0, a, u[%i], v0);" % (i, i))
        print("    NN_ADDMUL_S2_A2_1X1(res[%i], res[%i], 0, a, u[%i], v0);" % (n, n - 1, n - 1))
        print("    return res[%i];" % n)
        print("}")
    elif m == 2:
        print("mp_limb_t flint_mpn_mul_%ix%i(mp_ptr res, mp_srcptr u, mp_srcptr v)" % (n, m))
        print("{")
        print("    mp_limb_t b, a;")
        print("    mp_limb_t w[2];")
        print("    w[0] = v[0];")
        print("    w[1] = v[1];")
        print("    NN_MUL_1X1(a, res[0], u[0], w[0]);")
        print("    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, w, 2);")
        for i in range(2, m):
            print("    NN_DOTREV_S3_A3_1X1(b, a, res[%i], 0, b, a, u, w, %i);" % (i, i + 1))
        for i in range(m, n):
            print("    NN_DOTREV_S3_A3_1X1(b, a, res[%i], 0, b, a, u + %i, w, %i);" % (i, i - m + 1, m))
        for i in range(n, n+m-2):
            print("    NN_DOTREV_S3_A3_1X1(b, a, res[%i], 0, b, a, u + %i, w + %i, %i);" % (i, i - m + 1, i - n + 1, n + m - i - 1))
        print("    NN_ADDMUL_S2_A2_1X1(res[%i], res[%i], b, a, u[%i], w[%i]);" % (n + m - 1, n + m - 2, n - 1, m - 1))
        print("    return res[%i];" % (n + m - 1))
        print("}")
    else:
        print("mp_limb_t flint_mpn_mul_%ix%i(mp_ptr res, mp_srcptr u, mp_srcptr v)" % (n, m))
        print("{")
        print("    mp_limb_t b, a;")
        print("    NN_MUL_1X1(a, res[0], u[0], v[0]);")
        print("    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);")
        for i in range(2, m):
            print("    NN_DOTREV_S3_A3_1X1(b, a, res[%i], 0, b, a, u, v, %i);" % (i, i + 1))
        for i in range(m, n):
            print("    NN_DOTREV_S3_A3_1X1(b, a, res[%i], 0, b, a, u + %i, v, %i);" % (i, i - m + 1, m))
        for i in range(n, n+m-2):
            print("    NN_DOTREV_S3_A3_1X1(b, a, res[%i], 0, b, a, u + %i, v + %i, %i);" % (i, i - m + 1, i - n + 1, n + m - i - 1))
        print("    NN_ADDMUL_S2_A2_1X1(res[%i], res[%i], b, a, u[%i], v[%i]);" % (n + m - 1, n + m - 2, n - 1, m - 1))
        print("    return res[%i];" % (n + m - 1))
        print("}")

for n in range(2, 8):
    for m in range(1, n+1):
        if n >= 8 and m >= 5:
            print("void flint_mpn_mul_%ix%i(mp_ptr res, mp_srcptr u, mp_srcptr v)" % (n, m))
            print("{")
            print("    __gmpn_mul_basecase(res, u, %i, v, %i);" % (n, m))
            print("}")
        else:
            mulnm(n, m)
        print()
    print()

for n in range(8, 17):
    mulnm(n, 1)
    print()
    #mulnm(n, 2)
    #print()

*/


mp_limb_t flint_mpn_mul_1x1(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    NN_MUL_1X1(res[1], res[0], u[0], v[0]);
    return res[1];
}

mp_limb_t flint_mpn_mul_2x1(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t a, v0 = v[0];
    NN_MUL_1X1(a, res[0], u[0], v0);
    NN_ADDMUL_S2_A2_1X1(res[2], res[1], 0, a, u[1], v0);
    return res[2];
}

mp_limb_t flint_mpn_mul_2x2(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    mp_limb_t w[2];
    w[0] = v[0];
    w[1] = v[1];
    NN_MUL_1X1(a, res[0], u[0], w[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, w, 2);
    NN_ADDMUL_S2_A2_1X1(res[3], res[2], b, a, u[1], w[1]);
    return res[3];
}

mp_limb_t flint_mpn_mul_3x1(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t a, v0 = v[0];
    NN_MUL_1X1(a, res[0], u[0], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[1], 0, a, u[1], v0);
    NN_ADDMUL_S2_A2_1X1(res[3], res[2], 0, a, u[2], v0);
    return res[3];
}

mp_limb_t flint_mpn_mul_3x2(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    mp_limb_t w[2];
    w[0] = v[0];
    w[1] = v[1];
    NN_MUL_1X1(a, res[0], u[0], w[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, w, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u + 1, w, 2);
    NN_ADDMUL_S2_A2_1X1(res[4], res[3], b, a, u[2], w[1]);
    return res[4];
}

mp_limb_t flint_mpn_mul_3x3(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u + 1, v + 1, 2);
    NN_ADDMUL_S2_A2_1X1(res[5], res[4], b, a, u[2], v[2]);
    return res[5];
}


mp_limb_t flint_mpn_mul_4x1(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t a, v0 = v[0];
    NN_MUL_1X1(a, res[0], u[0], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[1], 0, a, u[1], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[2], 0, a, u[2], v0);
    NN_ADDMUL_S2_A2_1X1(res[4], res[3], 0, a, u[3], v0);
    return res[4];
}

mp_limb_t flint_mpn_mul_4x2(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    mp_limb_t w[2];
    w[0] = v[0];
    w[1] = v[1];
    NN_MUL_1X1(a, res[0], u[0], w[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, w, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u + 1, w, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u + 2, w, 2);
    NN_ADDMUL_S2_A2_1X1(res[5], res[4], b, a, u[3], w[1]);
    return res[5];
}

mp_limb_t flint_mpn_mul_4x3(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u + 1, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u + 2, v + 1, 2);
    NN_ADDMUL_S2_A2_1X1(res[6], res[5], b, a, u[3], v[2]);
    return res[6];
}

mp_limb_t flint_mpn_mul_4x4(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u + 1, v + 1, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u + 2, v + 2, 2);
    NN_ADDMUL_S2_A2_1X1(res[7], res[6], b, a, u[3], v[3]);
    return res[7];
}


mp_limb_t flint_mpn_mul_5x1(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t a, v0 = v[0];
    NN_MUL_1X1(a, res[0], u[0], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[1], 0, a, u[1], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[2], 0, a, u[2], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[3], 0, a, u[3], v0);
    NN_ADDMUL_S2_A2_1X1(res[5], res[4], 0, a, u[4], v0);
    return res[5];
}

mp_limb_t flint_mpn_mul_5x2(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    mp_limb_t w[2];
    w[0] = v[0];
    w[1] = v[1];
    NN_MUL_1X1(a, res[0], u[0], w[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, w, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u + 1, w, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u + 2, w, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u + 3, w, 2);
    NN_ADDMUL_S2_A2_1X1(res[6], res[5], b, a, u[4], w[1]);
    return res[6];
}

mp_limb_t flint_mpn_mul_5x3(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u + 1, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u + 2, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u + 3, v + 1, 2);
    NN_ADDMUL_S2_A2_1X1(res[7], res[6], b, a, u[4], v[2]);
    return res[7];
}

mp_limb_t flint_mpn_mul_5x4(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u + 1, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u + 2, v + 1, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u + 3, v + 2, 2);
    NN_ADDMUL_S2_A2_1X1(res[8], res[7], b, a, u[4], v[3]);
    return res[8];
}

mp_limb_t flint_mpn_mul_5x5(mp_ptr res, mp_srcptr u, mp_srcptr v)
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
    return res[9];
}


mp_limb_t flint_mpn_mul_6x1(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t a, v0 = v[0];
    NN_MUL_1X1(a, res[0], u[0], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[1], 0, a, u[1], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[2], 0, a, u[2], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[3], 0, a, u[3], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[4], 0, a, u[4], v0);
    NN_ADDMUL_S2_A2_1X1(res[6], res[5], 0, a, u[5], v0);
    return res[6];
}

mp_limb_t flint_mpn_mul_6x2(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    mp_limb_t w[2];
    w[0] = v[0];
    w[1] = v[1];
    NN_MUL_1X1(a, res[0], u[0], w[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, w, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u + 1, w, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u + 2, w, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u + 3, w, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u + 4, w, 2);
    NN_ADDMUL_S2_A2_1X1(res[7], res[6], b, a, u[5], w[1]);
    return res[7];
}

mp_limb_t flint_mpn_mul_6x3(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u + 1, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u + 2, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u + 3, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u + 4, v + 1, 2);
    NN_ADDMUL_S2_A2_1X1(res[8], res[7], b, a, u[5], v[2]);
    return res[8];
}

mp_limb_t flint_mpn_mul_6x4(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u + 1, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u + 2, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u + 3, v + 1, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[7], 0, b, a, u + 4, v + 2, 2);
    NN_ADDMUL_S2_A2_1X1(res[9], res[8], b, a, u[5], v[3]);
    return res[9];
}

mp_limb_t flint_mpn_mul_6x5(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u, v, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u + 1, v, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u + 2, v + 1, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[7], 0, b, a, u + 3, v + 2, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[8], 0, b, a, u + 4, v + 3, 2);
    NN_ADDMUL_S2_A2_1X1(res[10], res[9], b, a, u[5], v[4]);
    return res[10];
}

mp_limb_t flint_mpn_mul_6x6(mp_ptr res, mp_srcptr u, mp_srcptr v)
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
    return res[11];
}


mp_limb_t flint_mpn_mul_7x1(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t a, v0 = v[0];
    NN_MUL_1X1(a, res[0], u[0], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[1], 0, a, u[1], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[2], 0, a, u[2], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[3], 0, a, u[3], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[4], 0, a, u[4], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[5], 0, a, u[5], v0);
    NN_ADDMUL_S2_A2_1X1(res[7], res[6], 0, a, u[6], v0);
    return res[7];
}

mp_limb_t flint_mpn_mul_7x2(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    mp_limb_t w[2];
    w[0] = v[0];
    w[1] = v[1];
    NN_MUL_1X1(a, res[0], u[0], w[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, w, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u + 1, w, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u + 2, w, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u + 3, w, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u + 4, w, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u + 5, w, 2);
    NN_ADDMUL_S2_A2_1X1(res[8], res[7], b, a, u[6], w[1]);
    return res[8];
}

mp_limb_t flint_mpn_mul_7x3(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u + 1, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u + 2, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u + 3, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u + 4, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[7], 0, b, a, u + 5, v + 1, 2);
    NN_ADDMUL_S2_A2_1X1(res[9], res[8], b, a, u[6], v[2]);
    return res[9];
}

mp_limb_t flint_mpn_mul_7x4(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u + 1, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u + 2, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u + 3, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[7], 0, b, a, u + 4, v + 1, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[8], 0, b, a, u + 5, v + 2, 2);
    NN_ADDMUL_S2_A2_1X1(res[10], res[9], b, a, u[6], v[3]);
    return res[10];
}

mp_limb_t flint_mpn_mul_7x5(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u, v, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u + 1, v, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u + 2, v, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[7], 0, b, a, u + 3, v + 1, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[8], 0, b, a, u + 4, v + 2, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[9], 0, b, a, u + 5, v + 3, 2);
    NN_ADDMUL_S2_A2_1X1(res[11], res[10], b, a, u[6], v[4]);
    return res[11];
}

mp_limb_t flint_mpn_mul_7x6(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_DOTREV_S3_A3_1X1(b, a, res[2], 0, b, a, u, v, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[3], 0, b, a, u, v, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[4], 0, b, a, u, v, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[5], 0, b, a, u, v, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[6], 0, b, a, u + 1, v, 6);
    NN_DOTREV_S3_A3_1X1(b, a, res[7], 0, b, a, u + 2, v + 1, 5);
    NN_DOTREV_S3_A3_1X1(b, a, res[8], 0, b, a, u + 3, v + 2, 4);
    NN_DOTREV_S3_A3_1X1(b, a, res[9], 0, b, a, u + 4, v + 3, 3);
    NN_DOTREV_S3_A3_1X1(b, a, res[10], 0, b, a, u + 5, v + 4, 2);
    NN_ADDMUL_S2_A2_1X1(res[12], res[11], b, a, u[6], v[5]);
    return res[12];
}

mp_limb_t flint_mpn_mul_7x7(mp_ptr res, mp_srcptr u, mp_srcptr v)
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
    return res[13];
}

mp_limb_t flint_mpn_mul_8x1(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t a, v0 = v[0];
    NN_MUL_1X1(a, res[0], u[0], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[1], 0, a, u[1], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[2], 0, a, u[2], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[3], 0, a, u[3], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[4], 0, a, u[4], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[5], 0, a, u[5], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[6], 0, a, u[6], v0);
    NN_ADDMUL_S2_A2_1X1(res[8], res[7], 0, a, u[7], v0);
    return res[8];
}

mp_limb_t flint_mpn_mul_9x1(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t a, v0 = v[0];
    NN_MUL_1X1(a, res[0], u[0], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[1], 0, a, u[1], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[2], 0, a, u[2], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[3], 0, a, u[3], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[4], 0, a, u[4], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[5], 0, a, u[5], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[6], 0, a, u[6], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[7], 0, a, u[7], v0);
    NN_ADDMUL_S2_A2_1X1(res[9], res[8], 0, a, u[8], v0);
    return res[9];
}

mp_limb_t flint_mpn_mul_10x1(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t a, v0 = v[0];
    NN_MUL_1X1(a, res[0], u[0], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[1], 0, a, u[1], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[2], 0, a, u[2], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[3], 0, a, u[3], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[4], 0, a, u[4], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[5], 0, a, u[5], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[6], 0, a, u[6], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[7], 0, a, u[7], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[8], 0, a, u[8], v0);
    NN_ADDMUL_S2_A2_1X1(res[10], res[9], 0, a, u[9], v0);
    return res[10];
}

mp_limb_t flint_mpn_mul_11x1(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t a, v0 = v[0];
    NN_MUL_1X1(a, res[0], u[0], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[1], 0, a, u[1], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[2], 0, a, u[2], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[3], 0, a, u[3], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[4], 0, a, u[4], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[5], 0, a, u[5], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[6], 0, a, u[6], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[7], 0, a, u[7], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[8], 0, a, u[8], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[9], 0, a, u[9], v0);
    NN_ADDMUL_S2_A2_1X1(res[11], res[10], 0, a, u[10], v0);
    return res[11];
}

mp_limb_t flint_mpn_mul_12x1(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t a, v0 = v[0];
    NN_MUL_1X1(a, res[0], u[0], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[1], 0, a, u[1], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[2], 0, a, u[2], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[3], 0, a, u[3], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[4], 0, a, u[4], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[5], 0, a, u[5], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[6], 0, a, u[6], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[7], 0, a, u[7], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[8], 0, a, u[8], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[9], 0, a, u[9], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[10], 0, a, u[10], v0);
    NN_ADDMUL_S2_A2_1X1(res[12], res[11], 0, a, u[11], v0);
    return res[12];
}

mp_limb_t flint_mpn_mul_13x1(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t a, v0 = v[0];
    NN_MUL_1X1(a, res[0], u[0], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[1], 0, a, u[1], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[2], 0, a, u[2], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[3], 0, a, u[3], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[4], 0, a, u[4], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[5], 0, a, u[5], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[6], 0, a, u[6], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[7], 0, a, u[7], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[8], 0, a, u[8], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[9], 0, a, u[9], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[10], 0, a, u[10], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[11], 0, a, u[11], v0);
    NN_ADDMUL_S2_A2_1X1(res[13], res[12], 0, a, u[12], v0);
    return res[13];
}

mp_limb_t flint_mpn_mul_14x1(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t a, v0 = v[0];
    NN_MUL_1X1(a, res[0], u[0], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[1], 0, a, u[1], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[2], 0, a, u[2], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[3], 0, a, u[3], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[4], 0, a, u[4], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[5], 0, a, u[5], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[6], 0, a, u[6], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[7], 0, a, u[7], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[8], 0, a, u[8], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[9], 0, a, u[9], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[10], 0, a, u[10], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[11], 0, a, u[11], v0);
    NN_ADDMUL_S2_A2_1X1(a, res[12], 0, a, u[12], v0);
    NN_ADDMUL_S2_A2_1X1(res[14], res[13], 0, a, u[13], v0);
    return res[14];
}

typedef mp_limb_t (*flint_mpn_mul_func_t)(mp_ptr, mp_srcptr, mp_srcptr);

const flint_mpn_mul_func_t flint_mpn_mul_tab[8][8] = {
    { NULL, },
    { NULL, flint_mpn_mul_1x1, },
    { NULL, flint_mpn_mul_2x1, flint_mpn_mul_2x2, },
    { NULL, flint_mpn_mul_3x1, flint_mpn_mul_3x2, flint_mpn_mul_3x3, },
    { NULL, flint_mpn_mul_4x1, flint_mpn_mul_4x2, flint_mpn_mul_4x3, flint_mpn_mul_4x4, },
    { NULL, flint_mpn_mul_5x1, flint_mpn_mul_5x2, flint_mpn_mul_5x3, flint_mpn_mul_5x4, flint_mpn_mul_5x5, },
    { NULL, flint_mpn_mul_6x1, flint_mpn_mul_6x2, flint_mpn_mul_6x3, flint_mpn_mul_6x4, flint_mpn_mul_6x5, flint_mpn_mul_6x6, },
    { NULL, flint_mpn_mul_7x1, flint_mpn_mul_7x2, flint_mpn_mul_7x3, flint_mpn_mul_7x4, flint_mpn_mul_7x5, flint_mpn_mul_7x6, flint_mpn_mul_7x7, },
};

const flint_mpn_mul_func_t flint_mpn_mul_n_tab[8] = {
    NULL,
    flint_mpn_mul_1x1,
    flint_mpn_mul_2x2,
    flint_mpn_mul_3x3,
    flint_mpn_mul_4x4,
    flint_mpn_mul_5x5,
    flint_mpn_mul_6x6,
    flint_mpn_mul_7x7,
};

const flint_mpn_mul_func_t flint_mpn_mul_1_tab[15] = {
    NULL,
    flint_mpn_mul_1x1,
    flint_mpn_mul_2x1,
    flint_mpn_mul_3x1,
    flint_mpn_mul_4x1,
    flint_mpn_mul_5x1,
    flint_mpn_mul_6x1,
    flint_mpn_mul_7x1,
    flint_mpn_mul_8x1,
    flint_mpn_mul_9x1,
    flint_mpn_mul_10x1,
    flint_mpn_mul_11x1,
    flint_mpn_mul_12x1,
    flint_mpn_mul_13x1,
    flint_mpn_mul_14x1,
};

FLINT_FORCE_INLINE
mp_limb_t flint_mpn_mul_1x1_inline(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    NN_MUL_1X1(res[1], res[0], u[0], v[0]);
    return res[1];
}

FLINT_FORCE_INLINE
mp_limb_t flint_mpn_mul_2x1_inline(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t a, v0 = v[0];
    NN_MUL_1X1(a, res[0], u[0], v0);
    NN_ADDMUL_S2_A2_1X1(res[2], res[1], 0, a, u[1], v0);
    return res[2];
}

FLINT_FORCE_INLINE
mp_limb_t flint_mpn_mul_2x2_inline(mp_ptr res, mp_srcptr u, mp_srcptr v)
{
    mp_limb_t b, a;
    NN_MUL_1X1(a, res[0], u[0], v[0]);
    NN_DOTREV_S3_A3_1X1(b, a, res[1], 0, 0, a, u, v, 2);
    NN_ADDMUL_S2_A2_1X1(res[3], res[2], b, a, u[1], v[1]);
    return res[3];
}

#define MUL_STATS 0

#if MUL_STATS

slong mulcount[17][17] = { { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                           { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                           { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                           { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                           { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                           { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                           { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                           { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                           { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                           { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                           { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                           { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                           { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                           { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                           { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                           { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                           { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, };

slong mulcounter = 0;

#define T mulncount
#endif

void
flint_mpn_mul_n(mp_ptr z, mp_srcptr x, mp_srcptr y, mp_size_t n)
{
    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(z != x);
    FLINT_ASSERT(z != y);

#if MUL_STATS
    mulcounter++;
    mulcount[FLINT_MIN(n, 16)][FLINT_MIN(n,16)]++;
    if (mulcounter % 10000 == 0)
    {
        flint_printf("count %wd\n", mulcounter);
        int i, j;
        for (i = 1; i <= 16; i++)
        {
            for (j = 1; j <= i; j++)
            {
                flint_printf("%7wd ", mulcount[i][j]);
            }
            flint_printf("\n");
        }

    }
#endif

    /* Experimentally inline n = 2. Whether this is beneficial depends on the
       distribution of inputs. We do not waste a branch on checking for n = 1,
       following the assumption that many callers already handle single-limb
       input specially. */
    if (n == 2)
        flint_mpn_mul_2x2_inline(z, x, y);
    else if (n <= 7)
        flint_mpn_mul_n_tab[n](z, x, y);
    else if (n <= 14)
        __gmpn_mul_basecase(z, x, n, y, n);
    else if (n < FLINT_MPN_MUL_THRESHOLD)
        mpn_mul_n(z, x, y, n);
    else
        flint_mpn_mul_large(z, x, n, y, n);
}

mp_limb_t
flint_mpn_mul(mp_ptr z, mp_srcptr x, mp_size_t xn, mp_srcptr y, mp_size_t yn)
{
    FLINT_ASSERT(xn >= yn);
    FLINT_ASSERT(yn >= 1);
    FLINT_ASSERT(z != x);
    FLINT_ASSERT(z != y);

#if MUL_STATS
    mulcounter++;
    mulcount[FLINT_MIN(xn, 16)][FLINT_MIN(yn, 16)]++;
    if (mulcounter % 10000 == 0)
    {
        flint_printf("count %wd\n", mulcounter);
        int i, j;
        for (i = 1; i <= 16; i++)
        {
            for (j = 1; j <= i; j++)
            {
                flint_printf("%7wd ", mulcount[i][j]);
            }
            flint_printf("\n");
        }
    }
#endif

    /* Experimentally inline n = 2. Whether this is beneficial depends on the
       distribution of inputs. We do not waste a branch on checking for n = 1,
       following the assumption that many callers already handle single-limb
       input specially. */
    if (xn == 2)
    {
        if (yn == 1)
            return flint_mpn_mul_2x1_inline(z, x, y);
        else
            return flint_mpn_mul_2x2_inline(z, x, y);
    }
    else if (xn <= 7)
    {
        return flint_mpn_mul_tab[xn][yn](z, x, y);
    }
    else if (yn == 1)
    {
        if (xn <= 12)
            return flint_mpn_mul_1_tab[xn](z, x, y);
        else
            return (z[xn + yn - 1] = mpn_mul_1(z, x, xn, y[0]));
    }
    else if (xn <= 14)
    {
        __gmpn_mul_basecase(z, x, xn, y, yn);
        return z[xn + yn - 1];
    }
    else if (yn < FLINT_MPN_MUL_THRESHOLD)
        return mpn_mul(z, x, xn, y, yn);
    else
        return flint_mpn_mul_large(z, x, xn, y, yn);
}

