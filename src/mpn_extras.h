/*
    Copyright (C) 2010 Fredrik Johansson

    2x2 mul code taken from MPFR 2.3.0
    (Copyright (C) 1991-2007 Free Software Foundation, Inc.)

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MPN_EXTRAS_H
#define MPN_EXTRAS_H

#ifdef MPN_EXTRAS_INLINES_C
#define MPN_EXTRAS_INLINE
#else
#define MPN_EXTRAS_INLINE static inline
#endif

#include <gmp.h>
#include "longlong.h"

#ifdef __cplusplus
extern "C" {
#endif

/* mpn macros ****************************************************************/

FLINT_FORCE_INLINE
void flint_mpn_zero(mp_ptr xp, mp_size_t n)
{
    mp_size_t ix;
    for (ix = 0; ix < n; ix++)
        xp[ix] = UWORD(0);
}

FLINT_FORCE_INLINE
void flint_mpn_copyi(mp_ptr xp, mp_srcptr yp, mp_size_t n)
{
    mp_size_t ix;
    for (ix = 0; ix < n; ix++)
        xp[ix] = yp[ix];
}

FLINT_FORCE_INLINE
void flint_mpn_copyd(mp_ptr xp, mp_srcptr yp, mp_size_t n)
{
    mp_size_t ix;
    for (ix = n - 1; ix >= 0; ix--)
        xp[ix] = yp[ix];
}

FLINT_FORCE_INLINE
void flint_mpn_store(mp_ptr xp, mp_size_t n, mp_limb_t y)
{
    mp_size_t ix;
    for (ix = 0; ix < n; ix++)
        xp[ix] = y;
}

MPN_EXTRAS_INLINE
int flint_mpn_zero_p(mp_srcptr x, mp_size_t xsize)
{
    slong i;
    for (i = 0; i < xsize; i++)
    {
        if (x[i])
            return 0;
    }
    return 1;
}

FLINT_FORCE_INLINE
int flint_mpn_equal_p(mp_srcptr x, mp_srcptr y, mp_size_t xsize)
{
    slong i;
    for (i = 0; i < xsize; i++)
    {
        if (x[i] != y[i])
            return 0;
    }
    return 1;
}

FLINT_FORCE_INLINE void
flint_mpn_negmod_n(mp_ptr res, mp_srcptr x, mp_srcptr m, mp_size_t n)
{
    if (flint_mpn_zero_p(x, n))
        flint_mpn_zero(res, n);
    else
        mpn_sub_n(res, m, x, n);
}

FLINT_FORCE_INLINE void
flint_mpn_addmod_n(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_srcptr m, mp_size_t n)
{
    mp_limb_t cy;
    cy = mpn_add_n(res, x, y, n);
    if (cy || mpn_cmp(res, m, n) >= 0)
        mpn_sub_n(res, res, m, n);
}

FLINT_FORCE_INLINE void
flint_mpn_submod_n(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_srcptr m, mp_size_t n)
{
    int cmp = (mpn_cmp(x, y, n) < 0);
    mpn_sub_n(res, x, y, n);
    if (cmp)
        mpn_add_n(res, res, m, n);
}

/* assumes yn <= n and y < m */
FLINT_FORCE_INLINE void
flint_mpn_addmod_n_m(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_size_t yn, mp_srcptr m, mp_size_t n)
{
    mp_limb_t cy;
    cy = mpn_add(res, x, n, y, yn);
    if (cy || mpn_cmp(res, m, n) >= 0)
        mpn_sub_n(res, res, m, n);
}

/* assumes yn <= n and y < m */
FLINT_FORCE_INLINE void
flint_mpn_submod_n_m(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_size_t yn, mp_srcptr m, mp_size_t n)
{
    int cmp = (flint_mpn_zero_p(x + yn, n - yn) && mpn_cmp(x, y, yn) < 0);
    mpn_sub(res, x, n, y, yn);
    if (cmp)
        mpn_add_n(res, res, m, n);
}

FLINT_FORCE_INLINE void
flint_mpn_negmod_2(mp_ptr res, mp_srcptr x, mp_srcptr m)
{
    if (x[0] == 0 && x[1] == 0)
        res[1] = res[0] = 0;
    else
        sub_ddmmss(res[1], res[0], m[1], m[0], x[1], x[0]);
}

FLINT_FORCE_INLINE void
flint_mpn_addmod_2(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_srcptr m)
{
    mp_limb_t cy;
    mp_limb_t m1 = m[1], m0 = m[0];
    add_sssaaaaaa(cy, res[1], res[0], 0, x[1], x[0], 0, y[1], y[0]);
    if (cy || (res[1] > m1 || (res[1] == m1 && res[0] >= m0)))
        sub_ddmmss(res[1], res[0], res[1], res[0], m1, m0);
}

/* assumes msb of m is zero */
FLINT_FORCE_INLINE void
_flint_mpn_addmod_2(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_srcptr m)
{
    mp_limb_t m1 = m[1], m0 = m[0];
    add_ssaaaa(res[1], res[0], x[1], x[0], y[1], y[0]);
    if (res[1] > m1 || (res[1] == m1 && res[0] >= m0))
        sub_ddmmss(res[1], res[0], res[1], res[0], m1, m0);
}

FLINT_FORCE_INLINE void
flint_mpn_submod_2(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_srcptr m)
{
    int cmp;
    mp_limb_t m1 = m[1], m0 = m[0];
    mp_limb_t x1 = x[1], x0 = x[0];
    mp_limb_t y1 = y[1], y0 = y[0];
    cmp = (x1 < y1) || (x1 == y1 && x0 < y0);
    sub_ddmmss(res[1], res[0], x1, x0, y1, y0);
    if (cmp)
        add_ssaaaa(res[1], res[0], res[1], res[0], m1, m0);
}

FLINT_FORCE_INLINE int
flint_mpn_signed_sub_n(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_size_t n)
{
    if (mpn_cmp(x, y, n) >= 0)
    {
        mpn_sub_n(res, x, y, n);
        return 0;
    }
    else
    {
        mpn_sub_n(res, y, x, n);
        return 1;
    }
}

/* add without carry in or carry out */
#define NN_ADD_2(r, u, v) add_ssaaaa((r)[1], (r)[0], (u)[1], (u)[0], (v)[1], (v)[0])
#define NN_ADD_3(r, u, v) add_sssaaaaaa((r)[2], (r)[1], (r)[0], (u)[2], (u)[1], (u)[0], (v)[2], (v)[1], (v)[0])
#define NN_ADD_4(r, u, v) add_ssssaaaaaaaa((r)[3], (r)[2], (r)[1], (r)[0], (u)[3], (u)[2], (u)[1], (u)[0], (v)[3], (v)[2], (v)[1], (v)[0])
#define NN_ADD_5(r, u, v) add_sssssaaaaaaaaaa((r)[4], (r)[3], (r)[2], (r)[1], (r)[0], (u)[4], (u)[3], (u)[2], (u)[1], (u)[0], (v)[4], (v)[3], (v)[2], (v)[1], (v)[0])
#define NN_ADD_6(r, u, v) add_ssssssaaaaaaaaaaaa((r)[5], (r)[4], (r)[3], (r)[2], (r)[1], (r)[0], (u)[5], (u)[4], (u)[3], (u)[2], (u)[1], (u)[0], (v)[5], (v)[4], (v)[3], (v)[2], (v)[1], (v)[0])
#define NN_ADD_7(r, u, v) add_sssssssaaaaaaaaaaaaaa((r)[6], (r)[5], (r)[4], (r)[3], (r)[2], (r)[1], (r)[0], (u)[6], (u)[5], (u)[4], (u)[3], (u)[2], (u)[1], (u)[0], (v)[6], (v)[5], (v)[4], (v)[3], (v)[2], (v)[1], (v)[0])
#define NN_ADD_8(r, u, v) add_ssssssssaaaaaaaaaaaaaaaa((r)[7], (r)[6], (r)[5], (r)[4], (r)[3], (r)[2], (r)[1], (r)[0], (u)[7], (u)[6], (u)[5], (u)[4], (u)[3], (u)[2], (u)[1], (u)[0], (v)[7], (v)[6], (v)[5], (v)[4], (v)[3], (v)[2], (v)[1], (v)[0])

#define NN_SUB_2(r, u, v) sub_ddmmss((r)[1], (r)[0], (u)[1], (u)[0], (v)[1], (v)[0])
#define NN_SUB_3(r, u, v) sub_dddmmmsss((r)[2], (r)[1], (r)[0], (u)[2], (u)[1], (u)[0], (v)[2], (v)[1], (v)[0])
#define NN_SUB_4(r, u, v) sub_ddddmmmmssss((r)[3], (r)[2], (r)[1], (r)[0], (u)[3], (u)[2], (u)[1], (u)[0], (v)[3], (v)[2], (v)[1], (v)[0])
#define NN_SUB_5(r, u, v) sub_dddddmmmmmsssss((r)[4], (r)[3], (r)[2], (r)[1], (r)[0], (u)[4], (u)[3], (u)[2], (u)[1], (u)[0], (v)[4], (v)[3], (v)[2], (v)[1], (v)[0])
#define NN_SUB_6(r, u, v) sub_ddddddmmmmmmssssss((r)[5], (r)[4], (r)[3], (r)[2], (r)[1], (r)[0], (u)[5], (u)[4], (u)[3], (u)[2], (u)[1], (u)[0], (v)[5], (v)[4], (v)[3], (v)[2], (v)[1], (v)[0])
#define NN_SUB_7(r, u, v) sub_dddddddmmmmmmmsssssss((r)[6], (r)[5], (r)[4], (r)[3], (r)[2], (r)[1], (r)[0], (u)[6], (u)[5], (u)[4], (u)[3], (u)[2], (u)[1], (u)[0], (v)[6], (v)[5], (v)[4], (v)[3], (v)[2], (v)[1], (v)[0])
#define NN_SUB_8(r, u, v) sub_ddddddddmmmmmmmmssssssss((r)[7], (r)[6], (r)[5], (r)[4], (r)[3], (r)[2], (r)[1], (r)[0], (u)[7], (u)[6], (u)[5], (u)[4], (u)[3], (u)[2], (u)[1], (u)[0], (v)[7], (v)[6], (v)[5], (v)[4], (v)[3], (v)[2], (v)[1], (v)[0])

#define DEF_SIGNED_SUB(n) \
FLINT_FORCE_INLINE int \
flint_mpn_signed_sub_ ## n(mp_ptr res, mp_srcptr x, mp_srcptr y) \
{ \
    if (mpn_cmp(x, y, n) >= 0) \
    { \
        NN_SUB_ ## n(res, x, y); \
        return 0; \
    } \
    else \
    { \
        NN_SUB_ ## n(res, y, x); \
        return 1; \
    } \
}

DEF_SIGNED_SUB(2)
DEF_SIGNED_SUB(3)
DEF_SIGNED_SUB(4)
DEF_SIGNED_SUB(5)
DEF_SIGNED_SUB(6)
DEF_SIGNED_SUB(7)
DEF_SIGNED_SUB(8)

FLINT_FORCE_INLINE void
flint_mpn_signed_div2(mp_ptr res, mp_srcptr x, mp_size_t n)
{
    mp_limb_t s = x[n - 1] & (UWORD(1) << (FLINT_BITS - 1));
    mpn_rshift(res, x, n, 1);
    res[n - 1] |= s;
}

void flint_mpn_mulmod_preinvn_2(mp_ptr r,
        mp_srcptr a, mp_srcptr b,
        mp_srcptr d, mp_srcptr dinv, ulong norm);

char * _flint_mpn_get_str(mp_srcptr x, mp_size_t n);


#define MPN_NORM(a, an)                         \
    do {                                        \
        while ((an) != 0 && (a)[(an) - 1] == 0) \
           (an)--;                              \
    } while (0)

#define MPN_SWAP(a, an, b, bn) \
  do { \
    FLINT_SWAP(mp_ptr, a, b); \
    FLINT_SWAP(mp_size_t, an, bn); \
  } while (0)

#define BITS_TO_LIMBS(b) (((b) + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS)

/* mpn macros for few-limb arithmetic ****************************************/

#define MPN_LEFT_SHIFT_HI(hi, lo, n) \
    ((n) > 0 ? (((hi) << (n)) | ((lo) >> (GMP_LIMB_BITS - (n)))) : (hi))

#define MPN_RIGHT_SHIFT_LOW(hi, lo, n) \
    ((n) > 0 ? (((lo) >> (n)) | ((hi) << (GMP_LIMB_BITS - (n)))) : (lo))

#define FLINT_MPN_MUL_2X1(r2, r1, r0, a1, a0, b0)           \
    do {                                                    \
        mp_limb_t t1;                                       \
        umul_ppmm(r1, r0, a0, b0);                          \
        umul_ppmm(r2, t1, a1, b0);                          \
        add_ssaaaa(r2, r1, r2, r1, 0, t1);                  \
    } while (0)

#define FLINT_MPN_MUL_2X2(r3, r2, r1, r0, a1, a0, b1, b0)   \
    do {                                                                  \
        mp_limb_t __v1, __v2, __u1, __u2;                                 \
        mp_limb_t __r3, __r2, __r1, __r0;                                 \
        mp_limb_t __a1 = (a1), __a0 = (a0), __b1 = (b1), __b0 = (b0);     \
        umul_ppmm(__r1, __r0, __a0, __b0);                                \
        umul_ppmm(__r3, __r2, __a1, __b1);                                \
        umul_ppmm(__v2, __v1, __a0, __b1);                                \
        add_sssaaaaaa(__r3, __r2, __r1, __r3, __r2, __r1, 0, __v2, __v1); \
        umul_ppmm(__u2, __u1, __a1, __b0);                                \
        add_sssaaaaaa(__r3, __r2, __r1, __r3, __r2, __r1, 0, __u2, __u1); \
        (r0) = __r0; (r1) = __r1; (r2) = __r2; (r3) = __r3;               \
    } while (0)

/* Low three words of 2x2 product */
#define FLINT_MPN_MUL_3P2X2(r2, r1, r0, a1, a0, b1, b0)                \
    do {                                                               \
        mp_limb_t __t1, __t2, __u1, __u2;                              \
        mp_limb_t __r2, __r1, __r0;                                    \
        mp_limb_t __a1 = (a1), __a0 = (a0), __b1 = (b1), __b0 = (b0);  \
        umul_ppmm(__r1, __r0, __a0, __b0);                             \
        __r2 = __a1 * __b1;                                            \
        umul_ppmm(__t2, __t1, __a0, __b1);                             \
        add_ssaaaa(__r2, __r1, __r2, __r1, __t2, __t1);                \
        umul_ppmm(__u2, __u1, __a1, __b0);                             \
        add_ssaaaa(__r2, __r1, __r2, __r1, __u2, __u1);                \
        (r0) = __r0; (r1) = __r1; (r2) = __r2;                         \
    } while (0)

#define FLINT_MPN_SQR_2X2(r3, r2, r1, r0, a1, a0)   \
    do {                                                                     \
        mp_limb_t __u1, __u2, __u3;                                          \
        mp_limb_t __r3, __r2, __r1, __r0;                                    \
        mp_limb_t __a1 = (a1), __a0 = (a0);                                  \
        umul_ppmm(__u2, __u1, __a0, __a1);                                   \
        add_sssaaaaaa(__u3, __u2, __u1, 0, __u2, __u1, 0, __u2, __u1);       \
        umul_ppmm(__r1, __r0, __a0, __a0);                                   \
        umul_ppmm(__r3, __r2, __a1, __a1);                                   \
        add_sssaaaaaa(__r3, __r2, __r1, __r3, __r2, __r1, __u3, __u2, __u1); \
        (r0) = __r0; (r1) = __r1; (r2) = __r2; (r3) = __r3;                  \
    } while (0)


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
            umul_ppmm(__dt1, __dt0, (u)[__i], (v)[(n) - 1 - __i]); \
            add_sssaaaaaa(__ds2, __ds1, __ds0, __ds2, __ds1, __ds0, 0, __dt1, __dt0); \
        } \
        (s0) = __ds0; (s1) = __ds1; (s2) = __ds2; \
    } while (0) \

/* Like NN_DOTREV_S3_1X1 but summing only over the high parts of the products. */
#define NN_DOTREV_S3_1X1_HIGH(s2, s1, u, v, n) \
    do { \
        mp_limb_t __dt0, __dt1, __ds0, __ds1, __ds2; \
        slong __i; \
        FLINT_ASSERT((n) >= 2); \
        umul_ppmm(__ds1, __ds0, (u)[0], (v)[(n) - 1]); \
        umul_ppmm(__dt1, __dt0, (u)[1], (v)[(n) - 2]); \
        add_ssaaaa(__ds2, __ds1, 0, __ds1, 0, __dt1); \
        for (__i = 2; __i < (n); __i++) \
        { \
            umul_ppmm(__dt1, __dt0, (u)[__i], (v)[(n) - 1 - __i]); \
            add_ssaaaa(__ds2, __ds1, __ds2, __ds1, 0, __dt1); \
        } \
        (s1) = __ds1; (s2) = __ds2; \
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

#define flint_mpn_divrem21_preinv(q, a_hi, a_lo, dinv) \
   do { \
      mp_limb_t __q2, __q3, __q4; \
      umul_ppmm((q), __q2, (a_hi), (dinv)); \
      umul_ppmm(__q3, __q4, (a_lo), (dinv)); \
      add_ssaaaa((q), __q2, (q), __q2, 0, __q3); \
      add_ssaaaa((q), __q2, (q), __q2, (a_hi), (a_lo)); \
   } while (0)

/* addition ******************************************************************/

#if FLINT_HAVE_ASSEMBLY_x86_64_adx
/* Simultaneously adds two n-limbed integers onto result and returns carry. */
/* NOTE: Requires n >= 4 */
# define FLINT_HAVE_NATIVE_mpn_2add_n_inplace 1
mp_limb_t flint_mpn_2add_n_inplace(mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);
#endif

#if FLINT_HAVE_NATIVE_mpn_add_nc
# define mpn_add_nc __gmpn_add_nc
mp_limb_t mpn_add_nc(mp_ptr, mp_srcptr, mp_srcptr, mp_size_t, mp_limb_t);
#else
FLINT_FORCE_INLINE mp_limb_t
mpn_add_nc(mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n, mp_limb_t ci)
{
    mp_limb_t co;
    co = mpn_add_n(rp, up, vp, n);
    co += mpn_add_1(rp, rp, n, ci);
    return co;
}
#endif

#if FLINT_HAVE_NATIVE_mpn_sub_nc
# define mpn_sub_nc __gmpn_sub_nc
mp_limb_t mpn_sub_nc(mp_ptr, mp_srcptr, mp_srcptr, mp_size_t, mp_limb_t);
#else
FLINT_FORCE_INLINE mp_limb_t
mpn_sub_nc(mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n, mp_limb_t ci)
{
    mp_limb_t co;
    co = mpn_sub_n(rp, up, vp, n);
    co += mpn_sub_1(rp, rp, n, ci);
    return co;
}
#endif

#if FLINT_HAVE_NATIVE_mpn_add_n_sub_n
/* mpn_add_n_sub_n basically only exists for IA64 and certain PowerPC and s390
 * systems. We will assume that a native one does not exist. */
# undef FLINT_HAVE_NATIVE_mpn_add_n_sub_n
# define FLINT_HAVE_NATIVE_mpn_add_n_sub_n 0
#endif

mp_limb_t flint_mpn_sumdiff_n(mp_ptr s, mp_ptr d, mp_srcptr x, mp_srcptr y, mp_size_t n);

/* shifting ******************************************************************/

#if FLINT_HAVE_NATIVE_mpn_addlsh1_n
# define mpn_addlsh1_n __gmpn_addlsh1_n
mp_limb_t mpn_addlsh1_n(mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);
#endif

#if FLINT_HAVE_NATIVE_mpn_addlsh1_n_ip1
# define mpn_addlsh1_n_ip1 __gmpn_addlsh1_n_ip1
mp_limb_t mpn_addlsh1_n_ip1(mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);
#elif FLINT_HAVE_NATIVE_mpn_addlsh1_n
# define mpn_addlsh1_n_ip1(a,b,n) mpn_addlsh1_n(a,a,b,n)
# define FLINT_HAVE_NATIVE_mpn_addlsh1_n_ip1 2
#endif

#if FLINT_HAVE_NATIVE_mpn_rsh1add_n
# define mpn_rsh1add_n __gmpn_rsh1add_n
mp_limb_t mpn_rsh1add_n(mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);
#endif

#if FLINT_HAVE_NATIVE_mpn_rsh1sub_n
# define mpn_rsh1sub_n __gmpn_rsh1sub_n
mp_limb_t mpn_rsh1sub_n(mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);
#endif

/* multiplication (general) **************************************************/

#if FLINT_HAVE_ASSEMBLY_x86_64_adx
# define FLINT_MPN_MUL_FUNC_TAB_WIDTH 17
# define FLINT_MPN_SQR_FUNC_TAB_WIDTH 14

# define FLINT_HAVE_MUL_FUNC(n, m) ((n) <= 16)
# define FLINT_HAVE_MUL_N_FUNC(n) ((n) <= 16)
# define FLINT_HAVE_SQR_FUNC(n) ((n) <= FLINT_MPN_SQR_FUNC_TAB_WIDTH)

# define FLINT_MPN_MUL_HARD(rp, xp, xn, yp, yn) (flint_mpn_mul_func_tab[xn][yn](rp, xp, yp))
# define FLINT_MPN_MUL_N_HARD(rp, xp, yp, n) (flint_mpn_mul_n_func_tab[n](rp, xp, yp))
# define FLINT_MPN_SQR_HARD(rp, xp, n) (flint_mpn_sqr_func_tab[n](rp, xp))
#elif FLINT_HAVE_ASSEMBLY_armv8
# define FLINT_MPN_MUL_FUNC_N_TAB_WIDTH 15
# define FLINT_MPN_SQR_FUNC_TAB_WIDTH 9

# define FLINT_HAVE_MUL_FUNC(n, m) FLINT_HAVE_MUL_N_FUNC(n)
# define FLINT_HAVE_MUL_N_FUNC(n) ((n) <= FLINT_MPN_MUL_FUNC_N_TAB_WIDTH)
# define FLINT_HAVE_SQR_FUNC(n) ((n) <= FLINT_MPN_SQR_FUNC_TAB_WIDTH)

# define FLINT_MPN_MUL_HARD(rp, xp, xn, yp, yn) (flint_mpn_mul_func_n_tab[xn](rp, xp, yp, yn))
# define FLINT_MPN_MUL_N_HARD(rp, xp, yp, n) (flint_mpn_mul_func_n_tab[n](rp, xp, yp, n))
# define FLINT_MPN_SQR_HARD(rp, xp, n) (flint_mpn_sqr_func_tab[n](rp, xp))

# define FLINT_HAVE_NATIVE_mpn_mul_2 1
mp_limb_t flint_mpn_mul_2(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr);
#else
# define FLINT_MPN_MUL_FUNC_TAB_WIDTH 8
# define FLINT_MPN_SQR_FUNC_TAB_WIDTH 0

# define FLINT_HAVE_MUL_FUNC(n, m) ((n) <= 7 || ((n) <= 14 && (m) == 1))
# define FLINT_HAVE_MUL_N_FUNC(n) ((n) <= 7)
# define FLINT_HAVE_SQR_FUNC(n) (0)

# define FLINT_MPN_MUL_HARD(rp, xp, xn, yp, yn) (flint_mpn_mul_func_tab[xn][yn](rp, xp, yp))
# define FLINT_MPN_MUL_N_HARD(rp, xp, yp, n) (flint_mpn_mul_n_func_tab[n](rp, xp, yp))
# define FLINT_MPN_SQR_HARD(rp, xp, n) (flint_mpn_sqr_func_tab[n](rp, xp))
#endif

#define FLINT_MUL_USE_FUNC_TAB 1

typedef mp_limb_t (* flint_mpn_mul_func_t)(mp_ptr, mp_srcptr, mp_srcptr);
typedef mp_limb_t (* flint_mpn_mul_func_n_t)(mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);
typedef mp_limb_t (* flint_mpn_sqr_func_t)(mp_ptr, mp_srcptr);

#ifdef FLINT_MPN_MUL_FUNC_N_TAB_WIDTH
FLINT_DLL extern const flint_mpn_mul_func_n_t flint_mpn_mul_func_n_tab[];
#else
FLINT_DLL extern const flint_mpn_mul_func_t flint_mpn_mul_func_tab[][FLINT_MPN_MUL_FUNC_TAB_WIDTH];
FLINT_DLL extern const flint_mpn_mul_func_t flint_mpn_mul_n_func_tab[];
#endif

FLINT_DLL extern const flint_mpn_sqr_func_t flint_mpn_sqr_func_tab[];

void flint_mpn_mul_toom22(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t, mp_ptr);
void flint_mpn_mul_toom32(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t, mp_ptr);

mp_limb_t _flint_mpn_mul(mp_ptr r, mp_srcptr x, mp_size_t xn, mp_srcptr y, mp_size_t yn);
void _flint_mpn_mul_n(mp_ptr r, mp_srcptr x, mp_srcptr y, mp_size_t n);
mp_limb_t _flint_mpn_sqr(mp_ptr r, mp_srcptr x, mp_size_t n);

MPN_EXTRAS_INLINE mp_limb_t
flint_mpn_mul(mp_ptr r, mp_srcptr x, mp_size_t xn, mp_srcptr y, mp_size_t yn)
{
    FLINT_ASSERT(xn >= yn);
    FLINT_ASSERT(yn >= 1);
    FLINT_ASSERT(r != x);
    FLINT_ASSERT(r != y);

    if (FLINT_MUL_USE_FUNC_TAB && FLINT_HAVE_MUL_FUNC(xn, yn))
        return FLINT_MPN_MUL_HARD(r, x, xn, y, yn);
    else
        return _flint_mpn_mul(r, x, xn, y, yn);
}

MPN_EXTRAS_INLINE void
flint_mpn_mul_n(mp_ptr r, mp_srcptr x, mp_srcptr y, mp_size_t n)
{
    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(r != x);
    FLINT_ASSERT(r != y);

    if (FLINT_MUL_USE_FUNC_TAB && FLINT_HAVE_MUL_N_FUNC(n))
        FLINT_MPN_MUL_N_HARD(r, x, y, n);
    else
        _flint_mpn_mul_n(r, x, y, n);
}

MPN_EXTRAS_INLINE mp_limb_t
flint_mpn_sqr(mp_ptr r, mp_srcptr x, mp_size_t n)
{
    FLINT_ASSERT(n >= 1);

    if (FLINT_MUL_USE_FUNC_TAB && FLINT_HAVE_SQR_FUNC(n))
    {
        /* NOTE: Aliasing allowed */
        return FLINT_MPN_SQR_HARD(r, x, n);
    }
    else
    {
        FLINT_ASSERT(r != x);
        return _flint_mpn_sqr(r, x, n);
    }
}

/* Like flint_mpn_mul but allow operands in either order, completely
   inline some small products, and also check for squaring. */
#define FLINT_MPN_MUL_WITH_SPECIAL_CASES(_z, _x, _xn, _y, _yn) \
    if ((_xn) == (_yn)) \
    { \
        if ((_xn) == 1) \
        { \
            umul_ppmm((_z)[1], (_z)[0], (_x)[0], (_y)[0]); \
        } \
        else if ((_xn) == 2) \
        { \
            mp_limb_t __tt_x1, __tt_x0, __tt_y1, __tt_y0; \
            __tt_x0 = (_x)[0]; \
            __tt_x1 = (_x)[1]; \
            __tt_y0 = (_y)[0]; \
            __tt_y1 = (_y)[1]; \
            FLINT_MPN_MUL_2X2((_z)[3], (_z)[2], (_z)[1], (_z)[0], __tt_x1, __tt_x0, __tt_y1, __tt_y0); \
        } \
        else if ((_x) == (_y)) \
            flint_mpn_sqr((_z), (_x), (_xn)); \
        else \
            flint_mpn_mul_n((_z), (_x), (_y), (_xn)); \
    } \
    else if ((_xn) > (_yn)) \
    { \
        flint_mpn_mul((_z), (_x), (_xn), (_y), (_yn)); \
    } \
    else \
    { \
        flint_mpn_mul((_z), (_y), (_yn), (_x), (_xn)); \
    }

/* High and low multiplication *******************************************************/

#define FLINT_HAVE_MULLOW_FUNC(n) ((n) <= FLINT_MPN_MULLOW_FUNC_TAB_WIDTH)
#define FLINT_HAVE_MULHIGH_FUNC(n) ((n) <= FLINT_MPN_MULHIGH_FUNC_TAB_WIDTH)
#define FLINT_HAVE_SQRHIGH_FUNC(n) ((n) <= FLINT_MPN_SQRHIGH_FUNC_TAB_WIDTH)
#define FLINT_HAVE_MULHIGH_NORMALISED_FUNC(n) ((n) <= FLINT_MPN_MULHIGH_NORMALISED_FUNC_TAB_WIDTH)
#define FLINT_HAVE_SQRHIGH_NORMALISED_FUNC(n) ((n) <= FLINT_MPN_SQRHIGH_NORMALISED_FUNC_TAB_WIDTH)

typedef struct { mp_limb_t m1; mp_limb_t m2; } mp_limb_pair_t;
typedef mp_limb_pair_t (* flint_mpn_sqrhigh_normalised_func_t)(mp_ptr, mp_srcptr);
typedef mp_limb_pair_t (* flint_mpn_mulhigh_normalised_func_t)(mp_ptr, mp_srcptr, mp_srcptr);

FLINT_DLL extern const flint_mpn_mul_func_t flint_mpn_mullow_func_tab[];
FLINT_DLL extern const flint_mpn_mul_func_t flint_mpn_mulhigh_func_tab[];
FLINT_DLL extern const flint_mpn_sqr_func_t flint_mpn_sqrhigh_func_tab[];
FLINT_DLL extern const flint_mpn_mulhigh_normalised_func_t flint_mpn_mulhigh_normalised_func_tab[];
FLINT_DLL extern const flint_mpn_sqrhigh_normalised_func_t flint_mpn_sqrhigh_normalised_func_tab[];

#if FLINT_HAVE_ASSEMBLY_x86_64_adx
# define FLINT_MPN_MULLOW_FUNC_TAB_WIDTH 8
# define FLINT_MPN_MULHIGH_FUNC_TAB_WIDTH 13
/* n with best effective cycles/limb (and current largest assembly case) -- used by mulhigh_recursive */
# define FLINT_MPN_MULHIGH_BEST_TAB_N 9
# define FLINT_MPN_SQRHIGH_FUNC_TAB_WIDTH 8
# define FLINT_MPN_MULHIGH_NORMALISED_FUNC_TAB_WIDTH 9
# define FLINT_MPN_SQRHIGH_NORMALISED_FUNC_TAB_WIDTH 8

# define FLINT_HAVE_NATIVE_mpn_mullow_basecase 1
/* NOTE: This function only works for n >= 6 */
# define FLINT_HAVE_NATIVE_mpn_mulhigh_basecase 1
/* NOTE: This function only works for n >= 6 */
# define FLINT_HAVE_NATIVE_mpn_sqrhigh_basecase 1

#elif FLINT_HAVE_ASSEMBLY_armv8
# define FLINT_MPN_MULLOW_FUNC_TAB_WIDTH 0
# define FLINT_MPN_MULHIGH_FUNC_TAB_WIDTH 8
# define FLINT_MPN_MULHIGH_BEST_TAB_N 8
# define FLINT_MPN_SQRHIGH_FUNC_TAB_WIDTH 8
# define FLINT_MPN_MULHIGH_NORMALISED_FUNC_TAB_WIDTH 0
# define FLINT_MPN_SQRHIGH_NORMALISED_FUNC_TAB_WIDTH 0

/* NOTE: This function only works for n > 8 */
# define FLINT_HAVE_NATIVE_mpn_mulhigh_basecase 1

#else
/* TODO: generic hardcoded mullows */
# define FLINT_MPN_MULLOW_FUNC_TAB_WIDTH 0
# define FLINT_MPN_MULHIGH_FUNC_TAB_WIDTH 16
# define FLINT_MPN_MULHIGH_BEST_TAB_N 16
# define FLINT_MPN_SQRHIGH_FUNC_TAB_WIDTH 2
# define FLINT_MPN_MULHIGH_NORMALISED_FUNC_TAB_WIDTH 0
# define FLINT_MPN_SQRHIGH_NORMALISED_FUNC_TAB_WIDTH 0

#endif

/* FIXME: this tuning is for x86_64_adx with fft_small */
/* FIXME: we currently assume that the same parameters are optimal for both mulhigh and mullow */
#define FLINT_MPN_MULLOW_MULDERS_CUTOFF 50
#define FLINT_MPN_MULHIGH_MULDERS_CUTOFF 40
#define FLINT_MPN_MULHIGH_MUL_CUTOFF 2000
#define FLINT_MPN_MULHIGH_K_TAB_SIZE 2048

FLINT_DLL extern const signed short flint_mpn_mulhigh_k_tab[FLINT_MPN_MULHIGH_K_TAB_SIZE];

mp_limb_t flint_mpn_mullow_basecase(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n);
void _flint_mpn_mullow_n_mulders_recursive(mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n);
mp_limb_t _flint_mpn_mullow_n_mulders(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n);
mp_limb_t _flint_mpn_mullow_n_mul(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n);
mp_limb_t _flint_mpn_mullow_n(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n);

mp_limb_t _flint_mpn_mulhigh_basecase(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n);
void _flint_mpn_mulhigh_n_mulders_recursive(mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n);
mp_limb_t _flint_mpn_mulhigh_n_naive(mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n);
mp_limb_t _flint_mpn_mulhigh_n_recursive(mp_ptr r, mp_srcptr x, mp_srcptr y, mp_size_t n);
mp_limb_t _flint_mpn_mulhigh_n_mulders(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n);
mp_limb_t _flint_mpn_mulhigh_n_mul(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n);
mp_limb_t _flint_mpn_mulhigh_n(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n);

#if FLINT_HAVE_ASSEMBLY_x86_64_adx
MPN_EXTRAS_INLINE
mp_limb_t _flint_mpn_mulhigh_n_basecase2(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n)
{
    if (n <= 22)
        return _flint_mpn_mulhigh_n_recursive(rp, xp, yp, n);
    else
        return _flint_mpn_mulhigh_basecase(rp, xp, yp, n);
}
#else
#define _flint_mpn_mulhigh_n_basecase2 _flint_mpn_mulhigh_basecase
#endif

MPN_EXTRAS_INLINE
mp_limb_t flint_mpn_mullow_n(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n)
{
    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(rp != xp);

    if (FLINT_HAVE_MULLOW_FUNC(n))
        return flint_mpn_mullow_func_tab[n](rp, xp, yp);
    else
        return _flint_mpn_mullow_n(rp, xp, yp, n);
}

MPN_EXTRAS_INLINE
mp_limb_t flint_mpn_mulhigh_n(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n)
{
    FLINT_ASSERT(n >= 1);

    if (FLINT_HAVE_MULHIGH_FUNC(n)) /* NOTE: Aliasing allowed here */
        return flint_mpn_mulhigh_func_tab[n](rp, xp, yp);
    else
        return _flint_mpn_mulhigh_n(rp, xp, yp, n);
}

/* We just want the high or low n limbs, but rp has 2n limbs available
   which can be used for scratch space or for doing a full multiply
   without temporary allocations. TODO: exploit this in the Mulders range
   by calling Mulders directly. */
MPN_EXTRAS_INLINE
void flint_mpn_mul_or_mullow_n(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n)
{
    FLINT_ASSERT(n >= 1);

    if (FLINT_HAVE_MULLOW_FUNC(n))
        rp[n] = flint_mpn_mullow_func_tab[n](rp, xp, yp);
    else if (n < FLINT_MPN_MULHIGH_MUL_CUTOFF)
        rp[n] = _flint_mpn_mullow_n(rp, xp, yp, n);
    else
        flint_mpn_mul_n(rp, xp, yp, n);
}

MPN_EXTRAS_INLINE
void flint_mpn_mul_or_mulhigh_n(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n)
{
    FLINT_ASSERT(n >= 1);

    if (FLINT_HAVE_MULHIGH_FUNC(n))
        rp[n - 1] = flint_mpn_mulhigh_func_tab[n](rp + n, xp, yp);
    else if (n < FLINT_MPN_MULHIGH_MUL_CUTOFF)
        rp[n - 1] = _flint_mpn_mulhigh_n(rp + n, xp, yp, n);
    else
        flint_mpn_mul_n(rp, xp, yp, n);
}

#define FLINT_MPN_SQRHIGH_MULDERS_CUTOFF 90
#define FLINT_MPN_SQRHIGH_SQR_CUTOFF 2000
#define FLINT_MPN_SQRHIGH_K_TAB_SIZE 2048

#if FLINT_HAVE_ASSEMBLY_x86_64_adx
mp_limb_t _flint_mpn_sqrhigh_basecase_even(mp_ptr, mp_srcptr, mp_size_t);
mp_limb_t _flint_mpn_sqrhigh_basecase_odd(mp_ptr, mp_srcptr, mp_size_t);

MPN_EXTRAS_INLINE mp_limb_t _flint_mpn_sqrhigh_basecase(mp_ptr rp, mp_srcptr xp, mp_size_t n)
{
    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(rp != xp);

    if (n & 1)
        return _flint_mpn_sqrhigh_basecase_odd(rp, xp, n >> 1);
    else
        return _flint_mpn_sqrhigh_basecase_even(rp, xp, n >> 1);
}

#else
/* todo */
MPN_EXTRAS_INLINE mp_limb_t _flint_mpn_sqrhigh_basecase(mp_ptr res, mp_srcptr u, mp_size_t n)
{
    return _flint_mpn_mulhigh_basecase(res, u, u, n);
}
#endif

void _flint_mpn_sqrhigh_mulders_recursive(mp_ptr rp, mp_srcptr np, mp_size_t n);
mp_limb_t _flint_mpn_sqrhigh_mulders(mp_ptr res, mp_srcptr u, mp_size_t n);
mp_limb_t _flint_mpn_sqrhigh_sqr(mp_ptr res, mp_srcptr u, mp_size_t n);
mp_limb_t _flint_mpn_sqrhigh(mp_ptr, mp_srcptr, mp_size_t);

MPN_EXTRAS_INLINE
mp_limb_t flint_mpn_sqrhigh(mp_ptr rp, mp_srcptr xp, mp_size_t n)
{
    FLINT_ASSERT(n >= 1);

    if (FLINT_HAVE_SQRHIGH_FUNC(n)) /* NOTE: Aliasing allowed here */
        return flint_mpn_sqrhigh_func_tab[n](rp, xp);
    else
        return _flint_mpn_sqrhigh(rp, xp, n);
}

mp_limb_pair_t _flint_mpn_mulhigh_normalised(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n);

MPN_EXTRAS_INLINE
mp_limb_pair_t flint_mpn_mulhigh_normalised(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n)
{
    FLINT_ASSERT(n >= 1);

    if (FLINT_HAVE_MULHIGH_NORMALISED_FUNC(n))
        return flint_mpn_mulhigh_normalised_func_tab[n](rp, xp, yp);
    else
        return _flint_mpn_mulhigh_normalised(rp, xp, yp, n);
}

mp_limb_pair_t _flint_mpn_sqrhigh_normalised(mp_ptr rp, mp_srcptr xp, mp_size_t n);

MPN_EXTRAS_INLINE
mp_limb_pair_t flint_mpn_sqrhigh_normalised(mp_ptr rp, mp_srcptr xp, mp_size_t n)
{
    FLINT_ASSERT(n >= 1);

    if (FLINT_HAVE_SQRHIGH_NORMALISED_FUNC(n))
        return flint_mpn_sqrhigh_normalised_func_tab[n](rp, xp);
    else
        return _flint_mpn_sqrhigh_normalised(rp, xp, n);
}

/* division ******************************************************************/

#if FLINT_HAVE_NATIVE_mpn_modexact_1_odd
# define mpn_modexact_1_odd __gmpn_modexact_1_odd
mp_limb_t mpn_modexact_1_odd(mp_srcptr, mp_size_t, mp_limb_t);
#endif

#if FLINT_HAVE_NATIVE_mpn_invert_limb
# define mpn_invert_limb __gmpn_invert_limb
mp_limb_t mpn_invert_limb(mp_limb_t);
#endif

mp_limb_t flint_mpn_preinv1(mp_limb_t d, mp_limb_t d2);
void flint_mpn_preinvn(mp_ptr dinv, mp_srcptr d, mp_size_t n);

#if defined(mpn_modexact_1_odd)
MPN_EXTRAS_INLINE
int flint_mpn_divisible_1_odd(mp_srcptr x, mp_size_t xsize, mp_limb_t d)
{
    return mpn_modexact_1_odd(x, xsize, d) == 0;
}
#else
# include "gmpcompat.h"
MPN_EXTRAS_INLINE
int flint_mpn_divisible_1_odd(mp_srcptr x, mp_size_t xsize, mp_limb_t d)
{
    __mpz_struct s;
    s._mp_size = xsize;
    s._mp_d = (mp_ptr) x;
    return flint_mpz_divisible_ui_p(&s, d);
}
#endif

FLINT_FORCE_INLINE
void mpn_tdiv_q(mp_ptr qp, mp_srcptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn)
{
    mp_ptr _scratch;
    TMP_INIT;
    TMP_START;
    _scratch = (mp_ptr) TMP_ALLOC(dn * sizeof(mp_limb_t));
    mpn_tdiv_qr(qp, _scratch, 0, np, nn, dp, dn);
    TMP_END;
}

int flint_mpn_divides(mp_ptr q, mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn, mp_ptr scr);

void flint_mpn_mod_preinvn(mp_ptr r, mp_srcptr a, mp_size_t m, mp_srcptr d, mp_size_t n, mp_srcptr dinv);

mp_limb_t flint_mpn_divrem_preinv1(mp_ptr q, mp_ptr a, mp_size_t m, mp_srcptr b, mp_size_t n, mp_limb_t dinv);
mp_limb_t flint_mpn_divrem_preinvn(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t m, mp_srcptr d, mp_size_t n, mp_srcptr dinv);

/* composed arithmetic *******************************************************/

mp_size_t flint_mpn_fmms1(mp_ptr y, mp_limb_t a1, mp_srcptr x1, mp_limb_t a2, mp_srcptr x2, mp_size_t n);

/* debug *********************************************************************/

void flint_mpn_debug(mp_srcptr x, mp_size_t xsize);

/* factorisation *************************************************************/

mp_size_t flint_mpn_remove_2exp(mp_ptr x, mp_size_t xsize, flint_bitcnt_t * bits);

mp_size_t flint_mpn_remove_power_ascending(mp_ptr x, mp_size_t xsize, mp_ptr p, mp_size_t psize, ulong * exp);

int flint_mpn_factor_trial(mp_srcptr x, mp_size_t xsize, slong start, slong stop);

int flint_mpn_factor_trial_tree(slong * factors, mp_srcptr x, mp_size_t xsize, slong num_primes);

/* greatest common divisor ***************************************************/

mp_size_t flint_mpn_gcd_full2(mp_ptr gp, mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn, mp_ptr scr);
mp_size_t flint_mpn_gcd_full(mp_ptr gp, mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn);

/* modular arithmetic ********************************************************/

void flint_mpn_mulmod_preinv1(mp_ptr r, mp_srcptr a, mp_srcptr b, mp_size_t n, mp_srcptr d, mp_limb_t dinv, ulong norm);
void flint_mpn_mulmod_preinvn(mp_ptr r, mp_srcptr a, mp_srcptr b, mp_size_t n, mp_srcptr d, mp_srcptr dinv, ulong norm);

int flint_mpn_mulmod_2expp1_basecase(mp_ptr xp, mp_srcptr yp, mp_srcptr zp, int c, flint_bitcnt_t b, mp_ptr tp);

/* miscellaneous *************************************************************/

double flint_mpn_get_d(mp_srcptr ptr, mp_size_t size, mp_size_t sign, long exp);

/* random ********************************************************************/

void flint_mpn_rrandom(mp_ptr rp, flint_rand_t state, mp_size_t n);
void flint_mpn_urandomb(mp_ptr rp, flint_rand_t state, flint_bitcnt_t n);

/******************************************************************************
    Divisions where the quotient is expected to be small. All function do:
        input: n > d > 0
        output: q = n/d, r = n%d
    for various small sizes of n and d.
    Not in a function because compiler refuses to inline eudiv_qrrnndd.
    Each macro takes a prefix t for its local vars.
******************************************************************************/

#define eudiv_qrnd(q, r, n, d, t)           \
do {                                        \
    mp_limb_t t##q, t##a = n, t##b = d;     \
                                            \
    FLINT_ASSERT(t##a > t##b);              \
    FLINT_ASSERT(t##b > 0);                 \
                                            \
    t##a -= t##b;                           \
    for (t##q = 1; t##q < 5; t##q++)        \
    {                                       \
        if (t##a < t##b)                    \
            goto t##quotient_found;         \
        t##a -= t##b;                       \
    }                                       \
    t##q += t##a / t##b;                    \
    t##a = t##a % t##b;                     \
                                            \
t##quotient_found:                          \
                                            \
    q = t##q;                               \
    r = t##a;                               \
} while (0)


#define eudiv_qqrnnd(q1, q0, r0, n1, n0, d0, t)         \
do {                                                    \
    mp_limb_t t##a1 = n1, t##a0 = n0, t##b0 = d0;       \
    mp_limb_t t##q1, t##q0, t##r0, t##u;                \
                                                        \
    FLINT_ASSERT(t##a1 > 0 || t##a0 > t##b0);           \
                                                        \
    udiv_qrnnd(t##q1, t##u, 0, t##a1, t##b0);           \
    udiv_qrnnd(t##q0, t##r0, t##u, t##a0, t##b0);       \
                                                        \
    q1 = t##q1;                                         \
    q0 = t##q0;                                         \
    r0 = t##r0;                                         \
} while (0)

/* d must be normalized, i.e. d1 != 0 */
#define eudiv_qrrnndd(q0, r1, r0, n1, n0, d1, d0, t)                        \
do {                                                                        \
    int t##i;                                                               \
    mp_limb_t t##a1 = n1, t##a0 = n0, t##b1 = d1, t##b0 = d0;               \
    mp_limb_t t##r1, t##r0, t##u2, t##u1, t##u0, t##q, t##v1, t##v0;        \
                                                                            \
    FLINT_ASSERT(t##a1 != 0);                                               \
    FLINT_ASSERT(t##b1 != 0);                                               \
    FLINT_ASSERT(t##b1 < t##a1 || (t##b1 == t##a1 && t##b0 < t##a0));       \
                                                                            \
    t##q = 1;                                                               \
                                                                            \
    sub_ddmmss(t##r1,t##r0, t##a1,t##a0, t##b1,t##b0);                      \
                                                                            \
t##subtract:                                                                \
                                                                            \
    for (t##i = 2; t##i <= 4; t##i++)                                       \
    {                                                                       \
        sub_dddmmmsss(t##u2,t##u1,t##u0, 0,t##r1,t##r0, 0,t##b1,t##b0);     \
        if (t##u2 != 0)                                                     \
            goto t##quotient_found;                                         \
        t##q += 1;                                                          \
        t##r0 = t##u0;                                                      \
        t##r1 = t##u1;                                                      \
    }                                                                       \
                                                                            \
    if (t##r1 != 0)                                                         \
    {                                                                       \
        int t##ncnt, t##dcnt;                                               \
        mp_limb_t t##qq = 0;                                                \
                                                                            \
        t##ncnt = flint_clz(t##r1);                                \
        t##dcnt = flint_clz(t##b1);                                \
        t##dcnt -= t##ncnt;                                                 \
        if (t##dcnt <= 0)                                                   \
            goto t##subtract;                                               \
                                                                            \
        t##v1 = (t##b1 << t##dcnt) | (t##b0 >> (FLINT_BITS - t##dcnt));     \
        t##v0 = t##b0 << t##dcnt;                                           \
                                                                            \
        do {                                                                \
            sub_dddmmmsss(t##u2,t##u1,t##u0, 0,t##r1,t##r0, 0,t##v1,t##v0); \
            t##qq = 2*t##qq + 1 + t##u2;                                    \
            t##r1 = t##u2 ? t##r1 : t##u1;                                  \
            t##r0 = t##u2 ? t##r0 : t##u0;                                  \
            t##v0 = (t##v1 << (FLINT_BITS - 1)) | (t##v0 >> 1);             \
            t##v1 = t##v1 >> 1;                                             \
            t##dcnt--;                                                      \
        } while (t##dcnt >= 0);                                             \
                                                                            \
        t##q += t##qq;                                                      \
    }                                                                       \
                                                                            \
t##quotient_found:                                                          \
                                                                            \
    FLINT_ASSERT(t##r1 < t##b1 || (t##r1 == t##b1 && t##r0 < t##b0));       \
                                                                            \
    q0 = t##q;                                                              \
    r0 = t##r0;                                                             \
    r1 = t##r1;                                                             \
} while (0)


#ifdef __cplusplus
}
#endif

#endif
