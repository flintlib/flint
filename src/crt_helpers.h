/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef CRT_HELPERS_H
#define CRT_HELPERS_H

#if defined(__GNUC__) && defined(__AVX2__)
# include <immintrin.h>
#elif defined(_MSC_VER) && defined(__AVX2__)
# include <intrin.h>
#endif

#include "longlong.h"
#include "templates.h"

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__AVX2__)

FLINT_FORCE_INLINE unsigned char _addcarry_ulong(unsigned char cf, ulong x, ulong y, ulong* s)
{
    long long unsigned int _s;
    cf = _addcarry_u64(cf, (long long unsigned int)(x),
                           (long long unsigned int)(y),
                           &_s);
    *s = (ulong)(_s);
    return cf;
}

FLINT_FORCE_INLINE unsigned char _subborrow_ulong(unsigned char cf, ulong x, ulong y, ulong* s)
{
    long long unsigned int _s;
    cf = _subborrow_u64(cf, (long long unsigned int)(x),
                            (long long unsigned int)(y),
                           &_s);
    *s = (ulong)(_s);
    return cf;
}

#else

FLINT_FORCE_INLINE unsigned char _addcarry_ulong(unsigned char cf, ulong x, ulong y, ulong* s)
{
#if 0
    ulong cf2;
    *s = __builtin_addcl(x, y, cf, &cf2);
    return cf2;
#else
    ulong hi, lo;
    add_ssaaaa(hi, lo, 0, x, 0, y);
    add_ssaaaa(hi, lo, hi, lo, 0, (ulong) cf);
    *s = lo;
    return hi;
#endif
}

FLINT_FORCE_INLINE unsigned char _subborrow_ulong(unsigned char cf, ulong x, ulong y, ulong* s)
{
#if 0
    ulong cf2;
    *s = __builtin_subcl(x, y, cf, &cf2);
    return cf2;
#else
    ulong hi, lo;
    sub_ddmmss(hi, lo, 0, x, 0, y);
    sub_ddmmss(hi, lo, hi, lo, 0, (ulong) cf);
    *s = lo;
    return hi != 0;
#endif
}

#endif


#if 1

#if defined(__GNUC__) && defined(__AVX2__)


#elif defined(__GNUC__) && defined(__ARM_NEON)



#elif defined(_MSC_VER) && (defined(__AVX2__) || defined(_M_ARM64))



#else
# error crt_helpers.h requires AVX2 or Neon instructions
#endif

FLINT_FORCE_INLINE void multi_add_0(ulong FLINT_UNUSED(z[]), const ulong FLINT_UNUSED(a[]))
{
}

FLINT_FORCE_INLINE void multi_add_1(ulong z[], const ulong a[])
{
    z[0] += a[0];
}

FLINT_FORCE_INLINE void multi_add_2(ulong z[], const ulong a[])
{
    add_ssaaaa(z[1],z[0],
               z[1],z[0],
               a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_add_3(ulong z[], const ulong a[])
{
    add_sssaaaaaa(z[2],z[1],z[0],
                  z[2],z[1],z[0],
                  a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_add_4(ulong z[], const ulong a[])
{
    add_ssssaaaaaaaa(z[3],z[2],z[1],z[0],
                     z[3],z[2],z[1],z[0],
                     a[3],a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_add_5(ulong z[], const ulong a[])
{
    add_sssssaaaaaaaaaa(z[4],z[3],z[2],z[1],z[0],
                        z[4],z[3],z[2],z[1],z[0],
                        a[4],a[3],a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_add_6(ulong z[], const ulong a[])
{
    add_ssssssaaaaaaaaaaaa(z[5],z[4],z[3],z[2],z[1],z[0],
                           z[5],z[4],z[3],z[2],z[1],z[0],
                           a[5],a[4],a[3],a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_add_7(ulong z[], const ulong a[])
{
    add_sssssssaaaaaaaaaaaaaa(z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                              z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                              a[6],a[5],a[4],a[3],a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_add_8(ulong z[], const ulong a[])
{
    add_ssssssssaaaaaaaaaaaaaaaa(z[7],z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                                 z[7],z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                                 a[7],a[6],a[5],a[4],a[3],a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_sub_0(ulong FLINT_UNUSED(z[]), const ulong FLINT_UNUSED(a[]))
{
}

FLINT_FORCE_INLINE void multi_sub_1(ulong z[], const ulong a[])
{
    z[0] -= a[0];
}

FLINT_FORCE_INLINE void multi_sub_2(ulong z[], const ulong a[])
{
    sub_ddmmss(z[1],z[0],
               z[1],z[0],
               a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_sub_3(ulong z[], const ulong a[])
{
    sub_dddmmmsss(z[2],z[1],z[0],
                  z[2],z[1],z[0],
                  a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_sub_4(ulong z[], const ulong a[])
{
    sub_ddddmmmmssss(z[3],z[2],z[1],z[0],
                     z[3],z[2],z[1],z[0],
                     a[3],a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_sub_5(ulong z[], const ulong a[])
{
    sub_dddddmmmmmsssss(z[4],z[3],z[2],z[1],z[0],
                        z[4],z[3],z[2],z[1],z[0],
                        a[4],a[3],a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_sub_6(ulong z[], const ulong a[])
{
    sub_ddddddmmmmmmssssss(z[5],z[4],z[3],z[2],z[1],z[0],
                           z[5],z[4],z[3],z[2],z[1],z[0],
                           a[5],a[4],a[3],a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_sub_7(ulong z[], const ulong a[])
{
    sub_dddddddmmmmmmmsssssss(z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                              z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                              a[6],a[5],a[4],a[3],a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_sub_8(ulong z[], const ulong a[])
{
    sub_ddddddddmmmmmmmmssssssss(z[7],z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                                 z[7],z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                                 a[7],a[6],a[5],a[4],a[3],a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_rsub_0(ulong FLINT_UNUSED(z[]), const ulong FLINT_UNUSED(a[]))
{
}

FLINT_FORCE_INLINE void multi_rsub_1(ulong z[], const ulong a[])
{
    z[0] = a[0] - z[0];
}

FLINT_FORCE_INLINE void multi_rsub_2(ulong z[], const ulong a[])
{
    sub_ddmmss(z[1],z[0],
               a[1],a[0],
               z[1],z[0]);
}

FLINT_FORCE_INLINE void multi_rsub_3(ulong z[], const ulong a[])
{
    sub_dddmmmsss(z[2],z[1],z[0],
                  a[2],a[1],a[0],
                  z[2],z[1],z[0]);
}

FLINT_FORCE_INLINE void multi_rsub_4(ulong z[], const ulong a[])
{
    sub_ddddmmmmssss(z[3],z[2],z[1],z[0],
                     a[3],a[2],a[1],a[0],
                     z[3],z[2],z[1],z[0]);
}

FLINT_FORCE_INLINE void multi_rsub_5(ulong z[], const ulong a[])
{
    sub_dddddmmmmmsssss(z[4],z[3],z[2],z[1],z[0],
                        a[4],a[3],a[2],a[1],a[0],
                        z[4],z[3],z[2],z[1],z[0]);
}

FLINT_FORCE_INLINE void multi_rsub_6(ulong z[], const ulong a[])
{
    sub_ddddddmmmmmmssssss(z[5],z[4],z[3],z[2],z[1],z[0],
                           a[5],a[4],a[3],a[2],a[1],a[0],
                           z[5],z[4],z[3],z[2],z[1],z[0]);
}

FLINT_FORCE_INLINE void multi_rsub_7(ulong z[], const ulong a[])
{
    sub_dddddddmmmmmmmsssssss(z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                              a[6],a[5],a[4],a[3],a[2],a[1],a[0],
                              z[6],z[5],z[4],z[3],z[2],z[1],z[0]);
}

FLINT_FORCE_INLINE void multi_rsub_8(ulong z[], const ulong a[])
{
    sub_ddddddddmmmmmmmmssssssss(z[7],z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                                 a[7],a[6],a[5],a[4],a[3],a[2],a[1],a[0],
                                 z[7],z[6],z[5],z[4],z[3],z[2],z[1],z[0]);
}

#else


#define DEFINE_IT(n) \
FLINT_FORCE_INLINE void CAT(multi_add, n)(ulong z[], const ulong a[]) \
{ \
    unsigned char cf = 0; \
    for (ulong i = 0; i < n; i++) \
        cf = _addcarry_ulong(cf, z[i], a[i], &z[i]); \
}

DEFINE_IT(0)
DEFINE_IT(1)
DEFINE_IT(2)
DEFINE_IT(3)
DEFINE_IT(4)
DEFINE_IT(5)
DEFINE_IT(6)
DEFINE_IT(7)
DEFINE_IT(8)
#undef DEFINE_IT

#define DEFINE_IT(n) \
FLINT_FORCE_INLINE void CAT(multi_sub, n)(ulong z[], const ulong a[]) \
{ \
    unsigned char cf = 0; \
    for (ulong i = 0; i < n; i++) \
        cf = _subborrow_ulong(cf, z[i], a[i], &z[i]); \
}

DEFINE_IT(0)
DEFINE_IT(1)
DEFINE_IT(2)
DEFINE_IT(3)
DEFINE_IT(4)
DEFINE_IT(5)
DEFINE_IT(6)
DEFINE_IT(7)
#undef DEFINE_IT

#endif

#ifdef __GNUC__
FLINT_FORCE_INLINE void _mul(ulong* hi, ulong* lo, ulong y, ulong x)
{
    __uint128_t p = ((__uint128_t) x) * ((__uint128_t) y);
    *lo = (ulong) (p);
    *hi = (ulong) (p >> 64);
}

FLINT_FORCE_INLINE void _madd(ulong* hi, ulong* lo, ulong y, ulong x)
{
    __uint128_t p = ((__uint128_t) *lo) | (((__uint128_t) *hi) << 64);
    p += ((__uint128_t) x) * ((__uint128_t) y);
    *lo = (ulong) (p);
    *hi = (ulong) (p >> 64);
}
#else
FLINT_FORCE_INLINE void _mul(ulong* hi, ulong* lo, ulong y, ulong x)
{
    ulong r1, r0;
    umul_ppmm(r1, r0, x, y);
    *lo = r0;
    *hi = r1;
}

FLINT_FORCE_INLINE void _madd(ulong* hi, ulong* lo, ulong y, ulong x)
{
    ulong r1, r0;
    umul_ppmm(r1, r0, x, y);
    add_ssaaaa(*hi, *lo, r1, r0, *hi, *lo);
}
#endif

/* NOTE: Define these manually to avoid compiler warnings or errors */
FLINT_FORCE_INLINE void CAT3(_big_mul, 1, 0)(ulong r[], ulong FLINT_UNUSED(t[]), ulong FLINT_UNUSED(C[]), ulong FLINT_UNUSED(y))
{
    r[0] = 0;
}

FLINT_FORCE_INLINE void CAT3(_big_addmul, 1, 0)(ulong FLINT_UNUSED(r[]), ulong FLINT_UNUSED(t[]), ulong FLINT_UNUSED(C[]), ulong FLINT_UNUSED(y))
{
}

#define DEFINE_IT(n, m) \
FLINT_FORCE_INLINE void CAT3(_big_mul, n, m)(ulong r[], ulong t[], ulong C[], ulong y) \
{ \
    for (ulong k = 0; k < n; k += 2) \
    { \
        if (k + 1 < n) \
        { \
            FLINT_ASSERT(k < m); \
            _mul(&r[k+1],&r[k+0], C[k+0], y); \
        } \
        else \
        { \
            FLINT_ASSERT(k + 1 == n); \
            if (k < m) \
                r[k+0] = C[k+0]*y; \
            else \
                r[k+0] = 0; \
        } \
 \
        if (k + 2 < n) \
        { \
            FLINT_ASSERT(k + 1 < m); \
            _mul(&t[k+2],&t[k+1], C[k+1], y); \
        } \
        else if (k + 1 < n) \
        { \
            if (k + 1 < m) \
                t[k+1] = C[k+1]*y; \
            else \
                t[k+1] = 0; \
        } \
    } \
} \
FLINT_FORCE_INLINE void CAT3(_big_addmul, n, m)(ulong r[], ulong t[], ulong C[], ulong y) \
{ \
    for (ulong k = 0; k < n; k += 2) \
    { \
        if (k + 1 < n) \
        { \
            FLINT_ASSERT(k < m); \
            _madd(&r[k+1],&r[k+0], C[k+0], y); \
        } \
        else \
        { \
            FLINT_ASSERT(k + 1 == n); \
            if (k < m) \
                r[k+0] += C[k+0]*y; \
        } \
 \
        if (k + 2 < n) \
        { \
            FLINT_ASSERT(k + 1 < m); \
            _madd(&t[k+2],&t[k+1], C[k+1], y); \
        } \
        else if (k + 1 < n) \
        { \
            if (k + 1 < m) \
                t[k+1] += C[k+1]*y; \
        } \
    } \
}

DEFINE_IT(2, 1)
DEFINE_IT(3, 2)
DEFINE_IT(4, 3)
DEFINE_IT(4, 4)
DEFINE_IT(5, 4)
DEFINE_IT(6, 5)
DEFINE_IT(7, 6)
#undef DEFINE_IT



#define DEFINE_IT(n, n_minus_1) \
FLINT_FORCE_INLINE void CAT(_reduce_big_sum, n)(ulong r[], ulong t[], const ulong* limit) \
{ \
    CAT(multi_add, n_minus_1)(r+1, t+1); \
check: \
    for (ulong k = n; k > 1; k--) \
    { \
        if (FLINT_LIKELY(r[k-1] > limit[k-1])) \
            goto sub; \
        if (r[k-1] < limit[k-1]) \
            return; \
    } \
    if (r[0] < limit[0]) \
        return; \
sub: \
    CAT(multi_sub, n)(r, limit); \
    goto check; \
}

DEFINE_IT(1, 0)
DEFINE_IT(2, 1)
DEFINE_IT(3, 2)
DEFINE_IT(4, 3)
DEFINE_IT(5, 4)
DEFINE_IT(6, 5)
DEFINE_IT(7, 6)
#undef DEFINE_IT

#ifdef __cplusplus
}
#endif

#endif /* CRT_HELPERS_H */
