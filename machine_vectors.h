/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MACHINE_VECTORS_H
#define MACHINE_VECTORS_H

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <math.h>
#include <x86intrin.h>
#include <immintrin.h>
#undef ulong
#define ulong mp_limb_t

#include "flint.h"
#include "templates.h"

#define FLINT_INLINE static __inline__
#define FLINT_FORCE_INLINE static __attribute__((always_inline)) __inline__

#define UNLIKELY(x) __builtin_expect((x),0)
#define LIKELY(x)   __builtin_expect((x),1)

#ifdef __cplusplus
 extern "C" {
#endif

/*
    In general the machine vector types should either be passed by const ref or
    the whole function should be forced inline as some platforms have buggy
    pass by value.
*/

typedef ulong vec1ui;
typedef __m128i vec2ui;
typedef __m256i vec4ui;
typedef struct {__m256i e1, e2;} vec8ui;

typedef double vec1d;
typedef __m128d vec2d;
typedef __m256d vec4d;
typedef struct {__m256d e1, e2;} vec8d;


/* the max native size for this platform */
#define NATIVE 4
#define vecNATIVEd vec4d



/* vec1 **************************************************/

FLINT_FORCE_INLINE vec1d vec1d_load(const double* a) {
    return a[0];
}

FLINT_FORCE_INLINE vec1d vec1d_load_aligned(const double* a) {
    return a[0];
}

FLINT_FORCE_INLINE vec1d vec1d_load_unaligned(const double* a) {
    return a[0];
}

FLINT_FORCE_INLINE void vec1d_store(double* z, vec1d a) {
    z[0] = a;
}

FLINT_FORCE_INLINE void vec1d_store_aligned(double* z, vec1d a) {
    z[0] = a;
}

FLINT_FORCE_INLINE void vec1d_store_unaligned(double* z, vec1d a) {
    z[0] = a;
}

FLINT_FORCE_INLINE int vec1d_same(double a, double b) {
    return a == b;
}

FLINT_FORCE_INLINE vec1d vec1d_set_d(double a) {
    return a;
}

FLINT_FORCE_INLINE vec1d vec1d_round(vec1d a) {
    return rint(a);
}

FLINT_FORCE_INLINE vec1d vec1d_zero() {
    return 0.0;
}

FLINT_FORCE_INLINE vec1d vec1d_one() {
    return 1.0;
}

FLINT_FORCE_INLINE vec1d vec1d_add(vec1d a, vec1d b) {
    return a + b;
}

FLINT_FORCE_INLINE vec1d vec1d_sub(vec1d a, vec1d b) {
    return a - b;
}

FLINT_FORCE_INLINE vec1d vec1d_addsub(vec1d a, vec1d b) {
    return a - b;
}

FLINT_FORCE_INLINE vec1d vec1d_neg(vec1d a) {
    return -a;
}

FLINT_FORCE_INLINE vec1d vec1d_abs(vec1d a) {
    return fabs(a);
}

FLINT_FORCE_INLINE vec1d vec1d_max(vec1d a, vec1d b) {
    return fmax(a, b);
}

FLINT_FORCE_INLINE vec1d vec1d_min(vec1d a, vec1d b) {
    return fmin(a, b);
}

FLINT_FORCE_INLINE vec1d vec1d_mul(vec1d a, vec1d b) {
    return a*b;
}

FLINT_FORCE_INLINE vec1d vec1d_half(vec1d a) {
    return a*0.5;
}

FLINT_FORCE_INLINE vec1d vec1d_div(vec1d a, vec1d b) {
    return a/b;
}

FLINT_FORCE_INLINE vec1d vec1d_fmadd(vec1d a, vec1d b, vec1d c) {
    return fma(a, b, c);
}

FLINT_FORCE_INLINE vec1d vec1d_fmsub(vec1d a, vec1d b, vec1d c) {
    return fma(a, b, -c);
}

FLINT_FORCE_INLINE vec1d vec1d_fnmadd(vec1d a, vec1d b, vec1d c) {
    return fma(-a, b, c);
}

FLINT_FORCE_INLINE vec1d vec1d_fnmsub(vec1d a, vec1d b, vec1d c) {
    return fma(-a, b, -c);
}

FLINT_FORCE_INLINE vec1d vec1d_blendv(vec1d a, vec1d b, vec1d c) {
    return c >= 0 ? a : b;
}

/* [0,n] -> [-n/2, n/2] */
FLINT_FORCE_INLINE vec1d vec1d_reduce_0n_to_pmhn(vec1d a, vec1d n) {
    vec1d halfn = 0.5*n;
    return a > halfn ? a - n : a;
}

FLINT_FORCE_INLINE vec1d vec1d_reduce_pm1n_to_pmhn(vec1d a, vec1d n) {
    vec1d t = a + n;
    vec1d halfn = 0.5*n;
    if (a > halfn)
        return a - n;
    else if (t < halfn)
        return t;
    else
        return a;
}

/* vec4 *****************************************************/

FLINT_FORCE_INLINE vec4d vec4d_load(const double* a) {
    return _mm256_load_pd(a);
}

FLINT_FORCE_INLINE vec4d vec4d_load_aligned(const double* a) {
    return _mm256_load_pd(a);
}

FLINT_FORCE_INLINE vec4d vec4d_load_unaligned(const double* a) {
    return _mm256_loadu_pd(a);
}

FLINT_FORCE_INLINE void vec4d_store(double* z, vec4d a) {
    _mm256_store_pd(z, a);
}

FLINT_FORCE_INLINE void vec4d_store_aligned(double* z, vec4d a) {
    _mm256_store_pd(z, a);
}

FLINT_FORCE_INLINE void vec4d_store_unaligned(double* z, vec4d a) {
    _mm256_storeu_pd(z, a);
}

FLINT_FORCE_INLINE int vec4d_same(vec4d a, vec4d b) {
    return a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3];
}

FLINT_FORCE_INLINE vec4d vec4d_set_d(double a) {
    return _mm256_set1_pd(a);
}

FLINT_FORCE_INLINE vec4d vec4d_set_d4(double a0, double a1, double a2, double a3) {
    return _mm256_set_pd(a3, a2, a1, a0);
}

FLINT_FORCE_INLINE vec4d vec4d_round(vec4d a) {
    return _mm256_round_pd(a, 4);
}

FLINT_FORCE_INLINE vec4d vec4d_zero() {
    return _mm256_setzero_pd();
}

FLINT_FORCE_INLINE vec4d vec4d_one() {
    return vec4d_set_d(1);
}

FLINT_FORCE_INLINE vec4d vec4d_add(vec4d a, vec4d b) {
    return _mm256_add_pd(a, b);
}

FLINT_FORCE_INLINE vec4d vec4d_sub(vec4d a, vec4d b) {
    return _mm256_sub_pd(a, b);
}

FLINT_FORCE_INLINE vec4d vec4d_addsub(vec4d a, vec4d b) {
    return _mm256_addsub_pd(a, b);
}

FLINT_FORCE_INLINE vec4d vec4d_neg(vec4d a) {
    __m256d mask = _mm256_castsi256_pd(_mm256_set1_epi64x(UWORD(0x8000000000000000)));
    return _mm256_xor_pd(a, mask);
}

FLINT_FORCE_INLINE vec4d vec4d_abs(vec4d a) {
    __m256d mask = _mm256_castsi256_pd(_mm256_set1_epi64x(UWORD(0x7fffffffffffffff)));
    return _mm256_and_pd(a, mask);
}

FLINT_FORCE_INLINE vec4d vec4d_max(vec4d a, vec4d b) {
    return _mm256_max_pd(a, b);
}

FLINT_FORCE_INLINE vec4d vec4d_min(vec4d a, vec4d b) {
    return _mm256_min_pd(a, b);
}

FLINT_FORCE_INLINE vec4d vec4d_mul(vec4d a, vec4d b) {
    return _mm256_mul_pd(a, b);
}

FLINT_FORCE_INLINE vec4d vec4d_half(vec4d a) {
    return vec4d_mul(a, vec4d_set_d(0.5));
}

FLINT_FORCE_INLINE vec4d vec4d_div(vec4d a, vec4d b) {
    return _mm256_div_pd(a, b);
}

FLINT_FORCE_INLINE vec4d vec4d_fmadd(vec4d a, vec4d b, vec4d c) {
#ifndef AVOID_AVX2
    return _mm256_fmadd_pd(a, b, c);
#else
    return _mm256_macc_pd(a, b, c);
#endif
}

FLINT_FORCE_INLINE vec4d vec4d_fmsub(vec4d a, vec4d b, vec4d c) {
#ifndef AVOID_AVX2
    return _mm256_fmsub_pd(a, b, c);
#else
    return _mm256_msub_pd(a, b, c);
#endif
}

FLINT_FORCE_INLINE vec4d vec4d_fnmadd(vec4d a, vec4d b, vec4d c) {
#ifndef AVOID_AVX2
    return _mm256_fnmadd_pd(a, b, c);
#else
    return _mm256_nmacc_pd(a, b, c);
#endif
}

FLINT_FORCE_INLINE vec4d vec4d_fnmsub(vec4d a, vec4d b, vec4d c) {
#ifndef AVOID_AVX2
    return _mm256_fnmsub_pd(a, b, c);
#else
    return _mm256_nmsub_pd(a, b, c);
#endif
}


FLINT_FORCE_INLINE vec4d vec4d_blendv(vec4d a, vec4d b, vec4d c) {
    return _mm256_blendv_pd(a, b, c);
}

// return {a[0], b[0], a[2], b[2]}
FLINT_FORCE_INLINE vec4d vec4d_unpacklo(vec4d a, vec4d b) {
    return _mm256_unpacklo_pd(a, b);
}

// return {a[1], b[1], a[3], b[3]}
FLINT_FORCE_INLINE vec4d vec4d_unpackhi(vec4d a, vec4d b) {
    return _mm256_unpackhi_pd(a, b);
}

/* permute_i0_i1_i2_i3(a): return {a[i0], a[i1], a[i2], a[i3]} */
#ifndef AVOID_AVX2
#define DEFINE_IT(i0, i1, i2, i3) \
FLINT_FORCE_INLINE vec4d CAT6(vec4d, permute, i0, i1, i2, i3)(vec4d a) { \
    return _mm256_permute4x64_pd(a, i0 + 4*(i1 + 4*(i2 + 4*i3))); \
}
#else
#define DEFINE_IT(i0, i1, i2, i3) \
FLINT_FORCE_INLINE vec4d CAT6(vec4d, permute, i0, i1, i2, i3)(vec4d a) { \
    return vec4d_set_d4(a[i0], a[i1], a[i2], a[i3]); \
}
#endif
DEFINE_IT(0,2,1,3)
DEFINE_IT(3,1,2,0)
DEFINE_IT(3,2,1,0)
#undef DEFINE_IT

/* permute2_i0_i1(a):  return {v[i0], v[i1]}
                    |   v[0]     |    v[1]    |    v[2]    |   v[3]     |
                      a[0], a[1]   a[2], a[3]   b[0], b[1]   b[2], b[3]
*/
#define DEFINE_IT(i0, i1) \
FLINT_FORCE_INLINE vec4d CAT4(vec4d, permute2, i0, i1)(vec4d a, vec4d b) { \
    return _mm256_permute2f128_pd(a, b, i0 + 16*i1); \
}
DEFINE_IT(0,2)
DEFINE_IT(1,3)
#undef DEFINE_IT

/* view the 4 vectors as the rows of a 4x4 matrix */
#define VEC4D_TRANSPOSE(z0, z1, z2, z3, a0, a1, a2, a3) \
{ \
    vec4d _t0, _t1, _t2, _t3; \
    _t0 = vec4d_unpacklo(a0, a1); \
    _t1 = vec4d_unpackhi(a0, a1); \
    _t2 = vec4d_unpacklo(a2, a3); \
    _t3 = vec4d_unpackhi(a2, a3); \
    z0 = vec4d_permute2_0_2(_t0, _t2); \
    z1 = vec4d_permute2_0_2(_t1, _t3); \
    z2 = vec4d_permute2_1_3(_t0, _t2); \
    z3 = vec4d_permute2_1_3(_t1, _t3); \
}

FLINT_FORCE_INLINE vec4d vec4d_cmp_ge(vec4d a, vec4d b) {
    return _mm256_cmp_pd(a, b, _CMP_GE_OQ);
}

FLINT_FORCE_INLINE vec4d vec4d_cmp_gt(vec4d a, vec4d b) {
    return _mm256_cmp_pd(a, b, _CMP_GT_OQ);
}

FLINT_FORCE_INLINE vec4d vec4d_reduce_0n_to_pmhn(vec4d a, vec4d n) {
    vec4d halfn = vec4d_half(n);
    return vec4d_blendv(a, vec4d_sub(a, n), vec4d_cmp_gt(a, halfn));
}

FLINT_FORCE_INLINE vec4d vec4d_reduce_pm1n_to_pmhn(vec4d a, vec4d n) {
    vec4d halfn = vec4d_half(n);
    vec4d t = vec4d_blendv(n, vec4d_neg(n), a);
    a = vec4d_blendv(a, vec4d_sub(a, t), vec4d_cmp_gt(vec4d_abs(a), halfn));
    return a;
}


/* vec8 **********************************************************************/

FLINT_FORCE_INLINE vec8d vec8d_set_d(double a) {
    vec4d z1 = vec4d_set_d(a);
    vec8d z = {z1, z1};
    return z;
}

FLINT_FORCE_INLINE vec8d vec8d_set_d8(double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7) {
    vec4d z1 = vec4d_set_d4(a0, a1, a2, a3);
    vec4d z2 = vec4d_set_d4(a4, a5, a6, a7);
    vec8d z = {z1, z2};
    return z;
}

FLINT_FORCE_INLINE vec8d vec8d_load(const double* a) {
    vec8d z = {vec4d_load(a+0), vec4d_load(a+4)};
    return z;
}

FLINT_FORCE_INLINE vec8d vec8d_load_aligned(const double* a) {
    vec8d z = {vec4d_load_aligned(a+0), vec4d_load_aligned(a+4)};
    return z;
}

FLINT_FORCE_INLINE vec8d vec8d_load_unaligned(const double* a) {
    vec8d z = {vec4d_load_unaligned(a+0), vec4d_load_unaligned(a+4)};
    return z;
}

FLINT_FORCE_INLINE void vec8d_store(double* z, vec8d a) {
    vec4d_store(z+0, a.e1);
    vec4d_store(z+4, a.e2);
}

FLINT_FORCE_INLINE void vec8d_store_aligned(double* z, vec8d a) {
    vec4d_store_aligned(z+0, a.e1);
    vec4d_store_aligned(z+4, a.e2);
}

FLINT_FORCE_INLINE void vec8d_store_unaligned(double* z, vec8d a) {
    vec4d_store_unaligned(z+0, a.e1);
    vec4d_store_unaligned(z+4, a.e2);
}

FLINT_FORCE_INLINE int vec8d_same(vec8d a, vec8d b) {
    return vec4d_same(a.e1, b.e1) && vec4d_same(a.e2, b.e2);
}


/* reduce_pm1no_to_0n(a, n): return a mod n in [0,n) assuming a in (-n,n) */
#define DEFINE_IT(V) \
FLINT_FORCE_INLINE V V##_reduce_pm1no_to_0n(V a, V n) { \
    return V##_blendv(a, V##_add(a, n), a); \
}
DEFINE_IT(vec1d)
DEFINE_IT(vec4d)
#undef DEFINE_IT

/* reduce_to_pm1n(a, n, ninv): return a mod n in [-n,n] */
#define DEFINE_IT(V) \
FLINT_FORCE_INLINE V V##_reduce_to_pm1n(V a, V n, V ninv) { \
    return V##_fnmadd(V##_round(V##_mul(a, ninv)), n, a); \
}
DEFINE_IT(vec1d)
DEFINE_IT(vec4d)
#undef DEFINE_IT

/* reduce_to_pm1n(a, n, ninv): return a mod n in (-n,n) */
#define DEFINE_IT(V) \
FLINT_FORCE_INLINE V V##_reduce_to_pm1no(V a, V n, V ninv) { \
    return V##_fnmadd(V##_round(V##_mul(a, ninv)), n, a); \
}
DEFINE_IT(vec1d)
DEFINE_IT(vec4d)
#undef DEFINE_IT


/* reduce_to_0n(a, n, ninv): return a mod n in [0,n) */
#define DEFINE_IT(V) \
FLINT_FORCE_INLINE V V##_reduce_to_0n(V a, V n, V ninv) { \
    return V##_reduce_pm1no_to_0n(V##_reduce_to_pm1no(a, n, ninv), n); \
}
DEFINE_IT(vec1d)
DEFINE_IT(vec4d)
#undef DEFINE_IT

#define DEFINE_IT(V) \
FLINT_FORCE_INLINE int V##_same_mod(V a, V b, V n, V ninv) { \
    return V##_same(V##_reduce_to_0n(a, n, ninv), V##_reduce_to_0n(b, n, ninv)); \
}
DEFINE_IT(vec1d)
DEFINE_IT(vec4d)
#undef DEFINE_IT

/* mulmod2(a, b, n, ninv): return a*b mod n in [-n,n] with assumptions */
#define DEFINE_IT(V) \
FLINT_FORCE_INLINE V V##_mulmod2(V a, V b, V n, V ninv) { \
    V h = V##_mul(a, b); \
    V q = V##_round(V##_mul(h, ninv)); \
    V l = V##_fmsub(a, b, h); \
    return V##_add(V##_fnmadd(q, n, h), l); \
}
DEFINE_IT(vec1d)
DEFINE_IT(vec4d)
#undef DEFINE_IT




#define EXTEND_VEC_DEF0(U, V, f) \
FLINT_FORCE_INLINE V V##f() { \
    U z1 = U##f(); \
    U z2 = U##f(); \
    V z = {z1, z2}; \
    return z; \
}

#define EXTEND_VEC_DEF1(U, V, f) \
FLINT_FORCE_INLINE V V##f(V a) { \
    U z1 = U##f(a.e1); \
    U z2 = U##f(a.e2); \
    V z = {z1, z2}; \
    return z; \
}

#define EXTEND_VEC_DEF2(U, V, f) \
FLINT_FORCE_INLINE V V##f(V a, V b) { \
    U z1 = U##f(a.e1, b.e1); \
    U z2 = U##f(a.e2, b.e2); \
    V z = {z1, z2}; \
    return z; \
}

#define EXTEND_VEC_DEF3(U, V, f) \
FLINT_FORCE_INLINE V V##f(V a, V b, V c) { \
    U z1 = U##f(a.e1, b.e1, c.e1); \
    U z2 = U##f(a.e2, b.e2, c.e2); \
    V z = {z1, z2}; \
    return z; \
}

#define EXTEND_VEC_DEF4(U, V, f) \
FLINT_FORCE_INLINE V V##f(V a, V b, V c, V d) { \
    U z1 = U##f(a.e1, b.e1, c.e1, d.e1); \
    U z2 = U##f(a.e2, b.e2, c.e2, d.e2); \
    V z = {z1, z2}; \
    return z; \
}

EXTEND_VEC_DEF0(vec4d, vec8d, _zero)
EXTEND_VEC_DEF1(vec4d, vec8d, _neg)
EXTEND_VEC_DEF1(vec4d, vec8d, _round)
EXTEND_VEC_DEF2(vec4d, vec8d, _add)
EXTEND_VEC_DEF2(vec4d, vec8d, _sub)
EXTEND_VEC_DEF2(vec4d, vec8d, _min)
EXTEND_VEC_DEF2(vec4d, vec8d, _max)
EXTEND_VEC_DEF2(vec4d, vec8d, _mul)
EXTEND_VEC_DEF2(vec4d, vec8d, _div)
EXTEND_VEC_DEF2(vec4d, vec8d, _reduce_pm1n_to_pmhn)
EXTEND_VEC_DEF2(vec4d, vec8d, _reduce_pm1no_to_0n)
EXTEND_VEC_DEF3(vec4d, vec8d, _reduce_to_pm1n)
EXTEND_VEC_DEF3(vec4d, vec8d, _reduce_to_pm1no)
EXTEND_VEC_DEF3(vec4d, vec8d, _reduce_to_0n)
EXTEND_VEC_DEF3(vec4d, vec8d, _fmadd)
EXTEND_VEC_DEF3(vec4d, vec8d, _fmsub)
EXTEND_VEC_DEF3(vec4d, vec8d, _fnmadd)
EXTEND_VEC_DEF3(vec4d, vec8d, _fnmsub)
EXTEND_VEC_DEF3(vec4d, vec8d, _blendv)
EXTEND_VEC_DEF4(vec4d, vec8d, _mulmod2)

#undef EXTEND_VEC_DEF4
#undef EXTEND_VEC_DEF3
#undef EXTEND_VEC_DEF2
#undef EXTEND_VEC_DEF1
#undef EXTEND_VEC_DEF0


#ifdef __cplusplus
}
#endif

#endif

