/*
    Copyright (C) 2022 Daniel Schultz
    Copyright (C) 2023 Mathieu Gouttenoire

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MACHINE_VECTORS_H
#define MACHINE_VECTORS_H

#define ALIGN_STRUCT(x) __attribute__((aligned(x)))

#include <math.h>

#if defined(__GNUC__)
# if defined(__AVX2__)
#  include <immintrin.h>
# elif defined(__ARM_NEON)
#  include <arm_neon.h>
# endif
#elif defined(_MSC_VER)
# if defined(__AVX2__)
#  include <intrin.h>
# elif defined(_M_ARM64)
#  include <arm_neon.h>
# endif
#endif

#include "flint.h"
#include "templates.h"

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__AVX2__)
/*
    In general the machine vector types should either be passed by const ref or
    the whole function should be forced inline as some platforms have buggy
    pass by value.
*/

typedef ulong vec1n;
typedef __m128i vec2n;
typedef __m256i vec4n;
typedef struct {__m256i e1, e2;} vec8n;

typedef double vec1d;
typedef __m128d vec2d;
typedef __m256d vec4d;
typedef struct {__m256d e1, e2;} vec8d;



/* vec1d -- AVX2 ***********************************************************/

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

FLINT_FORCE_INLINE vec1d vec1d_zero(void) {
    return 0.0;
}

FLINT_FORCE_INLINE vec1d vec1d_one(void) {
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

/* [0,2n) to [0,n) */
FLINT_FORCE_INLINE vec1d vec1d_reduce_2n_to_n(vec1d a, vec1d n) {
    return a-n >= 0 ? a-n : a;
}


/* vec4d -- AVX2 ***********************************************************/

FLINT_FORCE_INLINE void vec4d_print(vec4d a)
{
#ifdef _MSC_VER
    double as[4];
    _mm256_storeu_pd(as, a);
    flint_printf("{%f, %f, %f, %f}", as[0], as[1], as[2], as[3]);
#else
    flint_printf("{%f, %f, %f, %f}", a[0], a[1], a[2], a[3]);
#endif
}

FLINT_FORCE_INLINE double vec4d_get_index(vec4d a, const int i) {
#ifdef _MSC_VER
    double as[4];
    _mm256_storeu_pd(as, a);
    return as[i];
#else
    return a[i];
#endif
}

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
#ifdef _MSC_VER
    double as[4], bs[4];
    _mm256_storeu_pd(as, a);
    _mm256_storeu_pd(bs, b);
    return as[0] == bs[0] && as[1] == bs[1] && as[2] == bs[2] && as[3] == bs[3];
#else
    return a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3];
#endif
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

FLINT_FORCE_INLINE vec4d vec4d_zero(void) {
    return _mm256_setzero_pd();
}

FLINT_FORCE_INLINE vec4d vec4d_one(void) {
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

FLINT_FORCE_INLINE vec4d vec4d_unpack_lo_permute_0_2_1_3(vec4d u, vec4d v) {
    return vec4d_permute_0_2_1_3(vec4d_unpacklo(u, v));
}

FLINT_FORCE_INLINE vec4d vec4d_unpack_hi_permute_0_2_1_3(vec4d u, vec4d v) {
    return vec4d_permute_0_2_1_3(vec4d_unpackhi(u, v));
}

FLINT_FORCE_INLINE vec4d vec4d_unpackhi_permute_3_1_2_0(vec4d u, vec4d v) {
    return vec4d_permute_3_1_2_0(vec4d_unpackhi(u, v));
}

FLINT_FORCE_INLINE vec4d vec4d_unpacklo_permute_3_1_2_0(vec4d u, vec4d v) {
    return vec4d_permute_3_1_2_0(vec4d_unpacklo(u, v));
}

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

/* [0,2n) to [0,n) */
FLINT_FORCE_INLINE vec4d vec4d_reduce_2n_to_n(vec4d a, vec4d n) {
    vec4d s = vec4d_sub(a, n);
    return vec4d_blendv(s, a, s);
}

FLINT_FORCE_INLINE vec4n vec4d_convert_limited_vec4n(vec4d a) {
    __m256d t = _mm256_set1_pd(0x1.0p52);
    return _mm256_castpd_si256(_mm256_xor_pd(_mm256_add_pd(a, t), t));
}

/* vec4n -- AVX2 ***********************************************************/

FLINT_FORCE_INLINE void vec4n_print(vec4n a)
{
    flint_printf("[hi %016llx_%016llx_%016llx_%016llx lo]",
        _mm256_extract_epi64(a, 3),
        _mm256_extract_epi64(a, 2),
        _mm256_extract_epi64(a, 1),
        _mm256_extract_epi64(a, 0));
}

FLINT_FORCE_INLINE vec4n vec4n_load_unaligned(const ulong* a) {
    return _mm256_loadu_si256((__m256i*) a);
}

FLINT_FORCE_INLINE void vec4n_store_unaligned(ulong* z, vec4n a) {
    _mm256_storeu_si256((__m256i*) z, a);
}

/* permute_i0_i1_i2_i3(a): return {a[i0], a[i1], a[i2], a[i3]} */
#ifndef AVOID_AVX2
#define DEFINE_IT(i0, i1, i2, i3)                                        \
FLINT_FORCE_INLINE vec4n CAT6(vec4n, permute, i0, i1, i2, i3)(vec4n a) { \
    return _mm256_permute4x64_epi64(a, i0 + 4*(i1 + 4*(i2 + 4*i3)));     \
}
#else
#define DEFINE_IT(i0, i1, i2, i3)                                        \
FLINT_FORCE_INLINE vec4n CAT6(vec4n, permute, i0, i1, i2, i3)(vec4n a) { \
    return vec4n_set_n4(a[i0], a[i1], a[i2], a[i3]);                     \
}
#endif
DEFINE_IT(3,2,1,0)
#undef DEFINE_IT

FLINT_FORCE_INLINE vec4n vec4n_zero()
{
    return _mm256_setzero_si256();
}

FLINT_FORCE_INLINE vec4n vec4n_set_n4(ulong a0, ulong a1, ulong a2, ulong a3) {
    return _mm256_set_epi64x(a3, a2, a1, a0);
}

FLINT_FORCE_INLINE vec4n vec4n_set_n(ulong a) {
  return _mm256_set1_epi64x(a);
}

FLINT_FORCE_INLINE vec4n vec4n_add(vec4n a, vec4n b)
{
    return _mm256_add_epi64(a, b);
}

FLINT_FORCE_INLINE vec4n vec4n_sub(vec4n a, vec4n b)
{
    return _mm256_sub_epi64(a, b);
}

/* for n < 2^63 */
FLINT_FORCE_INLINE vec4n vec4n_addmod_limited(vec4n a, vec4n b, vec4n n)
{
    vec4n s = vec4n_add(a, b);
    vec4n t = vec4n_sub(s, n);
    vec4n m = _mm256_srai_epi32(t, 31);
          m = _mm256_shuffle_epi32(m, 1 + 4*(1 + 4*(3 + 4*(3))));
    return _mm256_blendv_epi8(t, s, m);
}

FLINT_FORCE_INLINE vec4n vec4n_addmod(vec4n a, vec4n b, vec4n n)
{
    vec4n tt = vec4n_set_n(0x8000000000000000);
    vec4n s = vec4n_add(a, b);
#if 0
    vec4n  m = vec4n_sub(n, tt);
    vec4n t0 = vec4n_sub(m, a);
    vec4n t1 = vec4n_sub(b, tt);
    vec4n t2 = vec4n_sub(t1, t0);
    return _mm256_blendv_epi8(s, t2, _mm256_cmpgt_epi64(t1, t0));
#else
    vec4n t0 = vec4n_sub(s, n);
    vec4n t1 = vec4n_sub(a, tt);
    vec4n t2 = vec4n_sub(t0, tt);
    return _mm256_blendv_epi8(t0, s, _mm256_cmpgt_epi64(t2, t1));
#endif
}

FLINT_FORCE_INLINE vec4n vec4n_mul(vec4n u, vec4n v)
{
    return _mm256_mul_epu32(u, v);
}

// compilers accept non-(compile-time-)constant b with srli
// (still, faster if constant)
FLINT_FORCE_INLINE vec4n vec4n_bit_shift_right(vec4n a, ulong b) {
    return _mm256_srli_epi64(a, b);
}

#define vec4n_bit_shift_right_32(a) vec4n_bit_shift_right((a), 32)

FLINT_FORCE_INLINE vec4n vec4n_bit_and(vec4n a, vec4n b) {
    return _mm256_and_si256(a, b);
}

FLINT_FORCE_INLINE vec4d vec4n_convert_limited_vec4d(vec4n a) {
    __m256d t = _mm256_set1_pd(0x1.0p52);
    return _mm256_sub_pd(_mm256_or_pd(_mm256_castsi256_pd(a), t), t);
}

// horizontal sum
FLINT_FORCE_INLINE ulong vec4n_horizontal_sum(vec4n a) {
    vec4n a_hi = _mm256_shuffle_epi32(a, 14);  // 14 == 0b00001110
    vec4n sum_lo = _mm256_add_epi64(a, a_hi);
    vec2n sum_hi = _mm256_extracti128_si256(sum_lo, 1);
    vec2n sum = _mm_add_epi64(_mm256_castsi256_si128(sum_lo), sum_hi);
    return (ulong) _mm_cvtsi128_si64(sum);
}


/* vec8d -- AVX2 ***********************************************************/

FLINT_FORCE_INLINE double vec8d_get_index(vec8d a, int i) {
    return i < 4 ? vec4d_get_index(a.e1, i) : vec4d_get_index(a.e2, i - 4);
}

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

FLINT_FORCE_INLINE vec8d vec8n_convert_limited_vec8d(vec8n a) {
    vec8d z = {vec4n_convert_limited_vec4d(a.e1), vec4n_convert_limited_vec4d(a.e2)};
    return z;
}

/* !!! the outputs are also permuted !!! */
FLINT_FORCE_INLINE vec8d _vec8i32_convert_vec8d(__m256i a)
{
    __m256i mask = _mm256_set1_epi32(0x43300000);
    __m256i ak0 = _mm256_unpacklo_epi32(a, mask);
    __m256i ak1 = _mm256_unpackhi_epi32(a, mask);
    __m256d t = _mm256_set1_pd(0x1.0p52);
    vec8d z;
    z.e1 = _mm256_sub_pd(_mm256_castsi256_pd(ak0), t);
    z.e2 = _mm256_sub_pd(_mm256_castsi256_pd(ak1), t);
    return z;
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

/*
    mulmod(a, b, n, 1.0/n): return a*b mod n
    assuming fft_small_mulmod_satisfies_bounds(n):
     Products a*b in the range (-2*n^2, 2*n^2) are reduced to the range (-1*n, 1*n).
     Products a*b in the range (-4*n^2, 4*n^2) are reduced to the range (-3/2*n, 3/2*n).
    otherwise the more pessimistic:
     Products a*b in the range (-2*n^2, 2*n^2) are reduced to the range (-9/8*n, 9/8*n).
     Products a*b in the range (-4*n^2, 4*n^2) are reduced to the range (-7/4*n, 7/4*n).
    floating point round mode is always assumed to be set to nearest
*/
#define DEFINE_IT(V) \
FLINT_FORCE_INLINE V V##_mulmod(V a, V b, V n, V ninv) { \
    V h = V##_mul(a, b); \
    V q = V##_round(V##_mul(h, ninv)); \
    V l = V##_fmsub(a, b, h); \
    return V##_add(V##_fnmadd(q, n, h), l); \
} \
 \
FLINT_FORCE_INLINE V V##_nmulmod(V a, V b, V n, V ninv) { \
    V h = V##_mul(a, b); \
    V q = V##_round(V##_mul(h, ninv)); \
    V l = V##_fnmadd(a, b, h); \
    return V##_sub(l, V##_fnmadd(q, n, h)); \
}

DEFINE_IT(vec1d)
DEFINE_IT(vec4d)
#undef DEFINE_IT




#define EXTEND_VEC_DEF0(U, V, f) \
FLINT_FORCE_INLINE V V##f(void) { \
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
EXTEND_VEC_DEF2(vec4d, vec8d, _reduce_2n_to_n)
EXTEND_VEC_DEF3(vec4d, vec8d, _reduce_to_pm1n)
EXTEND_VEC_DEF3(vec4d, vec8d, _reduce_to_pm1no)
EXTEND_VEC_DEF3(vec4d, vec8d, _reduce_to_0n)
EXTEND_VEC_DEF3(vec4d, vec8d, _fmadd)
EXTEND_VEC_DEF3(vec4d, vec8d, _fmsub)
EXTEND_VEC_DEF3(vec4d, vec8d, _fnmadd)
EXTEND_VEC_DEF3(vec4d, vec8d, _fnmsub)
EXTEND_VEC_DEF3(vec4d, vec8d, _blendv)
EXTEND_VEC_DEF3(vec4n, vec8n, _addmod)
EXTEND_VEC_DEF3(vec4n, vec8n, _addmod_limited)
EXTEND_VEC_DEF4(vec4d, vec8d, _mulmod)
EXTEND_VEC_DEF4(vec4d, vec8d, _nmulmod)

#undef EXTEND_VEC_DEF4
#undef EXTEND_VEC_DEF3
#undef EXTEND_VEC_DEF2
#undef EXTEND_VEC_DEF1
#undef EXTEND_VEC_DEF0



/* vec8n -- AVX2 ***********************************************************/

FLINT_FORCE_INLINE vec8n vec8n_load_unaligned(const ulong* a) {
    vec8n z = {vec4n_load_unaligned(a+0), vec4n_load_unaligned(a+4)};
    return z;
}

FLINT_FORCE_INLINE vec8n vec8n_set_n(ulong a) {
    vec4n x = vec4n_set_n(a);
    vec8n z = {x, x};
    return z;
}

FLINT_FORCE_INLINE vec8n vec8n_bit_shift_right(vec8n a, ulong b) {
    vec8n z = {vec4n_bit_shift_right(a.e1, b), vec4n_bit_shift_right(a.e2, b)};
    return z;
}

#define vec8n_bit_shift_right_32(a) vec8n_bit_shift_right((a), 32)

FLINT_FORCE_INLINE vec8n vec8n_bit_and(vec8n a, vec8n b) {
    vec8n z = {vec4n_bit_and(a.e1, b.e1), vec4n_bit_and(a.e2, b.e2)};
    return z;
}



#elif defined(__ARM_NEON) || defined(_M_ARM64)

typedef ulong vec1n;
typedef uint64x2_t vec2n;
typedef struct {vec2n e1, e2;} vec4n;
typedef struct {vec4n e1, e2;} vec8n;

typedef double vec1d;
typedef float64x2_t vec2d;
typedef struct {vec2d e1, e2;} vec4d;
typedef struct {vec4d e1, e2;} vec8d;

#define EXTEND_VEC_DEF0(U, V, f) \
FLINT_FORCE_INLINE V V##f(void) { \
    U z1 = U##f(); \
    U z2 = U##f(); \
    V z = {z1, z2}; return z; \
}

#define EXTEND_VEC_DEF1(U, V, f) \
FLINT_FORCE_INLINE V V##f(V a) { \
    U z1 = U##f(a.e1); \
    U z2 = U##f(a.e2); \
    V z = {z1, z2}; return z; \
}

#define EXTEND_VEC_DEF2(U, V, f) \
FLINT_FORCE_INLINE V V##f(V a, V b) { \
    U z1 = U##f(a.e1, b.e1); \
    U z2 = U##f(a.e2, b.e2); \
    V z = {z1, z2}; return z; \
}

#define EXTEND_VEC_DEF3(U, V, f) \
FLINT_FORCE_INLINE V V##f(V a, V b, V c) { \
    U z1 = U##f(a.e1, b.e1, c.e1); \
    U z2 = U##f(a.e2, b.e2, c.e2); \
    V z = {z1, z2}; return z; \
}

#define EXTEND_VEC_DEF4(U, V, f) \
FLINT_FORCE_INLINE V V##f(V a, V b, V c, V d) { \
    U z1 = U##f(a.e1, b.e1, c.e1, d.e1); \
    U z2 = U##f(a.e2, b.e2, c.e2, d.e2); \
    V z = {z1, z2}; return z; \
}


/* floating point stuff ******************************************************/

/* vec1d -- NEON/ARM64 *********************************************/

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

FLINT_FORCE_INLINE vec1d vec1d_zero(void) {
    return 0.0;
}

FLINT_FORCE_INLINE vec1d vec1d_one(void) {
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

/* vec2d -- NEON/ARM64 *********************************************/

FLINT_FORCE_INLINE double vec2d_get_index(vec2d a, int i) {
    return a[i];
}

FLINT_FORCE_INLINE vec2d vec2d_set_d2(double a0, double a1)
{
    double ALIGN_STRUCT(16) data[2] = {a0, a1};
    return vld1q_f64((float64_t *) data);
}

FLINT_FORCE_INLINE vec2d vec2d_set_d(double a)
{
    return vdupq_n_f64(a);
}

FLINT_FORCE_INLINE vec2d vec2d_load(const double* a)
{
    return vld1q_f64(a);
}

FLINT_FORCE_INLINE vec2d vec2d_load_unaligned(const double* a)
{
    return vld1q_f64(a);
}

FLINT_FORCE_INLINE vec2d vec2d_load_aligned(const double* a)
{
    return vld1q_f64(a);
}

FLINT_FORCE_INLINE void vec2d_store(double* z, vec2d a)
{
    vst1q_f64(z, a);
}

FLINT_FORCE_INLINE void vec2d_store_unaligned(double* z, vec2d a)
{
    vst1q_f64(z, a);
}

FLINT_FORCE_INLINE void vec2d_store_aligned(double* z, vec2d a)
{
    vst1q_f64(z, a);
}

FLINT_FORCE_INLINE vec2d vec2d_zero(void)
{
    return vdupq_n_f64(0.0);
}

FLINT_FORCE_INLINE vec2d vec2d_one(void)
{
    return vdupq_n_f64(1.0);
}

FLINT_FORCE_INLINE vec2d vec2d_neg(vec2d a)
{
    return vnegq_f64(a);
}

FLINT_FORCE_INLINE vec2d vec2d_abs(vec2d a)
{
    return vabsq_f64(a);
}

FLINT_FORCE_INLINE vec2d vec2d_min(vec2d a, vec2d b)
{
    return vminq_f64(a, b);
}

FLINT_FORCE_INLINE vec2d vec2d_max(vec2d a, vec2d b)
{
    return vmaxq_f64(a, b);
}

FLINT_FORCE_INLINE vec2d vec2d_cmp_lt(vec2d a, vec2d b) {
    return vcvtq_f64_u64(vcltq_f64(a, b));
}

FLINT_FORCE_INLINE vec2d vec2d_cmp_gt(vec2d a, vec2d b) {
    return vcvtq_f64_u64(vcgtq_f64(a, b));
}

FLINT_FORCE_INLINE vec2d vec2d_blendv(vec2d a, vec2d b, vec2d c) {
    return vbslq_f64(vcvtq_u64_f64(c), b, a);
}

FLINT_FORCE_INLINE vec2d vec2d_add(vec2d a, vec2d b)
{
    return vaddq_f64(a, b);
}

FLINT_FORCE_INLINE vec2d vec2d_sub(vec2d a, vec2d b)
{
    return vsubq_f64(a, b);
}

FLINT_FORCE_INLINE vec2d vec2d_mul(vec2d a, vec2d b)
{
    return vmulq_f64(a, b);
}

FLINT_FORCE_INLINE vec2d vec2d_half(vec2d a)
{
    return vmulq_n_f64(a, 0.5);
}

FLINT_FORCE_INLINE vec2d vec2d_fmadd(vec2d a, vec2d b, vec2d c) {
    return vfmaq_f64(c, a, b);
}

FLINT_FORCE_INLINE vec2d vec2d_fnmadd(vec2d a, vec2d b, vec2d c) {
    return vfmsq_f64(c, a, b);
}

/* these two might be slow due to the extra negation */
FLINT_FORCE_INLINE vec2d vec2d_fmsub(vec2d a, vec2d b, vec2d c) {
    return -vfmsq_f64(c, a, b);
}

FLINT_FORCE_INLINE vec2d vec2d_fnmsub(vec2d a, vec2d b, vec2d c) {
    return -vfmaq_f64(c, a, b);
}

FLINT_FORCE_INLINE vec2d vec2d_div(vec2d a, vec2d b)
{
    return vdivq_f64(a, b);
}

FLINT_FORCE_INLINE vec2d vec2d_round(vec2d a)
{
    return vrndnq_f64(a);
}

FLINT_FORCE_INLINE vec2d vec2d_unpacklo(vec2d a, vec2d b) {
    return vtrn1q_f64(a, b);
}

FLINT_FORCE_INLINE vec2d vec2d_unpackhi(vec2d a, vec2d b) {
    return vtrn2q_f64(a, b);
}

FLINT_FORCE_INLINE vec1d vec1d_reduce_to_pm1no(vec1d a, vec1d n, vec1d ninv) {
    return vec1d_fnmadd(vec1d_round(vec1d_mul(a, ninv)), n, a);
}

FLINT_FORCE_INLINE vec2d vec2d_reduce_to_pm1no(vec2d a, vec2d n, vec2d ninv) {
    return vec2d_fnmadd(vec2d_round(vec2d_mul(a, ninv)), n, a);
}

FLINT_FORCE_INLINE vec1d vec1d_reduce_to_pm1n(vec1d a, vec1d n, vec1d ninv) {
    return vec1d_fnmadd(vec1d_round(vec1d_mul(a, ninv)), n, a);
}

FLINT_FORCE_INLINE vec2d vec2d_reduce_to_pm1n(vec2d a, vec2d n, vec2d ninv) {
    return vec2d_fnmadd(vec2d_round(vec2d_mul(a, ninv)), n, a);
}

FLINT_FORCE_INLINE vec2d vec2d_reduce_0n_to_pmhn(vec2d a, vec2d n) {
    vec2d halfn = vec2d_half(n);
    return vbslq_f64(vcgtq_f64(a, halfn), vec2d_sub(a, n), a);
}

FLINT_FORCE_INLINE vec2d vec2d_reduce_pm1n_to_pmhn(vec2d a, vec2d n) {
    vec2d halfn = vec2d_half(n);
    vec2d t = vec2d_add(a, n);

    vec2n condition_a = vcgtq_f64(a, halfn);
    vec2n condition_t = vcltq_f64(t, halfn);

    return vbslq_f64(condition_a, vec2d_sub(a, n), vbslq_f64(condition_t, t, a));
}

FLINT_FORCE_INLINE vec1d vec1d_reduce_pm1no_to_0n(vec1d a, vec1d n) {
    return a >= 0 ? a : a + n;
}

FLINT_FORCE_INLINE vec2d vec2d_reduce_pm1no_to_0n(vec2d a, vec2d n) {
    return vbslq_f64(vcgeq_f64(a, vec2d_zero()), a, vaddq_f64(a, n));
}

FLINT_FORCE_INLINE vec1d vec1d_reduce_to_0n(vec1d a, vec1d n, vec1d ninv) {
    return vec1d_reduce_pm1no_to_0n(vec1d_reduce_to_pm1no(a, n, ninv), n);
}

FLINT_FORCE_INLINE vec2d vec2d_reduce_to_0n(vec2d a, vec2d n, vec2d ninv) {
    return vec2d_reduce_pm1no_to_0n(vec2d_reduce_to_pm1no(a, n, ninv), n);
}


/* mulmod(a, b, n, ninv): return a*b mod n in [-n,n] with assumptions */
#define DEFINE_IT(V) \
FLINT_FORCE_INLINE V V##_mulmod(V a, V b, V n, V ninv) { \
    V h = V##_mul(a, b); \
    V q = V##_round(V##_mul(h, ninv)); \
    V l = V##_fnmadd(a, b, h); \
    return V##_sub(V##_fnmadd(q, n, h), l); \
} \
 \
FLINT_FORCE_INLINE V V##_nmulmod(V a, V b, V n, V ninv) { \
    V h = V##_mul(a, b); \
    V q = V##_round(V##_mul(h, ninv)); \
    V l = V##_fnmadd(a, b, h); \
    return V##_sub(l, V##_fnmadd(q, n, h)); \
}

DEFINE_IT(vec1d)
DEFINE_IT(vec2d)
#undef DEFINE_IT



/* vec4d -- NEON/ARM64 *********************************************/

FLINT_FORCE_INLINE vec4d vec4d_set_vec2d2(vec2d a, vec2d b) {
    vec4d z = {a, b}; return z;
}

FLINT_FORCE_INLINE double vec4d_get_index(vec4d a, int i) {
    return i < 2 ? vec2d_get_index(a.e1, i) : vec2d_get_index(a.e2, i - 2);
}

FLINT_FORCE_INLINE vec4d vec4d_set_d4(double a0, double a1, double a2, double a3)
{
    return vec4d_set_vec2d2(vec2d_set_d2(a0, a1), vec2d_set_d2(a2, a3));
}

FLINT_FORCE_INLINE vec4d vec4d_set_d(double a)
{
    vec2d z1 = vec2d_set_d(a);
    vec4d z = {z1, z1}; return z;
}

FLINT_FORCE_INLINE void vec4d_store(double* z, vec4d a)
{
    vec2d_store(z+0, a.e1);
    vec2d_store(z+2, a.e2);
}

FLINT_FORCE_INLINE void vec4d_store_aligned(double* z, vec4d a)
{
    vec2d_store_aligned(z+0, a.e1);
    vec2d_store_aligned(z+2, a.e2);
}

FLINT_FORCE_INLINE void vec4d_store_unaligned(double* z, vec4d a)
{
    vec2d_store_unaligned(z+0, a.e1);
    vec2d_store_unaligned(z+2, a.e2);
}

FLINT_FORCE_INLINE vec4d vec4d_load(const double* a)
{
    return vec4d_set_vec2d2(vec2d_load(a+0),
                            vec2d_load(a+2));
}

FLINT_FORCE_INLINE vec4d vec4d_load_aligned(const double* a)
{
    return vec4d_set_vec2d2(vec2d_load_aligned(a+0),
                            vec2d_load_aligned(a+2));
}

FLINT_FORCE_INLINE vec4d vec4d_load_unaligned(const double* a)
{
    return vec4d_set_vec2d2(vec2d_load_unaligned(a+0),
                            vec2d_load_unaligned(a+2));
}

FLINT_FORCE_INLINE vec4d vec4d_permute_0_2_1_3(vec4d a)
{
    return vec4d_set_vec2d2(vec2d_unpacklo(a.e1, a.e2),
                            vec2d_unpackhi(a.e1, a.e2));
}

FLINT_FORCE_INLINE vec4d vec4d_permute_3_2_1_0(vec4d a) {
    return vec4d_set_vec2d2(vextq_f64(a.e2, a.e2, 1),
                            vextq_f64(a.e1, a.e1, 1));
}

FLINT_FORCE_INLINE vec4d vec4d_permute_3_1_2_0(vec4d a) {
    vec2d z1 = vzip2q_f64(a.e2, a.e1);
    vec2d z2 = vzip1q_f64(a.e2, a.e1);
    return vec4d_set_vec2d2(z1, z2);
}

FLINT_FORCE_INLINE vec4d vec4d_permute2_0_2(vec4d a, vec4d b) {
    return vec4d_set_vec2d2(a.e1, b.e1);
}

FLINT_FORCE_INLINE vec4d vec4d_permute2_1_3(vec4d a, vec4d b) {
    return vec4d_set_vec2d2(a.e2, b.e2);
}

FLINT_FORCE_INLINE vec4d vec4d_unpack_lo_permute_0_2_1_3(vec4d u, vec4d v) {
    return vec4d_set_vec2d2(vec2d_unpacklo(u.e1, u.e2),
                            vec2d_unpacklo(v.e1, v.e2));
}

FLINT_FORCE_INLINE vec4d vec4d_unpack_hi_permute_0_2_1_3(vec4d u, vec4d v) {
    return vec4d_set_vec2d2(vec2d_unpackhi(u.e1, u.e2),
                            vec2d_unpackhi(v.e1, v.e2));
}

FLINT_FORCE_INLINE vec4d vec4d_unpackhi_permute_3_1_2_0(vec4d u, vec4d v) {
    return vec4d_set_vec2d2(vec2d_unpackhi(v.e2, v.e1),
                            vec2d_unpackhi(u.e2, u.e1));
}

FLINT_FORCE_INLINE vec4d vec4d_unpacklo_permute_3_1_2_0(vec4d u, vec4d v) {
    return vec4d_set_vec2d2(vec2d_unpacklo(v.e2, v.e1),
                            vec2d_unpacklo(u.e2, u.e1));
}

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

EXTEND_VEC_DEF0(vec2d, vec4d, _zero)
EXTEND_VEC_DEF0(vec2d, vec4d, _one)
EXTEND_VEC_DEF1(vec2d, vec4d, _neg)
EXTEND_VEC_DEF1(vec2d, vec4d, _round)
EXTEND_VEC_DEF2(vec2d, vec4d, _add)
EXTEND_VEC_DEF2(vec2d, vec4d, _sub)
EXTEND_VEC_DEF2(vec2d, vec4d, _min)
EXTEND_VEC_DEF2(vec2d, vec4d, _max)
EXTEND_VEC_DEF2(vec2d, vec4d, _mul)
EXTEND_VEC_DEF2(vec2d, vec4d, _div)
EXTEND_VEC_DEF2(vec2d, vec4d, _reduce_pm1n_to_pmhn)
EXTEND_VEC_DEF2(vec2d, vec4d, _reduce_pm1no_to_0n)
EXTEND_VEC_DEF2(vec2d, vec4d, _unpacklo)
EXTEND_VEC_DEF2(vec2d, vec4d, _unpackhi)
EXTEND_VEC_DEF3(vec2d, vec4d, _reduce_to_pm1n)
EXTEND_VEC_DEF3(vec2d, vec4d, _reduce_to_pm1no)
EXTEND_VEC_DEF3(vec2d, vec4d, _reduce_to_0n)
EXTEND_VEC_DEF3(vec2d, vec4d, _fmadd)
EXTEND_VEC_DEF3(vec2d, vec4d, _fmsub)
EXTEND_VEC_DEF3(vec2d, vec4d, _fnmadd)
EXTEND_VEC_DEF3(vec2d, vec4d, _fnmsub)
EXTEND_VEC_DEF3(vec2d, vec4d, _blendv)
EXTEND_VEC_DEF4(vec2d, vec4d, _mulmod)
EXTEND_VEC_DEF4(vec2d, vec4d, _nmulmod)


/* vec8d -- NEON/ARM64 *********************************************/

FLINT_FORCE_INLINE vec8d vec8d_set_vec4d2(vec4d a, vec4d b) {
    vec8d z = {a, b}; return z;
}

FLINT_FORCE_INLINE double vec8d_get_index(vec8d a, int i) {
    return i < 4 ? vec4d_get_index(a.e1, i) : vec4d_get_index(a.e2, i - 4);
}

FLINT_FORCE_INLINE vec8d vec8d_set_d8(double a0, double a1, double a2, double a3,
                                      double a4, double a5, double a6, double a7)
{
    vec4d z1 = vec4d_set_d4(a0, a1, a2, a3);
    vec4d z2 = vec4d_set_d4(a4, a5, a6, a7);
    vec8d z = {z1, z2}; return z;
}

FLINT_FORCE_INLINE vec8d vec8d_set_d(double a)
{
    vec4d z1 = vec4d_set_d(a);
    vec8d z = {z1, z1}; return z;
}

FLINT_FORCE_INLINE void vec8d_store(double* z, vec8d a)
{
    vec4d_store(z+0, a.e1);
    vec4d_store(z+4, a.e2);
}

FLINT_FORCE_INLINE void vec8d_store_aligned(double* z, vec8d a)
{
    vec4d_store_aligned(z+0, a.e1);
    vec4d_store_aligned(z+4, a.e2);
}

FLINT_FORCE_INLINE void vec8d_store_unaligned(double* z, vec8d a)
{
    vec4d_store_unaligned(z+0, a.e1);
    vec4d_store_unaligned(z+4, a.e2);
}


FLINT_FORCE_INLINE vec8d vec8d_load(const double* a)
{
    vec4d z1 = vec4d_load(a+0);
    vec4d z2 = vec4d_load(a+4);
    vec8d z = {z1, z2}; return z;
}

FLINT_FORCE_INLINE vec8d vec8d_load_aligned(const double* a)
{
    vec4d z1 = vec4d_load_aligned(a+0);
    vec4d z2 = vec4d_load_aligned(a+4);
    vec8d z = {z1, z2}; return z;
}

FLINT_FORCE_INLINE vec8d vec8d_load_unaligned(const double* a)
{
    vec4d z1 = vec4d_load_unaligned(a+0);
    vec4d z2 = vec4d_load_unaligned(a+4);
    vec8d z = {z1, z2}; return z;
}


EXTEND_VEC_DEF0(vec4d, vec8d, _zero)
EXTEND_VEC_DEF0(vec4d, vec8d, _one)
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
EXTEND_VEC_DEF2(vec4d, vec8d, _unpacklo)
EXTEND_VEC_DEF2(vec4d, vec8d, _unpackhi)
EXTEND_VEC_DEF3(vec4d, vec8d, _reduce_to_pm1n)
EXTEND_VEC_DEF3(vec4d, vec8d, _reduce_to_pm1no)
EXTEND_VEC_DEF3(vec4d, vec8d, _reduce_to_0n)
EXTEND_VEC_DEF3(vec4d, vec8d, _fmadd)
EXTEND_VEC_DEF3(vec4d, vec8d, _fmsub)
EXTEND_VEC_DEF3(vec4d, vec8d, _fnmadd)
EXTEND_VEC_DEF3(vec4d, vec8d, _fnmsub)
EXTEND_VEC_DEF3(vec4d, vec8d, _blendv)
EXTEND_VEC_DEF4(vec4d, vec8d, _mulmod)
EXTEND_VEC_DEF4(vec4d, vec8d, _nmulmod)




FLINT_FORCE_INLINE int vec2d_same(vec2d a, vec2d b) {
    return a[0] == b[0] && a[1] == b[1];
}


#define DEFINE_IT(V) \
FLINT_FORCE_INLINE int V##_same_mod(V a, V b, V n, V ninv) { \
    return V##_same(V##_reduce_to_0n(a, n, ninv), V##_reduce_to_0n(b, n, ninv)); \
}
DEFINE_IT(vec1d)
DEFINE_IT(vec2d)
#undef DEFINE_IT

/* integer stuff *************************************************************/

/* vec1n -- NEON/ARM64 *********************************************/

FLINT_FORCE_INLINE void vec1n_store_unaligned(ulong* z, vec1n a) {
    z[0] = a;
}

FLINT_FORCE_INLINE vec1n vec1d_convert_limited_vec1n(vec1d a) {
    return (slong)a;
}

// (a + b) % n
FLINT_FORCE_INLINE vec1n vec1n_addmod(vec1n a, vec1n b, vec1n n) {
    vec1n nmb = n - b;
    return nmb > a ? a + b : a - nmb;
}

/* vec2n -- NEON/ARM64 *********************************************/

FLINT_FORCE_INLINE vec2d vec2n_convert_limited_vec2d(vec2n a) {
    float64x2_t t = vdupq_n_f64(0x1.0p52);
    return vsubq_f64(vreinterpretq_f64_u64(vorrq_u64(a, vreinterpretq_u64_f64(t))), t);
}

FLINT_FORCE_INLINE void vec2n_store_unaligned(ulong* z, vec2n a) {
   vst1q_u64(z, a);
}

FLINT_FORCE_INLINE vec2n vec2d_convert_limited_vec2n(vec2d a) {
    return vcvtnq_u64_f64(a);
}

FLINT_FORCE_INLINE vec2n vec2n_set_n(ulong a) {
    vec2n x = vdupq_n_u64(a);
    return x;
}

FLINT_FORCE_INLINE vec2n vec2n_load_unaligned(const ulong* a) {
    return vld1q_u64(a);
}

// Right shift 32bits
FLINT_FORCE_INLINE vec2n vec2n_bit_shift_right_32(vec2n a) {
    return vshrq_n_u64(a, 32);
}

// AND operation
FLINT_FORCE_INLINE vec2n vec2n_bit_and(vec2n a, vec2n b) {
    return vandq_u64(a, b);
}

// Addition
FLINT_FORCE_INLINE vec2n vec2n_add(vec2n a, vec2n b) {
    return vaddq_u64(a, b);
}

// Substraction
FLINT_FORCE_INLINE vec2n vec2n_sub(vec2n a, vec2n b) {
    return vsubq_u64(a, b);
}

// (a + b) % n
FLINT_FORCE_INLINE vec2n vec2n_addmod(vec2n a, vec2n b, vec2n n) {
    vec2n nmb = vec2n_sub(n, b);
    vec2n sum = vec2n_sub(a, nmb);

    vec2n mask = vcgtq_u64(nmb, a);

    return vec2n_add(sum, vandq_u64(n, mask));
}

// (a + b) % n for n < 2^63
FLINT_FORCE_INLINE vec2n vec2n_addmod_limited(vec2n a, vec2n b, vec2n n) {
    vec2n s = vec2n_add(a, b);

    vec2n mask = vcgeq_u64(s, n);

    return vec2n_sub(s, vandq_u64(n, mask));
}

/* vec4n -- NEON/ARM64 *********************************************/

FLINT_FORCE_INLINE vec4d vec4n_convert_limited_vec4d(vec4n a) {
    vec2d z1 = vec2n_convert_limited_vec2d(a.e1);
    vec2d z2 = vec2n_convert_limited_vec2d(a.e2);
    vec4d z = {z1, z2}; return z;
}

FLINT_FORCE_INLINE void vec4n_store_unaligned(ulong* z, vec4n a) {
    vec2n_store_unaligned(z+0, a.e1);
    vec2n_store_unaligned(z+2, a.e2);
}

FLINT_FORCE_INLINE vec4n vec4d_convert_limited_vec4n(vec4d a) {
    vec2n z1 = vec2d_convert_limited_vec2n(a.e1);
    vec2n z2 = vec2d_convert_limited_vec2n(a.e2);
    vec4n z = {z1, z2}; return z;
}

FLINT_FORCE_INLINE vec4n vec4n_set_n(ulong a) {
    vec2n x = vec2n_set_n(a);
    vec4n z = {x, x};
    return z;
}

FLINT_FORCE_INLINE vec4n vec4n_load_unaligned(const ulong* a) {
    vec2n z1 = vec2n_load_unaligned(a+0);
    vec2n z2 = vec2n_load_unaligned(a+2);
    vec4n z = {z1, z2}; return z;
}

EXTEND_VEC_DEF1(vec2n, vec4n, _bit_shift_right_32)
EXTEND_VEC_DEF2(vec2n, vec4n, _bit_and)
EXTEND_VEC_DEF2(vec2n, vec4n, _add)
EXTEND_VEC_DEF2(vec2n, vec4n, _sub)
EXTEND_VEC_DEF3(vec2n, vec4n, _addmod)
EXTEND_VEC_DEF3(vec2n, vec4n, _addmod_limited)

/* vec8n -- NEON/ARM64 *********************************************/

FLINT_FORCE_INLINE vec8d vec8n_convert_limited_vec8d(vec8n a) {
    vec4d z1 = vec4n_convert_limited_vec4d(a.e1);
    vec4d z2 = vec4n_convert_limited_vec4d(a.e2);
    vec8d z = {z1, z2}; return z;
}

FLINT_FORCE_INLINE vec8n vec8n_set_n(ulong a) {
    vec4n x = vec4n_set_n(a);
    vec8n z = {x, x};
    return z;
}

FLINT_FORCE_INLINE vec8n vec8n_load_unaligned(const ulong* a)
{
    vec4n z1 = vec4n_load_unaligned(a+0);
    vec4n z2 = vec4n_load_unaligned(a+4);
    vec8n z = {z1, z2}; return z;
}


EXTEND_VEC_DEF1(vec4n, vec8n, _bit_shift_right_32)
EXTEND_VEC_DEF2(vec4n, vec8n, _bit_and)
EXTEND_VEC_DEF2(vec4n, vec8n, _add)
EXTEND_VEC_DEF2(vec4n, vec8n, _sub)
EXTEND_VEC_DEF3(vec4n, vec8n, _addmod)
EXTEND_VEC_DEF3(vec4n, vec8n, _addmod_limited)


#if 0
/*
    vshrq_n_u64(a, n) cannot be used because n must be a compile-time
    constant, and the compiler doesn't see that n is constant
    even if the function is forced inline.
    And vshlq_s64(a, vdupq_n_s64(-(slong) n)) cannot be used to emulate
    vshrq_n_u64(a, n) as it propagates the sign bit.
*/
FLINT_FORCE_INLINE vec2n vec2n_bit_shift_right(vec2n a, ulong n)
{
    return vshrq_n_u64(a, n);
}

FLINT_FORCE_INLINE vec4n vec4n_bit_shift_right(vec4n a, ulong n)
{
    vec4n z = {vec2n_bit_shift_right(a.e1, n), vec2n_bit_shift_right(a.e2, n)};
    return z;
}

FLINT_FORCE_INLINE vec8n vec8n_bit_shift_right(vec8n a, ulong n)
{
    vec8n z = {vec4n_bit_shift_right(a.e1, n), vec4n_bit_shift_right(a.e2, n)};
    return z;
}
#endif


#undef EXTEND_VEC_DEF4
#undef EXTEND_VEC_DEF3
#undef EXTEND_VEC_DEF2
#undef EXTEND_VEC_DEF1
#undef EXTEND_VEC_DEF0


#else

#error machine_vector.h not implemented

#endif


#ifdef __cplusplus
}
#endif

#endif
