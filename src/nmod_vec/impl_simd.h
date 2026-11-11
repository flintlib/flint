/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Private helpers of the SIMD kernels of nmod_vec (reduce.c, ...): the
    instruction sets they use, chosen at compile time from the compiler's
    macros as in nmod_vec.h (so no SIMD code in public headers), and their
    common operations.

    Notation: n of nu bits, 2^(nu - 1) <= n < 2^nu, and
    m = floor((2^(nu + 63) - 1) / n) = 2^63 + floor(ninv / 2),
    the constant of the Barrett step of NMOD_RED; its top half m >> 32 gives
    the quotients of the kernels for nu >= 33, and its shifts the constants of
    those below.

    _nmod_vec_add, _sub and _neg have no kernels: the compilers vectorize
    their loops, and hand-written kernels gained only 1.0-1.2x on long
    vectors and 1.2-1.5x on short ones (Nov 2026).
*/

#ifndef NMOD_VEC_IMPL_SIMD_H
#define NMOD_VEC_IMPL_SIMD_H

#include <stdint.h>
#include "flint.h"
#include "nmod.h"

#if FLINT_BITS == 64 \
        && !defined(FLINT_MACHINE_VECTORS_FORCE_GENERIC) \
        && !defined(FLINT_MACHINE_VECTORS_STRICT_C)
# if defined(__AVX512F__) && defined(__AVX512DQ__) && defined(__AVX512VL__)
#  define NMOD_VEC_SIMD_ZMM 1
# endif
# if defined(__AVX2__)
#  define NMOD_VEC_SIMD_YMM 1
# endif
# if defined(__ARM_NEON) && defined(__aarch64__)
#  define NMOD_VEC_SIMD_NEON 1
# endif
#endif

#ifndef NMOD_VEC_SIMD_ZMM
# define NMOD_VEC_SIMD_ZMM 0
#endif
#ifndef NMOD_VEC_SIMD_YMM
# define NMOD_VEC_SIMD_YMM 0
#endif
#ifndef NMOD_VEC_SIMD_NEON
# define NMOD_VEC_SIMD_NEON 0
#endif

/* AVX-512 extensions used where present */
#if NMOD_VEC_SIMD_ZMM && defined(__AVX512IFMA__)
# define NMOD_VEC_SIMD_IFMA 1
#else
# define NMOD_VEC_SIMD_IFMA 0
#endif
#if NMOD_VEC_SIMD_ZMM && defined(__AVX512VBMI2__)
# define NMOD_VEC_SIMD_VBMI2 1
#else
# define NMOD_VEC_SIMD_VBMI2 0
#endif

#if NMOD_VEC_SIMD_ZMM || NMOD_VEC_SIMD_YMM
# include <immintrin.h>
#endif
#if NMOD_VEC_SIMD_NEON
# include <arm_neon.h>
#endif

/*
    Below this length the operations run their scalar loops: inside FLINT
    functions, where the length varies from call to call, the kernels lost
    on vectors of fewer than 8 entries (0.43-0.94x, Nov 2026), although an
    isolated benchmark at one fixed length shows a gain from 4 entries.
*/
#define NMOD_VEC_SIMD_MIN_LEN 8

/* the constant m above */
#define NMOD_VEC_SIMD_M(mod) ((UWORD(1) << 63) | ((mod).ninv >> 1))

/* for the scalar tails of the AVX2 and NEON kernels: clang would vectorize
   them, which costs more than it saves on at most 3 entries */
#if defined(__clang__)
# define NMOD_VEC_NO_VECTORIZE _Pragma("clang loop vectorize(disable)")
#else
# define NMOD_VEC_NO_VECTORIZE
#endif

/* AVX-512 ****************************************************************/

#if NMOD_VEC_SIMD_ZMM

#define _nv_zmm_set1(x) _mm512_set1_epi64((long long) (x))
#define _nv_zmm_ld(p) _mm512_loadu_si512((const void *) (p))
#define _nv_zmm_st(p, x) _mm512_storeu_si512((void *) (p), (x))
/* masked load (zeros in the other lanes) and store of t < 8 entries */
#define _nv_zmm_mask(t) ((__mmask8) ((1u << (t)) - 1))
#define _nv_zmm_mld(k, p) _mm512_maskz_loadu_epi64((k), (const void *) (p))
#define _nv_zmm_mst(p, k, x) _mm512_mask_storeu_epi64((void *) (p), (k), (x))

/* a >= b ? a - b : a, for any a and b > 0 */
FLINT_FORCE_INLINE __m512i _nv_zmm_csub(__m512i a, __m512i b)
{
    return _mm512_min_epu64(a, _mm512_sub_epi64(a, b));
}

/*
    The loops of the AVX-512 kernels. From NMOD_VEC_SIMD_PEEL_MIN entries,
    when the output is not 64-byte aligned and every vector has the same
    offset modulo 64 bytes, one masked step first does the entries up to the
    next 64-byte boundary, then the pointers are advanced: the main loop
    then loads and stores whole cache lines (vectors 16 bytes past a
    boundary, Ice Lake, Cascade Lake and Zen 4, Nov 2026: 1.01-1.09x on
    average over the kernels, up to 1.4-1.7x). With different offsets
    nothing is peeled: aligning the output alone lost up to 16%. Advancing the
    pointers keeps the main loop as without peeling, and the unroll pragma has
    clang unroll it as it does without peeling (with the loop starting at the
    peeled index instead, clang 19 did not unroll it: 3-12% slower on aligned
    data). The length is tested first, so that short vectors skip all of this
    (4-32 entries: 0.88-0.96x otherwise). The last len mod 8 entries are one
    masked step.  STEP1(x) and STEP2(x, y) are expressions of the vectors of
    entries; the macros advance r, a, b and decrease len, so they take local
    copies.
*/
#define NMOD_VEC_SIMD_PEEL_MIN 128

#define _NV_ZMM_SAME64(p, q) (((((uintptr_t) (p)) ^ ((uintptr_t) (q))) & 63) == 0)
#define _NV_ZMM_PEEL_H(r) ((slong) (((-(uintptr_t) (r)) & 63) >> 3))

#if defined(__clang__)
# define _NV_ZMM_UNROLL _Pragma("clang loop unroll_count(4)")
#else
# define _NV_ZMM_UNROLL
#endif

#define _NV_ZMM_LOOP1(STEP1, r, a, len) \
    do { \
        slong _i; \
        if ((len) >= NMOD_VEC_SIMD_PEEL_MIN) \
        { \
            slong _h = _NV_ZMM_PEEL_H(r); \
            if (_h != 0 && _NV_ZMM_SAME64(r, a)) \
            { \
                __mmask8 _k = _nv_zmm_mask(_h); \
                _nv_zmm_mst((r), _k, STEP1(_nv_zmm_mld(_k, (a)))); \
                (r) += _h; (a) += _h; (len) -= _h; \
            } \
        } \
        _NV_ZMM_UNROLL \
        for (_i = 0; _i + 8 <= (len); _i += 8) \
            _nv_zmm_st((r) + _i, STEP1(_nv_zmm_ld((a) + _i))); \
        if (_i < (len)) \
        { \
            __mmask8 _k = _nv_zmm_mask((len) - _i); \
            _nv_zmm_mst((r) + _i, _k, STEP1(_nv_zmm_mld(_k, (a) + _i))); \
        } \
    } while (0)

#define _NV_ZMM_LOOP2(STEP2, r, a, b, len) \
    do { \
        slong _i; \
        if ((len) >= NMOD_VEC_SIMD_PEEL_MIN) \
        { \
            slong _h = _NV_ZMM_PEEL_H(r); \
            if (_h != 0 && _NV_ZMM_SAME64(r, a) && _NV_ZMM_SAME64(r, b)) \
            { \
                __mmask8 _k = _nv_zmm_mask(_h); \
                _nv_zmm_mst((r), _k, STEP2(_nv_zmm_mld(_k, (a)), _nv_zmm_mld(_k, (b)))); \
                (r) += _h; (a) += _h; (b) += _h; (len) -= _h; \
            } \
        } \
        _NV_ZMM_UNROLL \
        for (_i = 0; _i + 8 <= (len); _i += 8) \
            _nv_zmm_st((r) + _i, STEP2(_nv_zmm_ld((a) + _i), _nv_zmm_ld((b) + _i))); \
        if (_i < (len)) \
        { \
            __mmask8 _k = _nv_zmm_mask((len) - _i); \
            _nv_zmm_mst((r) + _i, _k, STEP2(_nv_zmm_mld(_k, (a) + _i), _nv_zmm_mld(_k, (b) + _i))); \
        } \
    } while (0)

/* the same for the products by a scalar: STEPS(x) is the product of the
   entries x of vec, and with add the entries s of r are added, the sum
   corrected by csub with SN (n) */
#define _NV_ZMM_SADD(t, r, add, SN) \
    ((add) ? _nv_zmm_csub(_mm512_add_epi64((t), _nv_zmm_ld(r)), (SN)) : (t))
#define _NV_ZMM_SADDM(t, r, k, add, SN) \
    ((add) ? _nv_zmm_csub(_mm512_add_epi64((t), _nv_zmm_mld((k), (r))), (SN)) : (t))

#define _NV_ZMM_SLOOP(STEPS, r, a, len, add, SN) \
    do { \
        slong _i; \
        if ((len) >= NMOD_VEC_SIMD_PEEL_MIN) \
        { \
            slong _h = _NV_ZMM_PEEL_H(r); \
            if (_h != 0 && _NV_ZMM_SAME64(r, a)) \
            { \
                __mmask8 _k = _nv_zmm_mask(_h); \
                __m512i _t = STEPS(_nv_zmm_mld(_k, (a))); \
                _nv_zmm_mst((r), _k, _NV_ZMM_SADDM(_t, (r), _k, (add), (SN))); \
                (r) += _h; (a) += _h; (len) -= _h; \
            } \
        } \
        _NV_ZMM_UNROLL \
        for (_i = 0; _i + 8 <= (len); _i += 8) \
        { \
            __m512i _t = STEPS(_nv_zmm_ld((a) + _i)); \
            _nv_zmm_st((r) + _i, _NV_ZMM_SADD(_t, (r) + _i, (add), (SN))); \
        } \
        if (_i < (len)) \
        { \
            __mmask8 _k = _nv_zmm_mask((len) - _i); \
            __m512i _t = STEPS(_nv_zmm_mld(_k, (a) + _i)); \
            _nv_zmm_mst((r) + _i, _k, _NV_ZMM_SADDM(_t, (r) + _i, _k, (add), (SN))); \
        } \
    } while (0)

/* the high word of x y, exact (4 products of 32-bit halves); yh = y with
   its halves swapped */
FLINT_FORCE_INLINE __m512i _nv_zmm_mulhi(__m512i x, __m512i y, __m512i yh)
{
    const __m512i lo = _mm512_set1_epi64(0xFFFFFFFF);
    __m512i xh = _mm512_shuffle_epi32(x, 0xB1);
    __m512i w0 = _mm512_mul_epu32(x, y);
    __m512i w1 = _mm512_mul_epu32(x, yh);
    __m512i w2 = _mm512_mul_epu32(xh, y);
    __m512i w3 = _mm512_mul_epu32(xh, yh);
    __m512i s1 = _mm512_add_epi64(w1, _mm512_srli_epi64(w0, 32));
    __m512i s2 = _mm512_add_epi64(w2, _mm512_and_si512(s1, lo));
    return _mm512_add_epi64(_mm512_add_epi64(w3, _mm512_srli_epi64(s1, 32)),
                            _mm512_srli_epi64(s2, 32));
}

/* the high word of x y minus 0, 1 or 2 (3 products: the product of the
   low halves and the carries of the low word are left out) */
FLINT_FORCE_INLINE __m512i _nv_zmm_mulhi_lazy(__m512i x, __m512i y, __m512i yh)
{
    __m512i xh = _mm512_shuffle_epi32(x, 0xB1);
    __m512i w1 = _mm512_mul_epu32(x, yh);
    __m512i w2 = _mm512_mul_epu32(xh, y);
    __m512i w3 = _mm512_mul_epu32(xh, yh);
    return _mm512_add_epi64(w3, _mm512_add_epi64(_mm512_srli_epi64(w1, 32),
                                                 _mm512_srli_epi64(w2, 32)));
}

/* the full product: returns the high word, the low one in *lo */
FLINT_FORCE_INLINE __m512i _nv_zmm_mul_full(__m512i x, __m512i y, __m512i * lo)
{
    const __m512i lm = _mm512_set1_epi64(0xFFFFFFFF);
    __m512i xh = _mm512_shuffle_epi32(x, 0xB1);
    __m512i yh = _mm512_shuffle_epi32(y, 0xB1);
    __m512i w0 = _mm512_mul_epu32(x, y);
    __m512i w1 = _mm512_mul_epu32(x, yh);
    __m512i w2 = _mm512_mul_epu32(xh, y);
    __m512i w3 = _mm512_mul_epu32(xh, yh);
    __m512i s1 = _mm512_add_epi64(w1, _mm512_srli_epi64(w0, 32));
    __m512i s2 = _mm512_add_epi64(w2, _mm512_and_si512(s1, lm));
    *lo = _mm512_or_si512(_mm512_slli_epi64(s2, 32), _mm512_and_si512(w0, lm));
    return _mm512_add_epi64(_mm512_add_epi64(w3, _mm512_srli_epi64(s1, 32)),
                            _mm512_srli_epi64(s2, 32));
}

#endif /* NMOD_VEC_SIMD_ZMM */

/* AVX2 *******************************************************************/

#if NMOD_VEC_SIMD_YMM

#define _nv_ymm_set1(x) _mm256_set1_epi64x((long long) (x))
#define _nv_ymm_ld(p) _mm256_loadu_si256((const __m256i *) (p))
#define _nv_ymm_st(p, x) _mm256_storeu_si256((__m256i *) (p), (x))

/*
    Two rules for clang, without which it makes some AVX2 kernels up to 2x
    slower than gcc (clang 18 and 19, Nov 2026):
    1. a constant that is only a 32-bit factor of _mm256_mul_epu32 is set
       from a 32-bit value and hidden from the optimizer (_nv_ymm_set1_32):
       otherwise clang may turn x - q n into x + q (-n) with a 64-bit factor
       and emulate a 64-bit product (three products instead of one);
    2. no test on nu or on a flag inside a loop: each kernel is a
       FLINT_FORCE_INLINE body with constant flags, called once per case.
*/
FLINT_FORCE_INLINE __m256i _nv_ymm_set1_32(ulong x)
{
    __m256i v = _mm256_set1_epi64x((long long) (uint32_t) x);
#if defined(__GNUC__)
    __asm__("" : "+x" (v));
#endif
    return v;
}

/* a >= b ? a - b : a, for a, b < 2^63 */
FLINT_FORCE_INLINE __m256i _nv_ymm_csub_s(__m256i a, __m256i b)
{
    __m256i lt = _mm256_cmpgt_epi64(b, a);
    return _mm256_sub_epi64(a, _mm256_andnot_si256(lt, b));
}

/* the same for any a and b, with bx = b ^ 2^63 */
FLINT_FORCE_INLINE __m256i _nv_ymm_csub_u(__m256i a, __m256i b, __m256i bx)
{
    __m256i ax = _mm256_xor_si256(a, _mm256_set1_epi64x(INT64_MIN));
    __m256i lt = _mm256_cmpgt_epi64(bx, ax);
    return _mm256_sub_epi64(a, _mm256_andnot_si256(lt, b));
}

/* the same on 32-bit lanes, for a, b < 2^32 */
FLINT_FORCE_INLINE __m256i _nv_ymm_csub_32(__m256i a, __m256i b)
{
    return _mm256_min_epu32(a, _mm256_sub_epi32(a, b));
}

#endif /* NMOD_VEC_SIMD_YMM */

/* NEON *******************************************************************/

#if NMOD_VEC_SIMD_NEON

#define _nv_neon_ld(p) vld1q_u64((const uint64_t *) (p))
#define _nv_neon_st(p, x) vst1q_u64((uint64_t *) (p), (x))

/* a >= b ? a - b : a, for any a and b */
FLINT_FORCE_INLINE uint64x2_t _nv_neon_csub(uint64x2_t a, uint64x2_t b)
{
    return vsubq_u64(a, vandq_u64(b, vcgeq_u64(a, b)));
}

/*
    On the M4, interleaving 4 entries in NEON with 4 entries in scalar code
    (the scalar multiplier runs beside the NEON one) nearly doubles the
    speed of the reduction and of the products below 2^32: from 20 entries,
    1.4-2.1x and 2.7-3.7x over the scalar code, against 1.2-1.6x and
    1.6-2.4x for the NEON loops alone (Nov 2026); at 4 entries the hybrid
    loops ran at 0.71-0.93x. They replace the NEON ones from this length.
    The NEON kernels follow the second clang rule of the AVX2 ones (clang
    is the compiler of macOS).
*/
#define NMOD_VEC_SIMD_HYBRID_MIN 16

#endif /* NMOD_VEC_SIMD_NEON */

#endif
