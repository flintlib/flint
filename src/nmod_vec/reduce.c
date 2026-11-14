/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "impl_simd.h"

/*
    The SIMD kernels, by the size nu of n, for any input words, with the
    constant m of NMOD_RED (impl_simd.h) and products of 32-bit halves.
    With 32-bit words, a is a two-word integer and, for nu <= 32, n a
    one-word modulus: the kernels are reductions of
    doc: algorithms/arithmetic/barrett at that word size. Bounds and
    proofs: doc: nmod_vec/simd#reduce.
      fullword (nu = 64): a < 2n, one correction, as NMOD_RED.
      top32 (33 <= nu <= 63): q from the top halves of a and m, as NMOD_RED
        with a truncated estimate; one correction (two at nu = 33).
      fold29 (6 <= nu <= 32): the bits from nu + 29 on folded in with
        2^(nu + 29) mod n, then the step of NMOD_RED2_BOUNDED with 32-bit
        words (doc: algorithms/arithmetic/barrett#thm-bounded); one
        correction.
      fold32 (2 <= nu <= 5): the high half folded in with 2^32 mod n, as in
        n_lll_mod_barrett (doc: algorithms/arithmetic/barrett#two-words),
        then a coarser Barrett step; two corrections (2n, then n).
        NMOD2_RED2 would instead reduce the high half, then take a
        Moller-Granlund step.
    _nmod_vec_reduce_lazy (n < 2^(FLINT_BITS - 2)) skips the last correction of each
    kernel, so its outputs are in [0, 2n): 1.11-1.64x on long vectors
    (Nov 2026).
    The Barrett step of NMOD_RED itself, with the 64 x 64 high product
    emulated in SIMD, ran 1.3-1.6x slower than these (Nov 2026).
*/

/* the scalar code, also for the tails of the AVX2 and NEON kernels */
static void
_nmod_vec_reduce_scalar(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
    slong i;

    /* In place, a loop that reads res itself (see _nmod_vec_add): 1.1-2.1x
       faster in place with clang, 1.07-1.17x with gcc (Nov 2026) */
    if (res == vec)
        for (i = 0 ; i < len; i++)
            NMOD_RED(res[i], res[i], mod);
    else
        for (i = 0 ; i < len; i++)
            NMOD_RED(res[i], vec[i], mod);
}

#if !NMOD_VEC_SIMD_ZMM && (NMOD_VEC_SIMD_YMM || NMOD_VEC_SIMD_NEON)
FLINT_FORCE_INLINE void
_reduce_tail(nn_ptr res, nn_srcptr vec, slong i, slong len, nmod_t mod)
{
    NMOD_VEC_NO_VECTORIZE
    for ( ; i < len; i++)
        NMOD_RED(res[i], vec[i], mod);
}
#endif

/* AVX-512 ****************************************************************/

#if NMOD_VEC_SIMD_ZMM

/*
    q n is one vpmullq. Measured and not kept (Nov 2026): q n as two
    vpmuludq (q < 2^32), 30-45% slower with clang and 27% slower on Zen 4;
    for 33 <= nu <= 51 with IFMA, a - q n as (a + q (2^52 - n)) mod 2^52,
    one madd52lo and a mask instead of vpmullq (exact while 2n < 2^52):
    about 1.27x on Ice Lake, where vpmullq is slow (3 uops), but 2-8%
    slower on Zen 4.
*/
FLINT_FORCE_INLINE __m512i
_reduce_top32_step_zmm(__m512i x, __m512i vm, __m512i vs, __m512i vn, int two,
                       int lazy)
{
    __m512i q = _mm512_srlv_epi64(_mm512_mul_epu32(_mm512_shuffle_epi32(x, 0xB1), vm), vs);
    x = _mm512_sub_epi64(x, _mm512_mullo_epi64(q, vn));
    if (two)
        x = _nv_zmm_csub(x, vn);
    return lazy ? x : _nv_zmm_csub(x, vn);
}

FLINT_FORCE_INLINE void
_reduce_top32_zmm(nn_ptr r, nn_srcptr a, slong len, nmod_t mod, int two, int lazy)
{
    const __m512i vn = _nv_zmm_set1(mod.n);
    const __m512i vm = _nv_zmm_set1(NMOD_VEC_SIMD_M(mod) >> 32);
    const __m512i vs = _nv_zmm_set1(FLINT_BITS - 1 - mod.norm);
#define STEP(x) _reduce_top32_step_zmm((x), vm, vs, vn, two, lazy)
    _NV_ZMM_LOOP1(STEP, r, a, len);
#undef STEP
}

FLINT_FORCE_INLINE void
_reduce_fullword_zmm(nn_ptr r, nn_srcptr a, slong len, nmod_t mod)
{
    const __m512i vn = _nv_zmm_set1(mod.n);
#define STEP(x) _nv_zmm_csub((x), vn)
    _NV_ZMM_LOOP1(STEP, r, a, len);
#undef STEP
}

FLINT_FORCE_INLINE __m512i
_reduce_fold29_step_zmm(__m512i x, __m512i vsh, __m512i vc, __m512i vmask,
                        __m512i vs, __m512i vm, __m512i vn, int lazy)
{
    __m512i t = _mm512_add_epi64(_mm512_mul_epu32(_mm512_srlv_epi64(x, vsh), vc),
                                 _mm512_and_si512(x, vmask));
    __m512i q = _mm512_srli_epi64(_mm512_mul_epu32(_mm512_srlv_epi64(t, vs), vm), 33);
    t = _mm512_sub_epi64(t, _mm512_mul_epu32(q, vn));
    return lazy ? t : _nv_zmm_csub(t, vn);
}

FLINT_FORCE_INLINE void
_reduce_fold29_zmm(nn_ptr r, nn_srcptr a, slong len, nmod_t mod, int lazy)
{
    const int nu = FLINT_BITS - mod.norm;
    ulong c;
    NMOD_RED(c, UWORD(1) << (nu + 29), mod);
    {
        const __m512i vn = _nv_zmm_set1(mod.n);
        const __m512i vc = _nv_zmm_set1(c);
        const __m512i vmask = _nv_zmm_set1((UWORD(1) << (nu + 29)) - 1);
        const __m512i vm = _nv_zmm_set1(NMOD_VEC_SIMD_M(mod) >> 32);
        const __m512i vsh = _nv_zmm_set1(nu + 29);
        const __m512i vs = _nv_zmm_set1(nu - 2);
#define STEP(x) _reduce_fold29_step_zmm((x), vsh, vc, vmask, vs, vm, vn, lazy)
        _NV_ZMM_LOOP1(STEP, r, a, len);
#undef STEP
    }
}

FLINT_FORCE_INLINE __m512i
_reduce_fold32_step_zmm(__m512i x, __m512i vc, __m512i vmp, __m512i vs, __m512i vn,
                        __m512i v2n, int lazy)
{
    const __m512i lo = _mm512_set1_epi64(0xFFFFFFFF);
    __m512i t = _mm512_add_epi64(_mm512_mul_epu32(_mm512_shuffle_epi32(x, 0xB1), vc),
                                 _mm512_and_si512(x, lo));
    __m512i s = _mm512_srlv_epi64(t, vs);
    __m512i q = _mm512_add_epi64(s, _mm512_srli_epi64(_mm512_mul_epu32(s, vmp), 32));
    t = _mm512_sub_epi64(t, _mm512_mul_epu32(q, vn));
    t = _nv_zmm_csub(t, v2n);
    return lazy ? t : _nv_zmm_csub(t, vn);
}

FLINT_FORCE_INLINE void
_reduce_fold32_zmm(nn_ptr r, nn_srcptr a, slong len, nmod_t mod, int lazy)
{
    ulong c;
    NMOD_RED(c, UWORD(1) << 32, mod);
    {
        const __m512i vn = _nv_zmm_set1(mod.n);
        const __m512i v2n = _nv_zmm_set1(2 * mod.n);
        const __m512i vc = _nv_zmm_set1(c);
        const __m512i vmp = _nv_zmm_set1((NMOD_VEC_SIMD_M(mod) >> 31) - (UWORD(1) << 32));
        const __m512i vs = _nv_zmm_set1(FLINT_BITS - mod.norm);
#define STEP(x) _reduce_fold32_step_zmm((x), vc, vmp, vs, vn, v2n, lazy)
        _NV_ZMM_LOOP1(STEP, r, a, len);
#undef STEP
    }
}

/* returns 0 when no kernel applies; lazy (nu <= 62): without the last correction */
FLINT_FORCE_INLINE int
_nmod_vec_reduce_simd(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod, int lazy)
{
    const int nu = FLINT_BITS - mod.norm;

    if (nu == 64)
        _reduce_fullword_zmm(res, vec, len, mod);
    else if (nu >= 34)
        _reduce_top32_zmm(res, vec, len, mod, 0, lazy);
    else if (nu == 33)
        _reduce_top32_zmm(res, vec, len, mod, 1, lazy);
    else if (nu >= 6)
        _reduce_fold29_zmm(res, vec, len, mod, lazy);
    else if (nu >= 2)
        _reduce_fold32_zmm(res, vec, len, mod, lazy);
    else
        return 0;
    return 1;
}

/* AVX2 *******************************************************************/

#elif NMOD_VEC_SIMD_YMM

/*
    nu >= 33: q n as two vpmuludq (q < 2^32); mode 0: 34 <= nu <= 62,
    1: nu = 33 (two corrections), 2: nu = 63 (the remainder may exceed
    2^63: unsigned comparison); lazy (modes 0 and 1): without the last
    correction
*/
FLINT_FORCE_INLINE slong
_reduce_top32_ymm(nn_ptr r, nn_srcptr a, slong len, nmod_t mod, int mode, int lazy)
{
    const __m256i vn = _nv_ymm_set1(mod.n);
    const __m256i vnl = _nv_ymm_set1_32(mod.n);
    const __m256i vnh = _nv_ymm_set1_32(mod.n >> 32);
    const __m256i vnx = _nv_ymm_set1(mod.n ^ (UWORD(1) << 63));
    const __m256i vm = _nv_ymm_set1_32(NMOD_VEC_SIMD_M(mod) >> 32);
    const __m256i vs = _nv_ymm_set1(FLINT_BITS - 1 - mod.norm);
    slong i;

    for (i = 0; i + 4 <= len; i += 4)
    {
        __m256i x = _nv_ymm_ld(a + i);
        __m256i q = _mm256_srlv_epi64(_mm256_mul_epu32(_mm256_shuffle_epi32(x, 0xB1), vm), vs);
        x = _mm256_sub_epi64(x, _mm256_add_epi64(_mm256_mul_epu32(q, vnl),
                                _mm256_slli_epi64(_mm256_mul_epu32(q, vnh), 32)));
        if (mode == 2)
            x = _nv_ymm_csub_u(x, vn, vnx);
        else
        {
            if (mode == 1)
                x = _nv_ymm_csub_s(x, vn);
            if (!lazy)
                x = _nv_ymm_csub_s(x, vn);
        }
        _nv_ymm_st(r + i, x);
    }
    return i;
}

/* nu = 64 */
FLINT_FORCE_INLINE slong
_reduce_fullword_ymm(nn_ptr r, nn_srcptr a, slong len, nmod_t mod)
{
    const __m256i vn = _nv_ymm_set1(mod.n);
    const __m256i vnx = _nv_ymm_set1(mod.n ^ (UWORD(1) << 63));
    slong i;

    for (i = 0; i + 4 <= len; i += 4)
        _nv_ymm_st(r + i, _nv_ymm_csub_u(_nv_ymm_ld(a + i), vn, vnx));
    return i;
}

/* 6 <= nu <= 32; c32: nu <= 31, the correction on 32-bit lanes */
FLINT_FORCE_INLINE slong
_reduce_fold29_ymm(nn_ptr r, nn_srcptr a, slong len, nmod_t mod, int c32, int lazy)
{
    const int nu = FLINT_BITS - mod.norm;
    ulong c;
    slong i;
    NMOD_RED(c, UWORD(1) << (nu + 29), mod);
    {
        const __m256i vn = _nv_ymm_set1_32(mod.n);
        const __m256i vc = _nv_ymm_set1_32(c);
        const __m256i vmask = _nv_ymm_set1((UWORD(1) << (nu + 29)) - 1);
        const __m256i vm = _nv_ymm_set1_32(NMOD_VEC_SIMD_M(mod) >> 32);
        const __m256i vsh = _nv_ymm_set1(nu + 29);
        const __m256i vs = _nv_ymm_set1(nu - 2);

        for (i = 0; i + 4 <= len; i += 4)
        {
            __m256i x = _nv_ymm_ld(a + i);
            __m256i t = _mm256_add_epi64(_mm256_mul_epu32(_mm256_srlv_epi64(x, vsh), vc),
                                         _mm256_and_si256(x, vmask));
            __m256i q = _mm256_srli_epi64(_mm256_mul_epu32(_mm256_srlv_epi64(t, vs), vm), 33);
            t = _mm256_sub_epi64(t, _mm256_mul_epu32(q, vn));
            if (!lazy)
                t = c32 ? _nv_ymm_csub_32(t, vn) : _nv_ymm_csub_s(t, vn);
            _nv_ymm_st(r + i, t);
        }
    }
    return i;
}

/* 2 <= nu <= 5 */
FLINT_FORCE_INLINE slong
_reduce_fold32_ymm(nn_ptr r, nn_srcptr a, slong len, nmod_t mod, int lazy)
{
    ulong c;
    slong i;
    NMOD_RED(c, UWORD(1) << 32, mod);
    {
        const __m256i vn = _nv_ymm_set1_32(mod.n);
        const __m256i v2n = _nv_ymm_set1(2 * mod.n);
        const __m256i vc = _nv_ymm_set1_32(c);
        const __m256i vmp = _nv_ymm_set1_32((NMOD_VEC_SIMD_M(mod) >> 31) - (UWORD(1) << 32));
        const __m256i lo = _nv_ymm_set1(0xFFFFFFFF);
        const __m256i vs = _nv_ymm_set1(FLINT_BITS - mod.norm);

        for (i = 0; i + 4 <= len; i += 4)
        {
            __m256i x = _nv_ymm_ld(a + i);
            __m256i t = _mm256_add_epi64(_mm256_mul_epu32(_mm256_shuffle_epi32(x, 0xB1), vc),
                                         _mm256_and_si256(x, lo));
            __m256i s = _mm256_srlv_epi64(t, vs);
            __m256i q = _mm256_add_epi64(s, _mm256_srli_epi64(_mm256_mul_epu32(s, vmp), 32));
            t = _mm256_sub_epi64(t, _mm256_mul_epu32(q, vn));
            t = _nv_ymm_csub_s(t, v2n);
            if (!lazy)
                t = _nv_ymm_csub_s(t, vn);
            _nv_ymm_st(r + i, t);
        }
    }
    return i;
}

FLINT_FORCE_INLINE int
_nmod_vec_reduce_simd(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod, int lazy)
{
    const int nu = FLINT_BITS - mod.norm;
    slong i;

    if (nu < 2)
        return 0;

    i = _nv_ymm_peel(res, vec, vec, len);
    if (i != 0)
    {
        _reduce_tail(res, vec, 0, i, mod);
        res += i;
        vec += i;
        len -= i;
    }

    if (nu == 64)
        i = _reduce_fullword_ymm(res, vec, len, mod);
    else if (nu == 63)
        i = _reduce_top32_ymm(res, vec, len, mod, 2, 0);
    else if (nu >= 34)
        i = _reduce_top32_ymm(res, vec, len, mod, 0, lazy);
    else if (nu == 33)
        i = _reduce_top32_ymm(res, vec, len, mod, 1, lazy);
    else if (nu == 32 || (lazy && nu >= 6))
        i = _reduce_fold29_ymm(res, vec, len, mod, 0, lazy);
    else if (nu >= 6)
        i = _reduce_fold29_ymm(res, vec, len, mod, 1, 0);
    else
        i = _reduce_fold32_ymm(res, vec, len, mod, lazy);
    _reduce_tail(res, vec, i, len, mod);
    return 1;
}

/* NEON *******************************************************************/

#elif NMOD_VEC_SIMD_NEON

/* 33 <= nu <= 63, four entries per step; two: nu = 33; lazy (nu <= 62):
   without the last correction, here and below */
FLINT_FORCE_INLINE slong
_reduce_top32_neon(nn_ptr r, nn_srcptr a, slong len, nmod_t mod, int two, int lazy)
{
    const uint64x2_t vn = vdupq_n_u64(mod.n);
    const uint32x2_t vnl = vdup_n_u32((uint32_t) mod.n);
    const uint32x2_t vnh = vdup_n_u32((uint32_t) (mod.n >> 32));
    const uint32x2_t vm = vdup_n_u32((uint32_t) (NMOD_VEC_SIMD_M(mod) >> 32));
    const int64x2_t vs = vdupq_n_s64(-(int64_t) (FLINT_BITS - 1 - mod.norm));
    slong i;

    for (i = 0; i + 4 <= len; i += 4)
    {
        uint64x2_t x0 = _nv_neon_ld(a + i), x1 = _nv_neon_ld(a + i + 2);
        uint32x2_t q0 = vmovn_u64(vshlq_u64(vmull_u32(vshrn_n_u64(x0, 32), vm), vs));
        uint32x2_t q1 = vmovn_u64(vshlq_u64(vmull_u32(vshrn_n_u64(x1, 32), vm), vs));
        x0 = vsubq_u64(vmlsl_u32(x0, q0, vnl), vshlq_n_u64(vmull_u32(q0, vnh), 32));
        x1 = vsubq_u64(vmlsl_u32(x1, q1, vnl), vshlq_n_u64(vmull_u32(q1, vnh), 32));
        if (two)
        {
            x0 = _nv_neon_csub(x0, vn);
            x1 = _nv_neon_csub(x1, vn);
        }
        if (!lazy)
        {
            x0 = _nv_neon_csub(x0, vn);
            x1 = _nv_neon_csub(x1, vn);
        }
        _nv_neon_st(r + i, x0);
        _nv_neon_st(r + i + 2, x1);
    }
    return i;
}

/* the same, by blocks of 8: 4 entries in NEON and 4 with NMOD_RED */
FLINT_FORCE_INLINE slong
_reduce_top32_hybrid_neon(nn_ptr r, nn_srcptr a, slong len, nmod_t mod, int two,
                          int lazy)
{
    const uint64x2_t vn = vdupq_n_u64(mod.n);
    const uint32x2_t vnl = vdup_n_u32((uint32_t) mod.n);
    const uint32x2_t vnh = vdup_n_u32((uint32_t) (mod.n >> 32));
    const uint32x2_t vm = vdup_n_u32((uint32_t) (NMOD_VEC_SIMD_M(mod) >> 32));
    const int64x2_t vs = vdupq_n_s64(-(int64_t) (FLINT_BITS - 1 - mod.norm));
    slong i, k;

    for (i = 0; i + 8 <= len; i += 8)
    {
        uint64x2_t x0 = _nv_neon_ld(a + i), x1 = _nv_neon_ld(a + i + 2);
        uint32x2_t q0 = vmovn_u64(vshlq_u64(vmull_u32(vshrn_n_u64(x0, 32), vm), vs));
        uint32x2_t q1 = vmovn_u64(vshlq_u64(vmull_u32(vshrn_n_u64(x1, 32), vm), vs));
        x0 = vsubq_u64(vmlsl_u32(x0, q0, vnl), vshlq_n_u64(vmull_u32(q0, vnh), 32));
        x1 = vsubq_u64(vmlsl_u32(x1, q1, vnl), vshlq_n_u64(vmull_u32(q1, vnh), 32));
        if (two)
        {
            x0 = _nv_neon_csub(x0, vn);
            x1 = _nv_neon_csub(x1, vn);
        }
        if (!lazy)
        {
            x0 = _nv_neon_csub(x0, vn);
            x1 = _nv_neon_csub(x1, vn);
        }
        _nv_neon_st(r + i, x0);
        _nv_neon_st(r + i + 2, x1);
        for (k = 4; k < 8; k++)
            NMOD_RED(r[i + k], a[i + k], mod);
    }
    return i;
}

/* nu = 64 */
FLINT_FORCE_INLINE slong
_reduce_fullword_neon(nn_ptr r, nn_srcptr a, slong len, nmod_t mod)
{
    const uint64x2_t vn = vdupq_n_u64(mod.n);
    slong i;

    for (i = 0; i + 2 <= len; i += 2)
        _nv_neon_st(r + i, _nv_neon_csub(_nv_neon_ld(a + i), vn));
    return i;
}

/* 6 <= nu <= 32, four entries per step (the quotient below is
   ((t >> (nu - 2)) (m >> 32)) >> 32, halved) */
#define _REDUCE_FOLD29_STEP_NEON(t, x) \
    do { \
        uint32x2_t _q; \
        (t) = vmlal_u32(vandq_u64((x), vmask), vmovn_u64(vshlq_u64((x), vsh)), vc); \
        _q = vshrn_n_u64(vmull_u32(vmovn_u64(vshlq_u64((t), vs)), vm), 32); \
        (t) = vmlsl_u32((t), vshr_n_u32(_q, 1), vnl); \
        if (!lazy) \
            (t) = _nv_neon_csub((t), vn); \
    } while (0)

#define _REDUCE_FOLD29_CONSTS_NEON \
    const int nu = FLINT_BITS - mod.norm; \
    ulong c; \
    uint64x2_t vn, vmask; \
    uint32x2_t vnl, vc, vm; \
    int64x2_t vsh, vs; \
    NMOD_RED(c, UWORD(1) << (nu + 29), mod); \
    vn = vdupq_n_u64(mod.n); \
    vnl = vdup_n_u32((uint32_t) mod.n); \
    vc = vdup_n_u32((uint32_t) c); \
    vmask = vdupq_n_u64((UWORD(1) << (nu + 29)) - 1); \
    vm = vdup_n_u32((uint32_t) (NMOD_VEC_SIMD_M(mod) >> 32)); \
    vsh = vdupq_n_s64(-(int64_t) (nu + 29)); \
    vs = vdupq_n_s64(-(int64_t) (nu - 2))

FLINT_FORCE_INLINE slong
_reduce_fold29_neon(nn_ptr r, nn_srcptr a, slong len, nmod_t mod, int lazy)
{
    slong i;
    _REDUCE_FOLD29_CONSTS_NEON;

    for (i = 0; i + 4 <= len; i += 4)
    {
        uint64x2_t x0 = _nv_neon_ld(a + i), x1 = _nv_neon_ld(a + i + 2), t0, t1;
        _REDUCE_FOLD29_STEP_NEON(t0, x0);
        _REDUCE_FOLD29_STEP_NEON(t1, x1);
        _nv_neon_st(r + i, t0);
        _nv_neon_st(r + i + 2, t1);
    }
    return i;
}

/* the same, by blocks of 8: 4 entries in NEON and 4 with NMOD_RED */
FLINT_FORCE_INLINE slong
_reduce_fold29_hybrid_neon(nn_ptr r, nn_srcptr a, slong len, nmod_t mod, int lazy)
{
    slong i, k;
    _REDUCE_FOLD29_CONSTS_NEON;

    for (i = 0; i + 8 <= len; i += 8)
    {
        uint64x2_t x0 = _nv_neon_ld(a + i), x1 = _nv_neon_ld(a + i + 2), t0, t1;
        _REDUCE_FOLD29_STEP_NEON(t0, x0);
        _REDUCE_FOLD29_STEP_NEON(t1, x1);
        _nv_neon_st(r + i, t0);
        _nv_neon_st(r + i + 2, t1);
        for (k = 4; k < 8; k++)
            NMOD_RED(r[i + k], a[i + k], mod);
    }
    return i;
}

/*
    No kernel for nu <= 5: the one of AVX2 ran at 0.92-0.93x on the M4
    (Nov 2026). The hybrid loops go first from NMOD_VEC_SIMD_HYBRID_MIN
    entries, then the NEON loops do the rest (the hybrid alone with the
    scalar tail ran at 0.71x at 4 entries).
*/
FLINT_FORCE_INLINE int
_nmod_vec_reduce_simd(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod, int lazy)
{
    const int nu = FLINT_BITS - mod.norm;
    slong i = 0;

    if (nu < NMOD_VEC_SIMD_REDUCE_MIN_BITS)
        return 0;
    else if (nu == 64)
        i = _reduce_fullword_neon(res, vec, len, mod);
    else if (nu >= 34)
    {
        if (len >= NMOD_VEC_SIMD_HYBRID_MIN)
            i = _reduce_top32_hybrid_neon(res, vec, len, mod, 0, lazy);
        i += _reduce_top32_neon(res + i, vec + i, len - i, mod, 0, lazy);
    }
    else if (nu == 33)
    {
        if (len >= NMOD_VEC_SIMD_HYBRID_MIN)
            i = _reduce_top32_hybrid_neon(res, vec, len, mod, 1, lazy);
        i += _reduce_top32_neon(res + i, vec + i, len - i, mod, 1, lazy);
    }
    else
    {
        if (len >= NMOD_VEC_SIMD_HYBRID_MIN)
            i = _reduce_fold29_hybrid_neon(res, vec, len, mod, lazy);
        i += _reduce_fold29_neon(res + i, vec + i, len - i, mod, lazy);
    }
    _reduce_tail(res, vec, i, len, mod);
    return 1;
}

#endif

#if NMOD_VEC_SIMD_ZMM || NMOD_VEC_SIMD_YMM || NMOD_VEC_SIMD_NEON

/* out of line (impl_simd.h); n has at least NMOD_VEC_SIMD_REDUCE_MIN_BITS bits */
FLINT_STATIC_NOINLINE void
_reduce_simd(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
    _nmod_vec_reduce_simd(res, vec, len, mod, 0);
}

FLINT_STATIC_NOINLINE void
_reduce_lazy_simd(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
    _nmod_vec_reduce_simd(res, vec, len, mod, 1);
}

# define _REDUCE_SIMD_APPLIES(len, mod) ((len) >= NMOD_VEC_SIMD_MIN_LEN \
        && NMOD_BITS(mod) >= NMOD_VEC_SIMD_REDUCE_MIN_BITS)

#endif

void _nmod_vec_reduce(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
#if NMOD_VEC_SIMD_ZMM || NMOD_VEC_SIMD_YMM || NMOD_VEC_SIMD_NEON
    if (_REDUCE_SIMD_APPLIES(len, mod))
        _reduce_simd(res, vec, len, mod);
    else
#endif
        _nmod_vec_reduce_scalar(res, vec, len, mod);
}

/* the scalar code reduces fully, which is also in [0, 2n) */
void _nmod_vec_reduce_lazy(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
    FLINT_ASSERT(NMOD_BITS(mod) <= FLINT_BITS - 2);
#if NMOD_VEC_SIMD_ZMM || NMOD_VEC_SIMD_YMM || NMOD_VEC_SIMD_NEON
    if (_REDUCE_SIMD_APPLIES(len, mod))
        _reduce_lazy_simd(res, vec, len, mod);
    else
#endif
        _nmod_vec_reduce_scalar(res, vec, len, mod);
}
