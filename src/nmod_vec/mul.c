/*
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
    The SIMD kernels, by the size nu of n, for reduced inputs, with the
    constant m of NMOD_RED (impl_simd.h); bounds: doc: nmod_vec/simd#mul.
      halflimb (nu <= 32): p = a b < 2^64, then a 32-bit Barrett step with
        m >> 32 on p >> s; for nu <= 30, the step of NMOD_RED2_BOUNDED with
        32-bit words (doc: algorithms/arithmetic/barrett#thm-bounded) and
        one correction; for nu = 31, 32, two corrections (2n, then n).
      ifma50 (33 <= nu <= 50, AVX-512 with IFMA): the product in 52-bit
        halves and the same step with 52-bit words (the kernel of
        nmod_mat_mul_u52, #2597).
      fp50 (33 <= nu <= 50, otherwise): in doubles, as vec4d_mulmod.
      red2bounded (51 <= nu <= 62, AVX-512): the full product, then the
        reduction of NMOD_RED2_BOUNDED, with a high product that may be 2
        too small: corrections 2n, then n.
    Above, and above 50 bits without AVX-512, the scalar code: an AVX2
    form of red2bounded, with its 64-bit products emulated, ran below the
    scalar code on Broadwell, and on Ice Lake below the AVX-512 one
    (Nov 2026).
    _nmod_vec_mul_lazy (n < 2^62) skips the last correction of each kernel
    (fp50: z + n instead of the correction of z in (-n, n)), so its outputs
    are in [0, 2n): 0.99-1.29x on long vectors (Nov 2026).
*/

/* the loops of #2950 for each size of n */
static void
_nmod_vec_mul_scalar(nn_ptr res, nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
{
    slong i;

    if (NMOD_IS_HALFLIMB(mod))
    {
        for (i = 0; i < len; i++)
            res[i] = _nmod_mul_halflimb(vec1[i], vec2[i], mod);
    }
    else if (NMOD_CAN_USE_MUL_BOUNDED(mod))
    {
        /* as _scalar_mul_bounded in scalar.c: clang vectorizes this loop
           into slower code (x86, AVX2) */
        NMOD_VEC_NO_VECTORIZE
        for (i = 0; i < len; i++)
            res[i] = _nmod_mul_bounded(vec1[i], vec2[i], mod);
    }
    else
    {
        for (i = 0; i < len; i++)
            res[i] = nmod_mul(vec1[i], vec2[i], mod);
    }
}

#if !NMOD_VEC_SIMD_ZMM && (NMOD_VEC_SIMD_YMM || NMOD_VEC_SIMD_NEON)
FLINT_FORCE_INLINE void
_mul_tail(nn_ptr res, nn_srcptr vec1, nn_srcptr vec2, slong i, slong len, nmod_t mod)
{
    NMOD_VEC_NO_VECTORIZE
    for ( ; i < len; i++)
        res[i] = nmod_mul(vec1[i], vec2[i], mod);
}
#endif

/*
    The shifts of the 32-bit Barrett step for nu <= 32: q = ((p >> s)
    (m >> 32)) >> t with s + t = nu + 31; for nu <= 30 the remainder is
    below 2n, for nu = 31, 32 below 4n (two corrections).
*/
#define _HALFLIMB_SHIFTS(s, t, two, nu) \
    do { \
        if ((nu) <= 30) { (s) = (nu) - 2; (two) = 0; } \
        else if ((nu) == 31) { (s) = 30; (two) = 1; } \
        else { (s) = 32; (two) = 1; } \
        (t) = (nu) + 31 - (s); \
    } while (0)

/* AVX-512 ****************************************************************/

#if NMOD_VEC_SIMD_ZMM

FLINT_FORCE_INLINE __m512i
_mul_halflimb_step_zmm(__m512i x, __m512i y, __m512i vs, __m512i vm, __m512i vt,
                       __m512i vn, __m512i v2n, int two, int lazy)
{
    __m512i p = _mm512_mul_epu32(x, y);
    __m512i q = _mm512_srlv_epi64(_mm512_mul_epu32(_mm512_srlv_epi64(p, vs), vm), vt);
    p = _mm512_sub_epi64(p, _mm512_mul_epu32(q, vn));
    if (two)
        p = _nv_zmm_csub(p, v2n);
    return lazy ? p : _nv_zmm_csub(p, vn);
}

FLINT_FORCE_INLINE void
_mul_halflimb_zmm(nn_ptr r, nn_srcptr a, nn_srcptr b, slong len, nmod_t mod,
                  int s, int t, int two, int lazy)
{
    const __m512i vn = _nv_zmm_set1(mod.n);
    const __m512i v2n = _nv_zmm_set1(2 * mod.n);
    const __m512i vm = _nv_zmm_set1(NMOD_VEC_SIMD_M(mod) >> 32);
    const __m512i vs = _nv_zmm_set1(s);
    const __m512i vt = _nv_zmm_set1(t);
#define STEP(x, y) _mul_halflimb_step_zmm((x), (y), vs, vm, vt, vn, v2n, two, lazy)
    _NV_ZMM_LOOP2(STEP, r, a, b, len);
#undef STEP
}

#if NMOD_VEC_SIMD_IFMA

/*
    3 <= nu <= 50: p = phi 2^52 + plo with two IFMA products, then
    s = floor(p / 2^(nu - 2)) < 2^(nu + 2), q = floor(s w / 2^52) with
    w = m >> 13 (q from m >> 12 halved for nu = 50, where w would not fit
    52 bits), and plo - q n modulo 2^52 as plo + q (2^52 - n), in [0, 2n)
*/
FLINT_FORCE_INLINE __m512i
_mul_ifma50_step_zmm(__m512i x, __m512i y, __m512i vs, __m512i vpow, __m512i vw,
                     __m512i vnneg, __m512i vn, int half, int lazy)
{
    const __m512i zero = _mm512_setzero_si512();
    const __m512i mask = _mm512_set1_epi64((UWORD(1) << 52) - 1);
    __m512i plo = _mm512_madd52lo_epu64(zero, x, y);
    __m512i phi = _mm512_madd52hi_epu64(zero, x, y);
    __m512i s = _mm512_madd52lo_epu64(_mm512_srlv_epi64(plo, vs), phi, vpow);
    __m512i q = _mm512_madd52hi_epu64(zero, s, vw);
    if (half)
        q = _mm512_srli_epi64(q, 1);
    plo = _mm512_and_si512(_mm512_madd52lo_epu64(plo, q, vnneg), mask);
    return lazy ? plo : _nv_zmm_csub(plo, vn);
}

FLINT_FORCE_INLINE void
_mul_ifma50_zmm(nn_ptr r, nn_srcptr a, nn_srcptr b, slong len, nmod_t mod,
                int half, int lazy)
{
    const int nu = FLINT_BITS - mod.norm;
    const __m512i vn = _nv_zmm_set1(mod.n);
    const __m512i vnneg = _nv_zmm_set1((UWORD(1) << 52) - mod.n);
    const __m512i vw = _nv_zmm_set1(NMOD_VEC_SIMD_M(mod) >> (half ? 12 : 13));
    const __m512i vs = _nv_zmm_set1(nu - 2);
    const __m512i vpow = _nv_zmm_set1(UWORD(1) << (54 - nu));
#define STEP(x, y) _mul_ifma50_step_zmm((x), (y), vs, vpow, vw, vnneg, vn, half, lazy)
    _NV_ZMM_LOOP2(STEP, r, a, b, len);
#undef STEP
}

#else

/* 33 <= nu <= 50: in doubles, the product and its error with an fma, the
   quotient rounded, the remainder z in (-n, n) and corrected; lazy: z + n,
   with n folded into the constant 2^52 of the conversion (t52n = 2^52 + n,
   one addition fewer: 1.03-1.11x, Nov 2026) */
FLINT_FORCE_INLINE __m512i
_mul_fp50_step_zmm(__m512i xi, __m512i yi, __m512d vn, __m512d vninv,
                   __m512d t52n, int lazy)
{
    const __m512d t52 = _mm512_set1_pd(0x1.0p52);
    __m512d x = _mm512_sub_pd(_mm512_castsi512_pd(_mm512_or_si512(xi, _mm512_castpd_si512(t52))), t52);
    __m512d y = _mm512_sub_pd(_mm512_castsi512_pd(_mm512_or_si512(yi, _mm512_castpd_si512(t52))), t52);
    __m512d h = _mm512_mul_pd(x, y);
    __m512d q = _mm512_roundscale_pd(_mm512_mul_pd(h, vninv), _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    __m512d l = _mm512_fmsub_pd(x, y, h);
    __m512d z = _mm512_add_pd(_mm512_fnmadd_pd(q, vn, h), l);
    if (lazy)
        return _mm512_castpd_si512(_mm512_xor_pd(_mm512_add_pd(z, t52n), t52));
    z = _mm512_mask_add_pd(z, _mm512_cmp_pd_mask(z, _mm512_setzero_pd(), _CMP_LT_OQ), z, vn);
    return _mm512_castpd_si512(_mm512_xor_pd(_mm512_add_pd(z, t52), t52));
}

FLINT_FORCE_INLINE void
_mul_fp50_zmm(nn_ptr r, nn_srcptr a, nn_srcptr b, slong len, nmod_t mod, int lazy)
{
    const __m512d vn = _mm512_set1_pd((double) mod.n);
    const __m512d vninv = _mm512_set1_pd(1.0 / (double) mod.n);
    const __m512d t52n = _mm512_set1_pd(0x1.0p52 + (double) mod.n);
#define STEP(x, y) _mul_fp50_step_zmm((x), (y), vn, vninv, t52n, lazy)
    _NV_ZMM_LOOP2(STEP, r, a, b, len);
#undef STEP
}

#endif

/*
    51 <= nu <= 62: x = floor(p / 2^(nu - 2)) from the full product p,
    q = floor(h / 2) with h the high word of x m or up to 2 less; with
    VBMI2, the low word by vpmullq and x by vpshrdvq
*/
FLINT_FORCE_INLINE __m512i
_mul_red2bounded_step_zmm(__m512i x, __m512i y, __m512i vs, __m512i vs2, __m512i vm,
                          __m512i vmh, __m512i vn, __m512i v2n, int lazy)
{
    __m512i lo, hi, s, q;
#if NMOD_VEC_SIMD_VBMI2
    lo = _mm512_mullo_epi64(x, y);
    hi = _nv_zmm_mulhi(x, y, _mm512_shuffle_epi32(y, 0xB1));
    s = _mm512_shrdv_epi64(lo, hi, vs);
    (void) vs2;
#else
    hi = _nv_zmm_mul_full(x, y, &lo);
    s = _mm512_or_si512(_mm512_sllv_epi64(hi, vs2), _mm512_srlv_epi64(lo, vs));
#endif
    q = _mm512_srli_epi64(_nv_zmm_mulhi_lazy(s, vm, vmh), 1);
    lo = _nv_zmm_csub(_mm512_sub_epi64(lo, _mm512_mullo_epi64(q, vn)), v2n);
    return lazy ? lo : _nv_zmm_csub(lo, vn);
}

FLINT_FORCE_INLINE void
_mul_red2bounded_zmm(nn_ptr r, nn_srcptr a, nn_srcptr b, slong len, nmod_t mod,
                     int lazy)
{
    const int nu = FLINT_BITS - mod.norm;
    const ulong m = NMOD_VEC_SIMD_M(mod);
    const __m512i vn = _nv_zmm_set1(mod.n);
    const __m512i v2n = _nv_zmm_set1(2 * mod.n);
    const __m512i vm = _nv_zmm_set1(m);
    const __m512i vmh = _nv_zmm_set1(m >> 32);
    const __m512i vs = _nv_zmm_set1(nu - 2);
    const __m512i vs2 = _nv_zmm_set1(66 - nu);
#define STEP(x, y) _mul_red2bounded_step_zmm((x), (y), vs, vs2, vm, vmh, vn, v2n, lazy)
    _NV_ZMM_LOOP2(STEP, r, a, b, len);
#undef STEP
}

/* returns 0 when no kernel applies; lazy: without the last correction */
FLINT_FORCE_INLINE int
_nmod_vec_mul_simd(nn_ptr res, nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod,
                   int lazy)
{
    const int nu = FLINT_BITS - mod.norm;

    if (nu > 62 || nu < 2)
        return 0;
    else if (nu <= 32)
    {
        int s, t, two;
        _HALFLIMB_SHIFTS(s, t, two, nu);
        if (two)
            _mul_halflimb_zmm(res, vec1, vec2, len, mod, s, t, 1, lazy);
        else
            _mul_halflimb_zmm(res, vec1, vec2, len, mod, s, t, 0, lazy);
    }
    else if (nu <= 50)
#if NMOD_VEC_SIMD_IFMA
    {
        if (nu == 50)
            _mul_ifma50_zmm(res, vec1, vec2, len, mod, 1, lazy);
        else
            _mul_ifma50_zmm(res, vec1, vec2, len, mod, 0, lazy);
    }
#else
        _mul_fp50_zmm(res, vec1, vec2, len, mod, lazy);
#endif
    else
        _mul_red2bounded_zmm(res, vec1, vec2, len, mod, lazy);
    return 1;
}

/* AVX2 *******************************************************************/

#elif NMOD_VEC_SIMD_YMM

FLINT_FORCE_INLINE slong
_mul_halflimb_ymm(nn_ptr r, nn_srcptr a, nn_srcptr b, slong len, nmod_t mod, int s, int t,
                  int two, int lazy)
{
    const __m256i vn = _nv_ymm_set1_32(mod.n);
    const __m256i v2n = _nv_ymm_set1(2 * mod.n);
    const __m256i vm = _nv_ymm_set1_32(NMOD_VEC_SIMD_M(mod) >> 32);
    const __m256i vs = _nv_ymm_set1(s);
    const __m256i vt = _nv_ymm_set1(t);
    slong i;

    for (i = 0; i + 4 <= len; i += 4)
    {
        __m256i p = _mm256_mul_epu32(_nv_ymm_ld(a + i), _nv_ymm_ld(b + i));
        __m256i q = _mm256_srlv_epi64(_mm256_mul_epu32(_mm256_srlv_epi64(p, vs), vm), vt);
        p = _mm256_sub_epi64(p, _mm256_mul_epu32(q, vn));
        if (two)
            p = _nv_ymm_csub_s(p, v2n);
        if (!lazy)
            p = two ? _nv_ymm_csub_s(p, vn) : _nv_ymm_csub_32(p, vn);
        _nv_ymm_st(r + i, p);
    }
    return i;
}

#endif

/* AVX2 and NEON: doubles through machine_vectors.h ***********************/

#if !NMOD_VEC_SIMD_ZMM && (NMOD_VEC_SIMD_YMM || NMOD_VEC_SIMD_NEON)

/* 33 <= nu <= 50, as _mul_fp50_step_zmm; a 4-wide step for 4 to 7 remaining
   entries (with the 8-wide loop alone, 0.94-0.97x at 4 entries on the M4);
   lazy: z + n by _nv_st4_lazy */
FLINT_FORCE_INLINE slong
_mul_fp50(nn_ptr r, nn_srcptr a, nn_srcptr b, slong len, nmod_t mod, int lazy)
{
    const vec8d vn = vec8d_set_d((double) mod.n);
    const vec8d vninv = vec8d_set_d(1.0 / (double) mod.n);
    const vec4d t52n = vec4d_set_d(0x1.0p52 + (double) mod.n);
    slong i;

    for (i = 0; i + 8 <= len; i += 8)
    {
        vec8d x = vec8n_convert_limited_vec8d(vec8n_load_unaligned(a + i));
        vec8d y = vec8n_convert_limited_vec8d(vec8n_load_unaligned(b + i));
        vec8d z = vec8d_mulmod(x, y, vn, vninv);
        if (lazy)
        {
            _nv_st4_lazy(r + i, z.e1, vn.e1, t52n);
            _nv_st4_lazy(r + i + 4, z.e2, vn.e1, t52n);
        }
        else
        {
            z = vec8d_reduce_pm1no_to_0n(z, vn);
            vec4n_store_unaligned(r + i, vec4d_convert_limited_vec4n(z.e1));
            vec4n_store_unaligned(r + i + 4, vec4d_convert_limited_vec4n(z.e2));
        }
    }
    if (i + 4 <= len)
    {
        vec4d x = vec4n_convert_limited_vec4d(vec4n_load_unaligned(a + i));
        vec4d y = vec4n_convert_limited_vec4d(vec4n_load_unaligned(b + i));
        vec4d z = vec4d_mulmod(x, y, vn.e1, vninv.e1);
        if (lazy)
            _nv_st4_lazy(r + i, z, vn.e1, t52n);
        else
            vec4n_store_unaligned(r + i, vec4d_convert_limited_vec4n(
                                      vec4d_reduce_pm1no_to_0n(z, vn.e1)));
        i += 4;
    }
    return i;
}

#endif

#if NMOD_VEC_SIMD_YMM && !NMOD_VEC_SIMD_ZMM

FLINT_FORCE_INLINE int
_nmod_vec_mul_simd(nn_ptr res, nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod,
                   int lazy)
{
    const int nu = FLINT_BITS - mod.norm;
    slong i;

    if (nu > 50 || nu < 2)
        return 0;
    else if (nu <= 32)
    {
        int s, t, two;
        _HALFLIMB_SHIFTS(s, t, two, nu);
        if (two)
            i = _mul_halflimb_ymm(res, vec1, vec2, len, mod, s, t, 1, lazy);
        else
            i = _mul_halflimb_ymm(res, vec1, vec2, len, mod, s, t, 0, lazy);
    }
    else
        i = _mul_fp50(res, vec1, vec2, len, mod, lazy);
    _mul_tail(res, vec1, vec2, i, len, mod);
    return 1;
}

/* NEON *******************************************************************/

#elif NMOD_VEC_SIMD_NEON

/* nu <= 32, four entries per step */
#define _MUL_HALFLIMB_STEP_NEON(p, x, y) \
    do { \
        (p) = vmull_u32(vmovn_u64(x), vmovn_u64(y)); \
        (p) = vmlsl_u32((p), vmovn_u64(vshlq_u64(vmull_u32(vmovn_u64(vshlq_u64((p), vs)), vm), vt)), vnl); \
        if (two) \
            (p) = _nv_neon_csub((p), v2n); \
        if (!lazy) \
            (p) = _nv_neon_csub((p), vn); \
    } while (0)

#define _MUL_HALFLIMB_CONSTS_NEON \
    const uint64x2_t vn = vdupq_n_u64(mod.n); \
    const uint64x2_t v2n = vdupq_n_u64(2 * mod.n); \
    const uint32x2_t vnl = vdup_n_u32((uint32_t) mod.n); \
    const uint32x2_t vm = vdup_n_u32((uint32_t) (NMOD_VEC_SIMD_M(mod) >> 32)); \
    const int64x2_t vs = vdupq_n_s64(-(int64_t) s); \
    const int64x2_t vt = vdupq_n_s64(-(int64_t) t)

FLINT_FORCE_INLINE slong
_mul_halflimb_neon(nn_ptr r, nn_srcptr a, nn_srcptr b, slong len, nmod_t mod, int s, int t,
                   int two, int lazy)
{
    _MUL_HALFLIMB_CONSTS_NEON;
    slong i;

    for (i = 0; i + 4 <= len; i += 4)
    {
        uint64x2_t p0, p1;
        _MUL_HALFLIMB_STEP_NEON(p0, _nv_neon_ld(a + i), _nv_neon_ld(b + i));
        _MUL_HALFLIMB_STEP_NEON(p1, _nv_neon_ld(a + i + 2), _nv_neon_ld(b + i + 2));
        _nv_neon_st(r + i, p0);
        _nv_neon_st(r + i + 2, p1);
    }
    return i;
}

/* the same, by blocks of 8: 4 entries in NEON and 4 with NMOD_RED (the
   product fits one word) */
FLINT_FORCE_INLINE slong
_mul_halflimb_hybrid_neon(nn_ptr r, nn_srcptr a, nn_srcptr b, slong len, nmod_t mod, int s,
                          int t, int two, int lazy)
{
    _MUL_HALFLIMB_CONSTS_NEON;
    slong i, k;

    for (i = 0; i + 8 <= len; i += 8)
    {
        uint64x2_t p0, p1;
        _MUL_HALFLIMB_STEP_NEON(p0, _nv_neon_ld(a + i), _nv_neon_ld(b + i));
        _MUL_HALFLIMB_STEP_NEON(p1, _nv_neon_ld(a + i + 2), _nv_neon_ld(b + i + 2));
        _nv_neon_st(r + i, p0);
        _nv_neon_st(r + i + 2, p1);
        for (k = 4; k < 8; k++)
            NMOD_RED(r[i + k], a[i + k] * b[i + k], mod);
    }
    return i;
}

/* the hybrid loop from NMOD_VEC_SIMD_HYBRID_MIN entries, then the NEON one */
FLINT_FORCE_INLINE slong
_mul_halflimb_neon_all(nn_ptr r, nn_srcptr a, nn_srcptr b, slong len, nmod_t mod, int s,
                       int t, int two, int lazy)
{
    slong i = 0;
    if (len >= NMOD_VEC_SIMD_HYBRID_MIN)
        i = _mul_halflimb_hybrid_neon(r, a, b, len, mod, s, t, two, lazy);
    return i + _mul_halflimb_neon(r + i, a + i, b + i, len - i, mod, s, t, two, lazy);
}

/* above 50 bits the scalar code (the bounded loop of #2950 up to 62 bits,
   1.1-1.2x faster than nmod_mul on the M4) */
FLINT_FORCE_INLINE int
_nmod_vec_mul_simd(nn_ptr res, nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod,
                   int lazy)
{
    const int nu = FLINT_BITS - mod.norm;
    slong i;

    if (nu > 50 || nu < 2)
        return 0;
    else if (nu <= 32)
    {
        int s, t, two;
        _HALFLIMB_SHIFTS(s, t, two, nu);
        if (two)
            i = _mul_halflimb_neon_all(res, vec1, vec2, len, mod, s, t, 1, lazy);
        else
            i = _mul_halflimb_neon_all(res, vec1, vec2, len, mod, s, t, 0, lazy);
    }
    else
        i = _mul_fp50(res, vec1, vec2, len, mod, lazy);
    _mul_tail(res, vec1, vec2, i, len, mod);
    return 1;
}

#endif

void _nmod_vec_mul(nn_ptr res, nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
{
#if NMOD_VEC_SIMD_ZMM || NMOD_VEC_SIMD_YMM || NMOD_VEC_SIMD_NEON
    if (len >= NMOD_VEC_SIMD_MIN_LEN && _nmod_vec_mul_simd(res, vec1, vec2, len, mod, 0))
        return;
#endif
    _nmod_vec_mul_scalar(res, vec1, vec2, len, mod);
}

/* the scalar code reduces fully, which is also in [0, 2n) */
void _nmod_vec_mul_lazy(nn_ptr res, nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
{
    FLINT_ASSERT(NMOD_BITS(mod) <= 62);
#if NMOD_VEC_SIMD_ZMM || NMOD_VEC_SIMD_YMM || NMOD_VEC_SIMD_NEON
    if (len >= NMOD_VEC_SIMD_MIN_LEN && _nmod_vec_mul_simd(res, vec1, vec2, len, mod, 1))
        return;
#endif
    _nmod_vec_mul_scalar(res, vec1, vec2, len, mod);
}
