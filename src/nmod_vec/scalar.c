/*
    Copyright (C) 2010, 2015 William Hart
    Copyright (C) 2015 Vladimir Glazachev
    Copyright (C) 2021 Fredrik Johansson
    Copyright (C) 2024, 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint-mparam.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "impl_simd.h"

/*
    SIMD kernels of the products by a scalar c (not 0, 1 or n - 1), by the
    size nu of n, with c_pr = floor(c 2^64 / n) (doc: nmod_vec/simd#scalar,
    and doc: algorithms/arithmetic/shoup#lem-quot for the radixes 2^32 and
    2^52):
      shoup32 (nu <= 32): Shoup's method in radix 2^32 with c_pr >> 32,
        every product 32 x 32 bits, one correction.
      shoup52 (33 <= nu <= 51, AVX-512 with IFMA): Shoup's method in radix
        2^52 with c_pr >> 12, three IFMA, a mask, one correction (Nov 2026:
        1.2-3.4x faster than fp50 and shoup64lazy on Ice Lake and Zen 4).
      fp50 (33 <= nu <= 50, otherwise): doubles (ahead of shoup64lazy on
        Cascade Lake and Ice Lake, even on Zen 4: Nov 2026).
      shoup64lazy (51 or 52 <= nu <= 62, AVX-512): Shoup with the high word
        of x c_pr from three 32-bit products (2 too small at most),
        corrections 2n then n (5-20% faster than with the exact high word).
      shoup64 (nu = 63, AVX-512): Shoup with the exact high word.
    For nu = 64, the scalar code (Montgomery for odd n). The sum of
    scalar_addmul is corrected once more. The lazy variants (n < 2^(FLINT_BITS - 2)) skip
    the last correction of the product (fp50: z + n), and scalar_addmul_lazy
    adds res in [0, 2n) and corrects the sum, below 4n, by 2n: 1.08-1.60x for
    scalar_mul and 1.04-1.35x for scalar_addmul on long vectors, 1.10-1.21x
    for shoup52 (Nov 2026). Measured and not kept (Nov 2026):
    on AVX2, Shoup above 32 bits with emulated 64-bit products (0.8x on
    Broadwell); on NEON, Shoup above 32 bits (0.48-0.82x) and fp50 for
    scalar_addmul (0.89-0.96x on long vectors).
*/

#if !NMOD_VEC_SIMD_ZMM && (NMOD_VEC_SIMD_YMM || NMOD_VEC_SIMD_NEON)
FLINT_FORCE_INLINE void
_smul_tail(nn_ptr res, nn_srcptr vec, slong i, slong len, ulong c, ulong c_pr,
           nmod_t mod, int add, int lazy)
{
    NMOD_VEC_NO_VECTORIZE
    for ( ; i < len; i++)
    {
        ulong t = n_mulmod_shoup(c, vec[i], c_pr, mod.n);
        if (add && lazy)
        {
            /* res[i] in [0, 2n): the sum is below 3n */
            t += res[i];
            res[i] = (t >= 2 * mod.n) ? t - 2 * mod.n : t;
        }
        else
            res[i] = add ? nmod_add(res[i], t, mod) : t;
    }
}
#endif

#if NMOD_VEC_SIMD_ZMM

FLINT_FORCE_INLINE __m512i
_smul_shoup32_step_zmm(__m512i x, __m512i vc, __m512i vcp, __m512i vn, int lazy)
{
    __m512i q = _mm512_srli_epi64(_mm512_mul_epu32(x, vcp), 32);
    __m512i t = _mm512_sub_epi64(_mm512_mul_epu32(x, vc), _mm512_mul_epu32(q, vn));
    return lazy ? t : _nv_zmm_csub(t, vn);
}

/* lazy, here and below: the products without their last correction, and
   for addmul the sum (below 4n) corrected by 2n */
FLINT_FORCE_INLINE void
_smul_shoup32_zmm(nn_ptr r, nn_srcptr a, slong len, ulong c, nmod_t mod, int add,
                  int lazy)
{
    const ulong c_pr = _nmod_mulmod_precomp_shoup(c, mod);
    const __m512i vn = _nv_zmm_set1(mod.n);
    const __m512i v2n = _nv_zmm_set1(2 * mod.n);
    const __m512i vc = _nv_zmm_set1(c);
    const __m512i vcp = _nv_zmm_set1(c_pr >> 32);
#define STEP(x) _smul_shoup32_step_zmm((x), vc, vcp, vn, lazy)
    _NV_ZMM_SLOOP(STEP, r, a, len, add, lazy ? v2n : vn);
#undef STEP
}

/* 52 <= nu <= 62 (51 without IFMA): the high word of x c_pr from three
   32-bit products, up to 2 too small, so x c - q n < 4n: corrections 2n,
   then n */
FLINT_FORCE_INLINE __m512i
_smul_shoup64lazy_step_zmm(__m512i x, __m512i vc, __m512i vcp, __m512i vcph,
                           __m512i vn, __m512i v2n, int lazy)
{
    __m512i q = _nv_zmm_mulhi_lazy(x, vcp, vcph);
    __m512i t = _mm512_sub_epi64(_mm512_mullo_epi64(x, vc), _mm512_mullo_epi64(q, vn));
    t = _nv_zmm_csub(t, v2n);
    return lazy ? t : _nv_zmm_csub(t, vn);
}

FLINT_FORCE_INLINE void
_smul_shoup64lazy_zmm(nn_ptr r, nn_srcptr a, slong len, ulong c, nmod_t mod, int add,
                      int lazy)
{
    const ulong c_pr = _nmod_mulmod_precomp_shoup(c, mod);
    const __m512i vn = _nv_zmm_set1(mod.n);
    const __m512i v2n = _nv_zmm_set1(2 * mod.n);
    const __m512i vc = _nv_zmm_set1(c);
    const __m512i vcp = _nv_zmm_set1(c_pr);
    const __m512i vcph = _nv_zmm_set1(c_pr >> 32);
#define STEP(x) _smul_shoup64lazy_step_zmm((x), vc, vcp, vcph, vn, v2n, lazy)
    _NV_ZMM_SLOOP(STEP, r, a, len, add, lazy ? v2n : vn);
#undef STEP
}

/* nu = 63: the exact high word of x c_pr, one correction */
FLINT_FORCE_INLINE __m512i
_smul_shoup64_step_zmm(__m512i x, __m512i vc, __m512i vcp, __m512i vcph, __m512i vn)
{
    __m512i q = _nv_zmm_mulhi(x, vcp, vcph);
    __m512i t = _mm512_sub_epi64(_mm512_mullo_epi64(x, vc), _mm512_mullo_epi64(q, vn));
    return _nv_zmm_csub(t, vn);
}

FLINT_FORCE_INLINE void
_smul_shoup64_zmm(nn_ptr r, nn_srcptr a, slong len, ulong c, nmod_t mod, int add)
{
    const ulong c_pr = _nmod_mulmod_precomp_shoup(c, mod);
    const __m512i vn = _nv_zmm_set1(mod.n);
    const __m512i vc = _nv_zmm_set1(c);
    const __m512i vcp = _nv_zmm_set1(c_pr);
    const __m512i vcph = _nv_zmm_set1(c_pr >> 32);
#define STEP(x) _smul_shoup64_step_zmm((x), vc, vcp, vcph, vn)
    _NV_ZMM_SLOOP(STEP, r, a, len, add, vn);
#undef STEP
}

#if NMOD_VEC_SIMD_IFMA

/* 33 <= nu <= 51, IFMA: with c' = c_pr >> 12 = floor(c 2^52 / n), q =
   floor(x c' / 2^52) is floor(x c / n) or one less, so x c - q n is in
   [0, 2n), below 2^52: computed as (x c + q (2^52 - n)) mod 2^52 */
FLINT_FORCE_INLINE __m512i
_smul_shoup52_step_zmm(__m512i x, __m512i vc, __m512i vcp, __m512i vnneg, __m512i vn,
                       int lazy)
{
    const __m512i zero = _mm512_setzero_si512();
    const __m512i mask = _mm512_set1_epi64((UWORD(1) << 52) - 1);
    __m512i q = _mm512_madd52hi_epu64(zero, x, vcp);
    __m512i t = _mm512_madd52lo_epu64(_mm512_madd52lo_epu64(zero, x, vc), q, vnneg);
    t = _mm512_and_si512(t, mask);
    return lazy ? t : _nv_zmm_csub(t, vn);
}

FLINT_FORCE_INLINE void
_smul_shoup52_zmm(nn_ptr r, nn_srcptr a, slong len, ulong c, nmod_t mod, int add,
                  int lazy)
{
    const ulong c_pr = _nmod_mulmod_precomp_shoup(c, mod);
    const __m512i vn = _nv_zmm_set1(mod.n);
    const __m512i v2n = _nv_zmm_set1(2 * mod.n);
    const __m512i vc = _nv_zmm_set1(c);
    const __m512i vcp = _nv_zmm_set1(c_pr >> 12);
    const __m512i vnneg = _nv_zmm_set1((UWORD(1) << 52) - mod.n);
#define STEP(x) _smul_shoup52_step_zmm((x), vc, vcp, vnneg, vn, lazy)
    _NV_ZMM_SLOOP(STEP, r, a, len, add, lazy ? v2n : vn);
#undef STEP
}

#else

/* 33 <= nu <= 50: doubles, as _mul_fp50_step_zmm (lazy: z + n) */
FLINT_FORCE_INLINE __m512i
_smul_fp50_step_zmm(__m512i xi, __m512d vcd, __m512d vn, __m512d vninv,
                    __m512d t52n, int lazy)
{
    const __m512d t52 = _mm512_set1_pd(0x1.0p52);
    __m512d x = _mm512_sub_pd(_mm512_castsi512_pd(_mm512_or_si512(xi, _mm512_castpd_si512(t52))), t52);
    __m512d h = _mm512_mul_pd(x, vcd);
    __m512d q = _mm512_roundscale_pd(_mm512_mul_pd(h, vninv), _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    __m512d l = _mm512_fmsub_pd(x, vcd, h);
    __m512d z = _mm512_add_pd(_mm512_fnmadd_pd(q, vn, h), l);
    if (lazy)
        return _mm512_castpd_si512(_mm512_xor_pd(_mm512_add_pd(z, t52n), t52));
    z = _mm512_mask_add_pd(z, _mm512_cmp_pd_mask(z, _mm512_setzero_pd(), _CMP_LT_OQ), z, vn);
    return _mm512_castpd_si512(_mm512_xor_pd(_mm512_add_pd(z, t52), t52));
}

FLINT_FORCE_INLINE void
_smul_fp50_zmm(nn_ptr r, nn_srcptr a, slong len, ulong c, nmod_t mod, int add, int lazy)
{
    const __m512d vn = _mm512_set1_pd((double) mod.n);
    const __m512d vninv = _mm512_set1_pd(1.0 / (double) mod.n);
    const __m512d vcd = _mm512_set1_pd((double) c);
    const __m512d t52n = _mm512_set1_pd(0x1.0p52 + (double) mod.n);
    const __m512i vnz = _nv_zmm_set1(mod.n);
    const __m512i v2nz = _nv_zmm_set1(2 * mod.n);
#define STEP(x) _smul_fp50_step_zmm((x), vcd, vn, vninv, t52n, lazy)
    _NV_ZMM_SLOOP(STEP, r, a, len, add, lazy ? v2nz : vnz);
#undef STEP
}

#endif

/* returns 0 when no kernel applies; lazy (nu <= 62): outputs in [0, 2n),
   and for addmul res in [0, 2n) on input */
FLINT_FORCE_INLINE int
_nmod_vec_scalar_mul_simd(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod,
                          int add, int lazy)
{
    const int nu = FLINT_BITS - mod.norm;

    if (nu > NMOD_VEC_SIMD_SCALAR_MAX_BITS(add) || nu < 2)
        return 0;
    else if (nu <= 32)
        _smul_shoup32_zmm(res, vec, len, c, mod, add, lazy);
#if NMOD_VEC_SIMD_IFMA
    else if (nu <= 51)
        _smul_shoup52_zmm(res, vec, len, c, mod, add, lazy);
#else
    else if (nu <= 50)
        _smul_fp50_zmm(res, vec, len, c, mod, add, lazy);
#endif
    else if (nu <= 62)
        _smul_shoup64lazy_zmm(res, vec, len, c, mod, add, lazy);
    else
        _smul_shoup64_zmm(res, vec, len, c, mod, add);
    return 1;
}

#elif NMOD_VEC_SIMD_YMM

/* nu <= 32; c32: nu <= 31, the corrections on 32-bit lanes; lazy: the
   product without its correction, and for addmul the sum (below 4n)
   corrected by 2n on 64-bit lanes */
FLINT_FORCE_INLINE slong
_smul_shoup32_ymm(nn_ptr r, nn_srcptr a, slong len, ulong c, ulong c_pr, nmod_t mod,
                  int add, int c32, int lazy)
{
    const __m256i vn = _nv_ymm_set1_32(mod.n);
    const __m256i v2n = _nv_ymm_set1(2 * mod.n);
    const __m256i vc = _nv_ymm_set1_32(c);
    const __m256i vcp = _nv_ymm_set1_32(c_pr >> 32);
    slong i;

    for (i = 0; i + 4 <= len; i += 4)
    {
        __m256i x = _nv_ymm_ld(a + i);
        __m256i q = _mm256_srli_epi64(_mm256_mul_epu32(x, vcp), 32);
        __m256i t = _mm256_sub_epi64(_mm256_mul_epu32(x, vc), _mm256_mul_epu32(q, vn));
        if (!lazy)
            t = c32 ? _nv_ymm_csub_32(t, vn) : _nv_ymm_csub_s(t, vn);
        if (add)
        {
            t = _mm256_add_epi64(t, _nv_ymm_ld(r + i));
            if (lazy)
                t = _nv_ymm_csub_s(t, v2n);
            else
                t = c32 ? _nv_ymm_csub_32(t, vn) : _nv_ymm_csub_s(t, vn);
        }
        _nv_ymm_st(r + i, t);
    }
    return i;
}

#elif NMOD_VEC_SIMD_NEON

/* nu <= 32; lazy as in _smul_shoup32_ymm */
FLINT_FORCE_INLINE slong
_smul_shoup32_neon(nn_ptr r, nn_srcptr a, slong len, ulong c, ulong c_pr, nmod_t mod,
                   int add, int lazy)
{
    const uint64x2_t vn = vdupq_n_u64(mod.n);
    const uint64x2_t vs = vdupq_n_u64(lazy ? 2 * mod.n : mod.n);
    const uint32x2_t vnl = vdup_n_u32((uint32_t) mod.n);
    const uint32x2_t vc = vdup_n_u32((uint32_t) c);
    const uint32x2_t vcp = vdup_n_u32((uint32_t) (c_pr >> 32));
    slong i;

    for (i = 0; i + 4 <= len; i += 4)
    {
        uint32x2_t x0 = vmovn_u64(_nv_neon_ld(a + i)), x1 = vmovn_u64(_nv_neon_ld(a + i + 2));
        uint32x2_t q0 = vshrn_n_u64(vmull_u32(x0, vcp), 32);
        uint32x2_t q1 = vshrn_n_u64(vmull_u32(x1, vcp), 32);
        uint64x2_t t0 = vmlsl_u32(vmull_u32(x0, vc), q0, vnl);
        uint64x2_t t1 = vmlsl_u32(vmull_u32(x1, vc), q1, vnl);
        if (!lazy)
        {
            t0 = _nv_neon_csub(t0, vn);
            t1 = _nv_neon_csub(t1, vn);
        }
        if (add)
        {
            t0 = _nv_neon_csub(vaddq_u64(t0, _nv_neon_ld(r + i)), vs);
            t1 = _nv_neon_csub(vaddq_u64(t1, _nv_neon_ld(r + i + 2)), vs);
        }
        _nv_neon_st(r + i, t0);
        _nv_neon_st(r + i + 2, t1);
    }
    return i;
}

#endif

#if !NMOD_VEC_SIMD_ZMM && (NMOD_VEC_SIMD_YMM || NMOD_VEC_SIMD_NEON)

/* 33 <= nu <= 50: doubles, with a 4-wide step for 4 to 7 remaining
   entries; for addmul, z + s - n is in [-n, n). Lazy: z + n by _nv_st4_lazy,
   and for addmul z + s - n in (-2n, 2n), with z in (-n, n) not corrected and
   s in [0, 2n), brought to [0, 2n) */
FLINT_FORCE_INLINE slong
_smul_fp50(nn_ptr r, nn_srcptr a, slong len, ulong c, nmod_t mod, int add, int lazy)
{
    const vec8d vn = vec8d_set_d((double) mod.n);
    const vec8d v2n = vec8d_set_d(2.0 * (double) mod.n);
    const vec8d vninv = vec8d_set_d(1.0 / (double) mod.n);
    const vec8d vc = vec8d_set_d((double) c);
    const vec4d t52n = vec4d_set_d(0x1.0p52 + (double) mod.n);
    slong i;

    for (i = 0; i + 8 <= len; i += 8)
    {
        vec8d x = vec8n_convert_limited_vec8d(vec8n_load_unaligned(a + i));
        vec8d z = vec8d_mulmod(x, vc, vn, vninv);
        if (lazy && !add)
        {
            _nv_st4_lazy(r + i, z.e1, vn.e1, t52n);
            _nv_st4_lazy(r + i + 4, z.e2, vn.e1, t52n);
            continue;
        }
        if (!lazy)
            z = vec8d_reduce_pm1no_to_0n(z, vn);
        if (add)
        {
            vec8d s = vec8n_convert_limited_vec8d(vec8n_load_unaligned(r + i));
            z = vec8d_reduce_pm1no_to_0n(vec8d_sub(vec8d_add(z, s), vn), lazy ? v2n : vn);
        }
        vec4n_store_unaligned(r + i, vec4d_convert_limited_vec4n(z.e1));
        vec4n_store_unaligned(r + i + 4, vec4d_convert_limited_vec4n(z.e2));
    }
    if (i + 4 <= len)
    {
        vec4d x = vec4n_convert_limited_vec4d(vec4n_load_unaligned(a + i));
        vec4d z = vec4d_mulmod(x, vc.e1, vn.e1, vninv.e1);
        if (lazy && !add)
            _nv_st4_lazy(r + i, z, vn.e1, t52n);
        else
        {
            if (!lazy)
                z = vec4d_reduce_pm1no_to_0n(z, vn.e1);
            if (add)
            {
                vec4d s = vec4n_convert_limited_vec4d(vec4n_load_unaligned(r + i));
                z = vec4d_reduce_pm1no_to_0n(vec4d_sub(vec4d_add(z, s), vn.e1),
                                             lazy ? v2n.e1 : vn.e1);
            }
            vec4n_store_unaligned(r + i, vec4d_convert_limited_vec4n(z));
        }
        i += 4;
    }
    return i;
}

FLINT_FORCE_INLINE int
_nmod_vec_scalar_mul_simd(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod,
                          int add, int lazy)
{
    const int nu = FLINT_BITS - mod.norm;
    ulong c_pr;
    slong i;

    if (nu > NMOD_VEC_SIMD_SCALAR_MAX_BITS(add) || nu < 2)
        return 0;
    c_pr = _nmod_mulmod_precomp_shoup(c, mod);
#if NMOD_VEC_SIMD_YMM
    i = _nv_ymm_peel(res, vec, vec, len);
    if (i != 0)
    {
        _smul_tail(res, vec, 0, i, c, c_pr, mod, add, lazy);
        res += i;
        vec += i;
        len -= i;
    }
#endif
    if (nu <= 32)
#if NMOD_VEC_SIMD_YMM
        i = (nu <= 31) ? _smul_shoup32_ymm(res, vec, len, c, c_pr, mod, add, 1, lazy)
                       : _smul_shoup32_ymm(res, vec, len, c, c_pr, mod, add, 0, lazy);
#else
        i = _smul_shoup32_neon(res, vec, len, c, c_pr, mod, add, lazy);
#endif
    else
        i = _smul_fp50(res, vec, len, c, mod, add, lazy);
    _smul_tail(res, vec, i, len, c, c_pr, mod, add, lazy);
    return 1;
}

#endif

#if NMOD_VEC_SIMD_ZMM || NMOD_VEC_SIMD_YMM || NMOD_VEC_SIMD_NEON

/* out of line (impl_simd.h); n has 2 to NMOD_VEC_SIMD_SCALAR_MAX_BITS(add) bits */
FLINT_STATIC_NOINLINE void
_scalar_mul_simd(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)
{
    _nmod_vec_scalar_mul_simd(res, vec, len, c, mod, 0, 0);
}

FLINT_STATIC_NOINLINE void
_scalar_addmul_simd(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)
{
    _nmod_vec_scalar_mul_simd(res, vec, len, c, mod, 1, 0);
}

FLINT_STATIC_NOINLINE void
_scalar_mul_lazy_simd(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)
{
    _nmod_vec_scalar_mul_simd(res, vec, len, c, mod, 0, 1);
}

FLINT_STATIC_NOINLINE void
_scalar_addmul_lazy_simd(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)
{
    _nmod_vec_scalar_mul_simd(res, vec, len, c, mod, 1, 1);
}

/* n > 1 holds there, since c is not 0 */
# define _SCALAR_SIMD_APPLIES(len, mod, add) ((len) >= NMOD_VEC_SIMD_MIN_LEN \
        && NMOD_BITS(mod) <= NMOD_VEC_SIMD_SCALAR_MAX_BITS(add))

#endif

void
_nmod_vec_scalar_addmul_nmod_generic(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)
{
    slong i;

    for (i = 0; i < len; i++)
        res[i] = nmod_addmul(res[i], vec[i], c, mod);
}

void _nmod_vec_scalar_addmul_nmod_shoup(nn_ptr res, nn_srcptr vec,
				             slong len, ulong c, nmod_t mod)
{
    const ulong c_pr = _nmod_mulmod_precomp_shoup(c, mod);

    for (slong i = 0; i < len; i++)
    {
        ulong t = n_mulmod_shoup(c, vec[i], c_pr, mod.n);
        res[i] = _nmod_add(res[i], t, mod);
    }
}

/* the generic version for n < 2^(FLINT_BITS / 2), where res[i] + vec[i] c
   < n^2 fits one word, and for NMOD_CAN_USE_MUL_BOUNDED(mod); the dispatcher
   tests the modulus, so that the generic version keeps its code */
static void
_scalar_addmul_halflimb(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)
{
    slong i;

    /* clang vectorizes this loop into slower code at short lengths
       (measured on x86) */
#if defined(__clang__)
#pragma clang loop vectorize(disable)
#endif
    for (i = 0; i < len; i++)
        NMOD_RED(res[i], res[i] + vec[i] * c, mod);
}

static void
_scalar_addmul_bounded(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)
{
    slong i;

    for (i = 0; i < len; i++)
        res[i] = _nmod_addmul_bounded(res[i], vec[i], c, mod);
}

/* For n < 2^(FLINT_BITS / 2), the one reduction of res[i] + vec[i] c is
   faster than Shoup's product and addition. Built with clang, Shoup's loop is
   vectorized and becomes faster from 32 terms on (on Ice Lake; from 256 on
   Apple M4), so that the one-word loop goes first only below 32 terms. */
#if defined(__clang__)
#define HALFLIMB_BEFORE_SHOUP(len) ((len) < 32)
#else
#define HALFLIMB_BEFORE_SHOUP(len) 1
#endif

void _nmod_vec_scalar_addmul_nmod(nn_ptr res, nn_srcptr vec,
				             slong len, ulong c, nmod_t mod)
{
    if (c == UWORD(0))
        return;
    else if (c == UWORD(1))
        _nmod_vec_add(res, res, vec, len, mod);
    else if (c == mod.n - UWORD(1))
        _nmod_vec_sub(res, res, vec, len, mod);
#if NMOD_VEC_SIMD_ZMM || NMOD_VEC_SIMD_YMM || NMOD_VEC_SIMD_NEON
    else if (_SCALAR_SIMD_APPLIES(len, mod, 1))
        _scalar_addmul_simd(res, vec, len, c, mod);
#endif
    else if (NMOD_BITS(mod) == FLINT_BITS)
        _nmod_vec_scalar_addmul_nmod_generic(res, vec, len, c, mod);
    else if (mod.norm >= FLINT_BITS / 2 && HALFLIMB_BEFORE_SHOUP(len))
        _scalar_addmul_halflimb(res, vec, len, c, mod);
    else if (len >= FLINT_MULMOD_SHOUP_THRESHOLD)
        _nmod_vec_scalar_addmul_nmod_shoup(res, vec, len, c, mod);
    else if (mod.norm >= FLINT_BITS / 2)
        _scalar_addmul_halflimb(res, vec, len, c, mod);
    else if (NMOD_CAN_USE_MUL_BOUNDED(mod))
        _scalar_addmul_bounded(res, vec, len, c, mod);
    else
        _nmod_vec_scalar_addmul_nmod_generic(res, vec, len, c, mod);
}

void
_nmod_vec_scalar_mul_nmod_redc(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)
{
    nmod_redc_ctx_t ctx;
    slong i;
    ulong c_redc;
    nmod_redc_ctx_init_nmod(ctx, mod);
    c_redc = nmod_redc_set_nmod(c, ctx);
    for (i = 0; i < len; i++)
        res[i] = nmod_redc_mul(vec[i], c_redc, ctx);
}

static void _nmod_vec_scalar_mul_nmod_fullword(nn_ptr res, nn_srcptr vec,
                               slong len, ulong c, nmod_t mod)
{
    slong i;

    for (i = 0; i < len; i++)
        NMOD_MUL_FULLWORD(res[i], vec[i], c, mod);
}

void _nmod_vec_scalar_mul_nmod_generic(nn_ptr res, nn_srcptr vec,
                               slong len, ulong c, nmod_t mod)
{
    slong i;

    for (i = 0; i < len; i++)
        NMOD_MUL_PRENORM(res[i], vec[i], c << mod.norm, mod);
}

void _nmod_vec_scalar_mul_nmod_shoup(nn_ptr res, nn_srcptr vec,
                               slong len, ulong c, nmod_t mod)
{
    const ulong c_pr = _nmod_mulmod_precomp_shoup(c, mod);

    for (slong i = 0; i < len; i++)
        res[i] = n_mulmod_shoup(c, vec[i], c_pr, mod.n);
}

/* as _scalar_addmul_halflimb and _scalar_addmul_bounded; clang vectorizes
   the bounded loop into code slower than the scalar loop (measured on x86) */
static void
_scalar_mul_halflimb(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)
{
    slong i;

    for (i = 0; i < len; i++)
        res[i] = _nmod_mul_halflimb(vec[i], c, mod);
}

static void
_scalar_mul_bounded(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)
{
    slong i;

#if defined(__clang__)
#pragma clang loop vectorize(disable)
#endif
    for (i = 0; i < len; i++)
        res[i] = _nmod_mul_bounded(vec[i], c, mod);
}

void _nmod_vec_scalar_mul_nmod(nn_ptr res, nn_srcptr vec,
                               slong len, ulong c, nmod_t mod)
{
    if (c == UWORD(0))
        _nmod_vec_zero(res, len);
    else if (c == UWORD(1))
        _nmod_vec_set(res, vec, len);
    else if (c == mod.n - UWORD(1))
        _nmod_vec_neg(res, vec, len, mod);
#if NMOD_VEC_SIMD_ZMM || NMOD_VEC_SIMD_YMM || NMOD_VEC_SIMD_NEON
    else if (_SCALAR_SIMD_APPLIES(len, mod, 0))
        _scalar_mul_simd(res, vec, len, c, mod);
#endif
    else if (NMOD_BITS(mod) == FLINT_BITS)
    {
        if (len >= 8 && mod.n % 2 != 0)
            _nmod_vec_scalar_mul_nmod_redc(res, vec, len, c, mod);
        else
            _nmod_vec_scalar_mul_nmod_fullword(res, vec, len, c, mod);
    }
    else if (len >= FLINT_MULMOD_SHOUP_THRESHOLD)
        _nmod_vec_scalar_mul_nmod_shoup(res, vec, len, c, mod);
    else if (NMOD_IS_HALFLIMB(mod))
        _scalar_mul_halflimb(res, vec, len, c, mod);
    else if (NMOD_CAN_USE_MUL_BOUNDED(mod))
        _scalar_mul_bounded(res, vec, len, c, mod);
    else
        _nmod_vec_scalar_mul_nmod_generic(res, vec, len, c, mod);
}

/*
    The lazy variants. Below the kernels' length, and where no kernel
    applies (no SIMD, AVX2 above 50 bits, NEON above 50 bits and addmul
    above 32 bits):
      above 32 bits, from 8 entries, Shoup's product without its last
      correction, in [0, 2n) (doc: algorithms/arithmetic/shoup#lem-quot),
      and for addmul the sum with res, below 4n, corrected by 2n. Against
      the exact code: 1.05-1.66x at 8-12 entries (0.92 once, Zen 4 with
      AVX2 only), 1.06-1.53x from 16 (Nov 2026).
      scalar_addmul, n < 2^32: the one-word loop of the exact code, since
      res + vec c < 2n + (n - 1)^2 <= n^2 fits one word for res in [0, 2n)
      (with a reduction of res first: 0.56-0.77x the exact code).
      Otherwise the exact code, whose outputs are in [0, 2n) too;
      scalar_addmul first brings res to [0, n) there.
    The cases c = 1 and c = n - 1 of scalar_addmul go to the kernels or to
    the loops above, since _nmod_vec_add and _nmod_vec_sub need reduced
    inputs.
*/
#define NMOD_VEC_LAZY_SHOUP_MIN_LEN 8

static void
_scalar_mul_lazy_shoup(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)
{
    const ulong c_pr = _nmod_mulmod_precomp_shoup(c, mod), n = mod.n;
    slong i;

    for (i = 0; i < len; i++)
        res[i] = c * vec[i] - n_mulhi(vec[i], c_pr) * n;
}

static void
_scalar_addmul_lazy_shoup(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)
{
    const ulong c_pr = _nmod_mulmod_precomp_shoup(c, mod), n = mod.n, n2 = 2 * mod.n;
    slong i;

    for (i = 0; i < len; i++)
    {
        ulong t = c * vec[i] - n_mulhi(vec[i], c_pr) * n + res[i];
        res[i] = (t >= n2) ? t - n2 : t;
    }
}

void _nmod_vec_scalar_mul_nmod_lazy(nn_ptr res, nn_srcptr vec,
                                    slong len, ulong c, nmod_t mod)
{
    FLINT_ASSERT(NMOD_BITS(mod) <= FLINT_BITS - 2);
    if (c <= UWORD(1) || c == mod.n - UWORD(1))
        _nmod_vec_scalar_mul_nmod(res, vec, len, c, mod);
#if NMOD_VEC_SIMD_ZMM || NMOD_VEC_SIMD_YMM || NMOD_VEC_SIMD_NEON
    else if (_SCALAR_SIMD_APPLIES(len, mod, 0))
        _scalar_mul_lazy_simd(res, vec, len, c, mod);
#endif
    else if (len >= NMOD_VEC_LAZY_SHOUP_MIN_LEN && NMOD_BITS(mod) > 32)
        _scalar_mul_lazy_shoup(res, vec, len, c, mod);
    else
        _nmod_vec_scalar_mul_nmod(res, vec, len, c, mod);
}

void _nmod_vec_scalar_addmul_nmod_lazy(nn_ptr res, nn_srcptr vec,
                                       slong len, ulong c, nmod_t mod)
{
    FLINT_ASSERT(NMOD_BITS(mod) <= FLINT_BITS - 2);
    if (c == UWORD(0))
        return;
#if NMOD_VEC_SIMD_ZMM || NMOD_VEC_SIMD_YMM || NMOD_VEC_SIMD_NEON
    else if (_SCALAR_SIMD_APPLIES(len, mod, 1))
        _scalar_addmul_lazy_simd(res, vec, len, c, mod);
#endif
    else if (NMOD_BITS(mod) <= FLINT_BITS / 2)
        _scalar_addmul_halflimb(res, vec, len, c, mod);
    else if (len >= NMOD_VEC_LAZY_SHOUP_MIN_LEN)
        _scalar_addmul_lazy_shoup(res, vec, len, c, mod);
    else
    {
        _nmod_vec_reduce(res, res, len, mod);
        _nmod_vec_scalar_addmul_nmod(res, vec, len, c, mod);
    }
}
