/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Double precision vector primitives for modular arithmetic, in the style
    of fft_small (machine_vectors.h): a vector fpv of FPV_VL doubles and the
    operations the matrix multiplication kernels of mul_fp50.c (whole kernel
    in floating point) and mul_k52.c (integer kernel, floating point
    reduction at the end of a block) need.

    Backends: AVX-512 (F + DQ, 8 lanes), AVX2 + FMA (4 lanes), AArch64 NEON
    (2 lanes), plain C (a struct of 4 doubles the compiler may vectorize).
    The selection follows machine_vectors.h (FLINT_MACHINE_VECTORS_GENERIC
    forces the plain C path); NMOD_MAT_FP_FORCE_AVX2 / _FORCE_GENERIC select
    a lesser backend on purpose, for testing and profiling. The includer may
    test FPV_AVX512 / FPV_AVX2 / FPV_NEON / FPV_GENERIC to define matching
    integer primitives (mul_k52.c does), and gets fpv_i, the integer vector
    type of the same width.

    The primitives are not taken from machine_vectors.h because the kernels
    want a fixed lane count per ISA (vec8d is only an emulation on AVX2),
    conversions between 64-bit integers and doubles that it does not
    provide, and a NEON tier. Folding these back into machine_vectors.h is
    a possible future perspective.

    Rounding to the nearest integer. Every rounding here is of a product,
    through fpv_rint_mul(x, y) = rint(x*y), which requires |x*y| <= 2^51
    and is a fused multiply-add of the constant 3*2^51 followed by a
    subtraction of it: the sum lands in the binade [2^52, 2^53) where the
    unit in the last place is 1, hence rounds to the nearest integer (ties
    to even), with a quotient error of at most 1/2 rather than the
    1/2 + 2^-(B+1) of a separately rounded product (this is what the
    analysis in fft_small/mulmod_satisfies_bounds.c calls the optimal
    q_error). The rounding instructions (vroundpd, vrndscalepd) have no
    such restriction but cost one operation more here and are slower on
    every Intel core measured: a standalone microkernel probe of the fp50
    kernel gains 50% on Arrow Lake (AVX2), 10-12% on Ice Lake and 16-22%
    on Cascade Lake (AVX-512) from dropping them, and is unchanged on Zen
    4, Zen 5 and Apple M4. NEON keeps frintn, which is one cheap
    instruction, has no restriction, and leaves the operation count equal.

    All the uses here stay inside |x*y| <= 2^51: mul_fp50.c rounds a*b/n
    with |a*b| <= n^2/4 and n < 2^50, and acc/n with |acc| <= (1 + 9/8 KC) n;
    mul_k52.c rounds hi*c/n with |hi| < 2^31, and r*c/n with |r| <= n/2 and
    c < n <= 2^52, whence at most (n-1)/2 < 2^51.

    Modular reduction (the bounds are those of fft_small, see
    fft_small/mulmod_satisfies_bounds.c for the analysis, with prec = 53):

      fpv_mulmod(a, b, n, ninv)     a*b - q*n for an integer q near a*b/n:
                                    for n < 2^50 and |a*b| < 2 n^2 the result
                                    lies in (-9/8 n, 9/8 n) (fft_small); for
                                    any n <= 2^52, |a*b| <= 2^31 n gives
                                    (-0.51 n, 0.51 n) and |a*b| <= 0.51 n^2
                                    gives (-1.26 n, 1.26 n). (The quotient is
                                    off by at most |a*b| 2^-53 / n from the
                                    error of ninv, 2^(e-54) / n from the low
                                    part when |a*b| < 2^e, and 1/2 + 2^-(B+1)
                                    from the rounding when |h*ninv| < 2^(52-B).)
      fpv_reduce_pm1n(x, n, ninv)   x - n*rint(x/n) up to a quotient error
                                    below 2^-12, so in (-0.51 n, 0.51 n),
                                    for integer |x| < min(2^53, 2^40 n)
      fpv_rint_mul(x, y)            rint(x*y), for |x*y| <= 2^51
      fpv_pm1n_to_pmhn(x, n)        (-3/2 n, 3/2 n) -> [-n/2, n/2]
      fpv_reduce_0n(x, n)           [-n, n] -> [0, n)

    all of them for integer-valued doubles; ninv is the double nearest 1/n.
*/

#ifndef NMOD_MAT_MUL_FP_VEC_H
#define NMOD_MAT_MUL_FP_VEC_H

#include "flint.h"
#include "machine_vectors.h"

/* backend selection *********************************************************/

#if defined(FLINT_MACHINE_VECTORS_GENERIC) \
        || defined(NMOD_MAT_FP_FORCE_GENERIC)
# define FPV_GENERIC 1
#elif defined(__AVX512F__) && defined(__AVX512DQ__) \
        && !defined(NMOD_MAT_FP_FORCE_AVX2)
# define FPV_AVX512 1
#elif defined(__AVX2__) && defined(__FMA__)
# define FPV_AVX2 1
#elif (defined(__ARM_NEON) && defined(__aarch64__)) || defined(_M_ARM64)
# define FPV_NEON 1
#else
# define FPV_GENERIC 1
#endif

/* the rounding constant of fpv_rint_mul, 3*2^51 */
#define FPV_RND_MAGIC 0x1.8p52

#if defined(FPV_AVX512)

/* AVX-512 (F + DQ) *********************************************************/

#define FPV_VL 8
typedef __m512d fpv;
typedef __m512i fpv_i;

FLINT_FORCE_INLINE fpv fpv_zero(void) { return _mm512_setzero_pd(); }
FLINT_FORCE_INLINE fpv fpv_set1(double a) { return _mm512_set1_pd(a); }
FLINT_FORCE_INLINE fpv fpv_load(const double * p) { return _mm512_load_pd(p); }
FLINT_FORCE_INLINE fpv fpv_loadu(const double * p) { return _mm512_loadu_pd(p); }
FLINT_FORCE_INLINE void fpv_storeu(double * p, fpv a) { _mm512_storeu_pd(p, a); }
FLINT_FORCE_INLINE fpv fpv_add(fpv a, fpv b) { return _mm512_add_pd(a, b); }
FLINT_FORCE_INLINE fpv fpv_sub(fpv a, fpv b) { return _mm512_sub_pd(a, b); }
FLINT_FORCE_INLINE fpv fpv_mul(fpv a, fpv b) { return _mm512_mul_pd(a, b); }
/* a*b + c and -a*b + c */
FLINT_FORCE_INLINE fpv fpv_fmadd(fpv a, fpv b, fpv c) { return _mm512_fmadd_pd(a, b, c); }
FLINT_FORCE_INLINE fpv fpv_fnmadd(fpv a, fpv b, fpv c) { return _mm512_fnmadd_pd(a, b, c); }

/* rint(a*b) for |a*b| <= 2^51; see the header comment */
FLINT_FORCE_INLINE fpv
fpv_rint_mul(fpv a, fpv b)
{
    const fpv magic = _mm512_set1_pd(FPV_RND_MAGIC);
    return _mm512_sub_pd(_mm512_fmadd_pd(a, b, magic), magic);
}

/* x < 0 ? x + n : x, then x >= n ? x - n : x */
FLINT_FORCE_INLINE fpv
fpv_reduce_0n(fpv x, fpv n)
{
    __mmask8 m;

    m = _mm512_cmp_pd_mask(x, _mm512_setzero_pd(), _CMP_LT_OQ);
    x = _mm512_mask_add_pd(x, m, x, n);
    m = _mm512_cmp_pd_mask(x, n, _CMP_GE_OQ);
    x = _mm512_mask_sub_pd(x, m, x, n);

    return x;
}

/* x > n/2 ? x - n : (x < -n/2 ? x + n : x) */
FLINT_FORCE_INLINE fpv
fpv_pm1n_to_pmhn(fpv x, fpv n)
{
    fpv halfn = _mm512_mul_pd(n, _mm512_set1_pd(0.5));
    __mmask8 m;

    m = _mm512_cmp_pd_mask(x, halfn, _CMP_GT_OQ);
    x = _mm512_mask_sub_pd(x, m, x, n);
    m = _mm512_cmp_pd_mask(x, _mm512_sub_pd(_mm512_setzero_pd(), halfn), _CMP_LT_OQ);
    x = _mm512_mask_add_pd(x, m, x, n);

    return x;
}

/* unsigned 64-bit lanes below 2^53 <-> doubles */
FLINT_FORCE_INLINE fpv
fpv_load_u64(const ulong * p)
{
    return _mm512_cvtepu64_pd(_mm512_loadu_si512((const void *) p));
}

FLINT_FORCE_INLINE void
fpv_store_u64(ulong * p, fpv a)
{
    _mm512_storeu_si512((void *) p, _mm512_cvtpd_epu64(a));
}

/* signed 64-bit lanes x = hi*2^32 + lo, hi signed, 0 <= lo < 2^32, as two
   exact doubles */
FLINT_FORCE_INLINE void
fpv_split_i64(fpv_i x, fpv * hd, fpv * ld)
{
    __m512i hi = _mm512_srai_epi64(x, 32);
    __m512i lo = _mm512_and_si512(x, _mm512_set1_epi64(UWORD(0xFFFFFFFF)));

    *hd = _mm512_cvtepi64_pd(hi);
    *ld = _mm512_cvtepi64_pd(lo);
}

#elif defined(FPV_AVX2)

/* AVX2 + FMA ***************************************************************/

#define FPV_VL 4
typedef __m256d fpv;
typedef __m256i fpv_i;

FLINT_FORCE_INLINE fpv fpv_zero(void) { return _mm256_setzero_pd(); }
FLINT_FORCE_INLINE fpv fpv_set1(double a) { return _mm256_set1_pd(a); }
FLINT_FORCE_INLINE fpv fpv_load(const double * p) { return _mm256_load_pd(p); }
FLINT_FORCE_INLINE fpv fpv_loadu(const double * p) { return _mm256_loadu_pd(p); }
FLINT_FORCE_INLINE void fpv_storeu(double * p, fpv a) { _mm256_storeu_pd(p, a); }
FLINT_FORCE_INLINE fpv fpv_add(fpv a, fpv b) { return _mm256_add_pd(a, b); }
FLINT_FORCE_INLINE fpv fpv_sub(fpv a, fpv b) { return _mm256_sub_pd(a, b); }
FLINT_FORCE_INLINE fpv fpv_mul(fpv a, fpv b) { return _mm256_mul_pd(a, b); }
FLINT_FORCE_INLINE fpv fpv_fmadd(fpv a, fpv b, fpv c) { return _mm256_fmadd_pd(a, b, c); }
FLINT_FORCE_INLINE fpv fpv_fnmadd(fpv a, fpv b, fpv c) { return _mm256_fnmadd_pd(a, b, c); }

/* rint(a*b) for |a*b| <= 2^51 */
FLINT_FORCE_INLINE fpv
fpv_rint_mul(fpv a, fpv b)
{
    const fpv magic = _mm256_set1_pd(FPV_RND_MAGIC);
    return _mm256_sub_pd(_mm256_fmadd_pd(a, b, magic), magic);
}

FLINT_FORCE_INLINE fpv
fpv_reduce_0n(fpv x, fpv n)
{
    fpv m;

    m = _mm256_cmp_pd(x, _mm256_setzero_pd(), _CMP_LT_OQ);
    x = _mm256_add_pd(x, _mm256_and_pd(m, n));
    m = _mm256_cmp_pd(x, n, _CMP_GE_OQ);
    x = _mm256_sub_pd(x, _mm256_and_pd(m, n));

    return x;
}

FLINT_FORCE_INLINE fpv
fpv_pm1n_to_pmhn(fpv x, fpv n)
{
    fpv halfn = _mm256_mul_pd(n, _mm256_set1_pd(0.5));
    fpv m;

    m = _mm256_cmp_pd(x, halfn, _CMP_GT_OQ);
    x = _mm256_sub_pd(x, _mm256_and_pd(m, n));
    m = _mm256_cmp_pd(x, _mm256_sub_pd(_mm256_setzero_pd(), halfn), _CMP_LT_OQ);
    x = _mm256_add_pd(x, _mm256_and_pd(m, n));

    return x;
}

/*
    AVX2 has no conversion between 64-bit integers and doubles. For
    0 <= u < 2^52 the bit pattern 2^52 | u is the double 2^52 + u, so u
    is recovered by one subtraction; the other direction adds 2^52 to a
    double in [0, 2^52) and drops the exponent bits.
*/
#define FPV_MAGIC52 UWORD(0x4330000000000000)

FLINT_FORCE_INLINE fpv
fpv_load_u64(const ulong * p)
{
    const __m256i magic = _mm256_set1_epi64x(FPV_MAGIC52);
    __m256i u = _mm256_loadu_si256((const __m256i *) p);

    return _mm256_sub_pd(_mm256_castsi256_pd(_mm256_or_si256(u, magic)),
                         _mm256_castsi256_pd(magic));
}

FLINT_FORCE_INLINE void
fpv_store_u64(ulong * p, fpv a)
{
    const __m256i magic = _mm256_set1_epi64x(FPV_MAGIC52);
    __m256d t = _mm256_add_pd(a, _mm256_castsi256_pd(magic));

    _mm256_storeu_si256((__m256i *) p,
                        _mm256_xor_si256(_mm256_castpd_si256(t), magic));
}

/* the high word is signed: offset it by 2^31 to make it a 32-bit unsigned
   value, convert, and take the offset back */
FLINT_FORCE_INLINE void
fpv_split_i64(fpv_i x, fpv * hd, fpv * ld)
{
    const __m256i magic = _mm256_set1_epi64x(FPV_MAGIC52);
    const __m256i sign32 = _mm256_set1_epi64x(UWORD(0x80000000));
    __m256i hi = _mm256_srli_epi64(x, 32);
    __m256i lo = _mm256_and_si256(x, _mm256_set1_epi64x(UWORD(0xFFFFFFFF)));

    *ld = _mm256_sub_pd(_mm256_castsi256_pd(_mm256_or_si256(lo, magic)),
                        _mm256_castsi256_pd(magic));
    *hd = _mm256_sub_pd(_mm256_castsi256_pd(_mm256_or_si256(
                            _mm256_xor_si256(hi, sign32), magic)),
                        _mm256_set1_pd(0x1.0p52 + 0x1.0p31));
}

#elif defined(FPV_NEON)

/* AArch64 NEON *************************************************************/

#define FPV_VL 2
typedef float64x2_t fpv;
typedef int64x2_t fpv_i;

FLINT_FORCE_INLINE fpv fpv_zero(void) { return vdupq_n_f64(0.0); }
FLINT_FORCE_INLINE fpv fpv_set1(double a) { return vdupq_n_f64(a); }
FLINT_FORCE_INLINE fpv fpv_load(const double * p) { return vld1q_f64(p); }
FLINT_FORCE_INLINE fpv fpv_loadu(const double * p) { return vld1q_f64(p); }
FLINT_FORCE_INLINE void fpv_storeu(double * p, fpv a) { vst1q_f64(p, a); }
FLINT_FORCE_INLINE fpv fpv_add(fpv a, fpv b) { return vaddq_f64(a, b); }
FLINT_FORCE_INLINE fpv fpv_sub(fpv a, fpv b) { return vsubq_f64(a, b); }
FLINT_FORCE_INLINE fpv fpv_mul(fpv a, fpv b) { return vmulq_f64(a, b); }
/* vfmaq_f64(c, a, b) = c + a*b */
FLINT_FORCE_INLINE fpv fpv_fmadd(fpv a, fpv b, fpv c) { return vfmaq_f64(c, a, b); }
FLINT_FORCE_INLINE fpv fpv_fnmadd(fpv a, fpv b, fpv c) { return vfmsq_f64(c, a, b); }

/* frintn is one instruction and unrestricted, so no rounding constant here */
FLINT_FORCE_INLINE fpv
fpv_rint_mul(fpv a, fpv b) { return vrndnq_f64(vmulq_f64(a, b)); }

FLINT_FORCE_INLINE fpv
fpv_reduce_0n(fpv x, fpv n)
{
    uint64x2_t m;

    m = vcltzq_f64(x);
    x = vaddq_f64(x, vreinterpretq_f64_u64(vandq_u64(m, vreinterpretq_u64_f64(n))));
    m = vcgeq_f64(x, n);
    x = vsubq_f64(x, vreinterpretq_f64_u64(vandq_u64(m, vreinterpretq_u64_f64(n))));

    return x;
}

FLINT_FORCE_INLINE fpv
fpv_pm1n_to_pmhn(fpv x, fpv n)
{
    fpv halfn = vmulq_n_f64(n, 0.5);
    uint64x2_t m;

    m = vcgtq_f64(x, halfn);
    x = vsubq_f64(x, vreinterpretq_f64_u64(vandq_u64(m, vreinterpretq_u64_f64(n))));
    m = vcltq_f64(x, vnegq_f64(halfn));
    x = vaddq_f64(x, vreinterpretq_f64_u64(vandq_u64(m, vreinterpretq_u64_f64(n))));

    return x;
}

FLINT_FORCE_INLINE fpv
fpv_load_u64(const ulong * p)
{
    return vcvtq_f64_u64(vld1q_u64((const uint64_t *) p));
}

FLINT_FORCE_INLINE void
fpv_store_u64(ulong * p, fpv a)
{
    vst1q_u64((uint64_t *) p, vcvtq_u64_f64(a));
}

FLINT_FORCE_INLINE void
fpv_split_i64(fpv_i x, fpv * hd, fpv * ld)
{
    int64x2_t hi = vshrq_n_s64(x, 32);
    int64x2_t lo = vandq_s64(x, vdupq_n_s64((int64_t) UWORD(0xFFFFFFFF)));

    *hd = vcvtq_f64_s64(hi);
    *ld = vcvtq_f64_s64(lo);
}

#else

/* plain C ******************************************************************/

#include <math.h>

#define FPV_VL 4
typedef struct { double v[FPV_VL]; } fpv;
typedef struct { slong v[FPV_VL]; } fpv_i;

FLINT_FORCE_INLINE fpv
fpv_set1(double a)
{
    fpv r;
    slong i;
    for (i = 0; i < FPV_VL; i++)
        r.v[i] = a;
    return r;
}

FLINT_FORCE_INLINE fpv fpv_zero(void) { return fpv_set1(0.0); }

FLINT_FORCE_INLINE fpv
fpv_loadu(const double * p)
{
    fpv r;
    slong i;
    for (i = 0; i < FPV_VL; i++)
        r.v[i] = p[i];
    return r;
}

FLINT_FORCE_INLINE fpv fpv_load(const double * p) { return fpv_loadu(p); }

FLINT_FORCE_INLINE void
fpv_storeu(double * p, fpv a)
{
    slong i;
    for (i = 0; i < FPV_VL; i++)
        p[i] = a.v[i];
}

#define FPV_LANEWISE(name, expr)                                            \
FLINT_FORCE_INLINE fpv name(fpv a, fpv b)                                   \
{                                                                           \
    fpv r;                                                                  \
    slong i;                                                                \
    for (i = 0; i < FPV_VL; i++)                                            \
        r.v[i] = expr;                                                      \
    return r;                                                               \
}
FPV_LANEWISE(fpv_add, a.v[i] + b.v[i])
FPV_LANEWISE(fpv_sub, a.v[i] - b.v[i])
FPV_LANEWISE(fpv_mul, a.v[i] * b.v[i])
#undef FPV_LANEWISE

/* the exactness of fpv_mulmod's low part relies on real fused operations:
   the plain C tier uses fma() from the C library (which the compiler turns
   into an instruction on any target that has one) */
#define FPV_FMA3(name, expr)                                                \
FLINT_FORCE_INLINE fpv name(fpv a, fpv b, fpv c)                            \
{                                                                           \
    fpv r;                                                                  \
    slong i;                                                                \
    for (i = 0; i < FPV_VL; i++)                                            \
        r.v[i] = expr;                                                      \
    return r;                                                               \
}
FPV_FMA3(fpv_fmadd, fma(a.v[i], b.v[i], c.v[i]))
FPV_FMA3(fpv_fnmadd, fma(-a.v[i], b.v[i], c.v[i]))
#undef FPV_FMA3

FLINT_FORCE_INLINE fpv
fpv_rint_mul(fpv a, fpv b)
{
    fpv r;
    slong i;
    for (i = 0; i < FPV_VL; i++)
        r.v[i] = fma(a.v[i], b.v[i], FPV_RND_MAGIC) - FPV_RND_MAGIC;
    return r;
}

FLINT_FORCE_INLINE fpv
fpv_reduce_0n(fpv x, fpv n)
{
    slong i;
    for (i = 0; i < FPV_VL; i++)
    {
        if (x.v[i] < 0.0) x.v[i] += n.v[i];
        if (x.v[i] >= n.v[i]) x.v[i] -= n.v[i];
    }
    return x;
}

FLINT_FORCE_INLINE fpv
fpv_pm1n_to_pmhn(fpv x, fpv n)
{
    slong i;
    for (i = 0; i < FPV_VL; i++)
    {
        double halfn = 0.5 * n.v[i];
        if (x.v[i] > halfn) x.v[i] -= n.v[i];
        if (x.v[i] < -halfn) x.v[i] += n.v[i];
    }
    return x;
}

FLINT_FORCE_INLINE fpv
fpv_load_u64(const ulong * p)
{
    fpv r;
    slong i;
    for (i = 0; i < FPV_VL; i++)
        r.v[i] = (double) p[i];
    return r;
}

FLINT_FORCE_INLINE void
fpv_store_u64(ulong * p, fpv a)
{
    slong i;
    for (i = 0; i < FPV_VL; i++)
        p[i] = (ulong) a.v[i];
}

FLINT_FORCE_INLINE void
fpv_split_i64(fpv_i x, fpv * hd, fpv * ld)
{
    slong i;
    for (i = 0; i < FPV_VL; i++)
    {
        hd->v[i] = (double) (x.v[i] >> 32);                   /* arithmetic */
        ld->v[i] = (double) (ulong) (x.v[i] & (slong) UWORD(0xFFFFFFFF));
    }
}

#endif

/* shared across backends ****************************************************/

/*
    The fft_small formula, with h = rounded a*b and q = rint(h/n):
    a*b - q*n = (h - q*n) - (h - a*b), where both parentheses are exact
    fused operations (the second is the rounding error of h). Written
    with two fnmadd rather than an fnmadd and an fmsub so that NEON,
    which has no fused multiply-subtract of that shape, needs no extra
    negation.
*/
FLINT_FORCE_INLINE fpv
fpv_mulmod(fpv a, fpv b, fpv n, fpv ninv)
{
    fpv h = fpv_mul(a, b);
    fpv q = fpv_rint_mul(h, ninv);
    return fpv_sub(fpv_fnmadd(q, n, h), fpv_fnmadd(a, b, h));
}

FLINT_FORCE_INLINE fpv
fpv_reduce_pm1n(fpv x, fpv n, fpv ninv)
{
    return fpv_fnmadd(fpv_rint_mul(x, ninv), n, x);
}

#endif
