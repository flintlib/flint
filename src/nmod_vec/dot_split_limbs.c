/*
    Copyright (C) 2026 Marie Bonboire
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Dot products for moduli above 2^32 on x86-64 with AVX2 or AVX-512
    (without IFMA, or with IFMA for moduli of 53 to 58 bits), through the
    32 x 32 -> 64 bit products of vpmuludq. Let b be the number of bits of
    n - 1.

    For b <= 61 the entries are split as x = x1 2^s + x0 with s = ceil(b/2)
    (limbs of s and b - s <= 31 bits), and three accumulators receive
    L = sum x0 y0, M = sum x0 y1 + x1 y0, H = sum x1 y1 (schoolbook: 4
    products per term). After every chunk of F iterations, F the largest
    number of terms that a lane below 2^32 can take without overflow, each
    lane is folded: its high 32 bits go to a second accumulator. At the end

      sum x y = (L + 2^32 Lh) + 2^s (M + 2^32 Mh) + 2^(2s) (H + 2^32 Hh),

    with L, Lh, ... standing for the sums of the lanes.

    For b = 63 or 64 the entries are split into 32-bit halves, and each of
    the four products p (below 2^64) goes as X += p (modulo 2^64) and
    Xh += p >> 32; the sum of the low 32-bit halves of the products is then
    X - 2^32 Xh modulo 2^64, exact while a lane takes less than 2^32
    products (the accumulators are emptied every H32_CHUNK iterations).

    The unreduced total is formed on three limbs and reduced at the end:
    the result is exact for any length and both the two-limb and the
    three-limb bands of _nmod_vec_dot_params are covered. The three
    variants (vec1[i] * vec2[i], vec1[i] * vec2[len-1-i] and
    vec1[i] * vec2[i][offset]) differ only in how the entries of vec2 are
    loaded.
*/

#include "nmod.h"
#include "nmod_vec.h"

#if NMOD_VEC_HAVE_DOT_SPLIT_LIMBS

#include <immintrin.h>

#if defined(__AVX512F__)

#define VT __m512i
#define VL 8
#define V_ZERO() _mm512_setzero_si512()
#define V_SET1(x) _mm512_set1_epi64(x)
#define V_LOAD(p) _mm512_loadu_si512((const void *) (p))
#define V_ADD(a, b) _mm512_add_epi64(a, b)
#define V_SUB(a, b) _mm512_sub_epi64(a, b)
#define V_AND(a, b) _mm512_and_si512(a, b)
#define V_MUL(a, b) _mm512_mul_epu32(a, b)
#define V_SRLI(a, c) _mm512_srli_epi64(a, c)
#define V_SLLI(a, c) _mm512_slli_epi64(a, c)
#define V_SRLV(a, c) _mm512_srlv_epi64(a, c)
#define V_HSUM(a) ((ulong) _mm512_reduce_add_epi64(a))

#else  /* AVX2 */

#define VT __m256i
#define VL 4
#define V_ZERO() _mm256_setzero_si256()
#define V_SET1(x) _mm256_set1_epi64x(x)
#define V_LOAD(p) _mm256_loadu_si256((const __m256i *) (p))
#define V_ADD(a, b) _mm256_add_epi64(a, b)
#define V_SUB(a, b) _mm256_sub_epi64(a, b)
#define V_AND(a, b) _mm256_and_si256(a, b)
#define V_MUL(a, b) _mm256_mul_epu32(a, b)
#define V_SRLI(a, c) _mm256_srli_epi64(a, c)
#define V_SLLI(a, c) _mm256_slli_epi64(a, c)
#define V_SRLV(a, c) _mm256_srlv_epi64(a, c)
#define V_HSUM(a) split_limbs_hsum256(a)

FLINT_FORCE_INLINE ulong
split_limbs_hsum256(__m256i a)
{
    __m128i s = _mm_add_epi64(_mm256_castsi256_si128(a),
                              _mm256_extracti128_si256(a, 1));
    return (ulong) _mm_cvtsi128_si64(s) + (ulong) _mm_extract_epi64(s, 1);
}

#endif

/* chunks between two emptyings of the high accumulators (b <= 61): a lane
   of Xh stays below 2^45 */
#define SPLIT_DUMP 8192
/* iterations between two emptyings of the accumulators (b >= 63): 2^14
   products per lane of M, sums of the lanes below 2^49 */
#define H32_CHUNK 8192

/* t2:t1:t0 += x * 2^sh, 0 <= sh < 128 */
FLINT_FORCE_INLINE void
split_limbs_add_shifted(ulong * t2, ulong * t1, ulong * t0, ulong x, int sh)
{
    ulong a2, a1, a0;

    if (sh == 0)
        a2 = 0, a1 = 0, a0 = x;
    else if (sh < 64)
        a2 = 0, a1 = x >> (64 - sh), a0 = x << sh;
    else if (sh == 64)
        a2 = 0, a1 = x, a0 = 0;
    else
        a2 = x >> (128 - sh), a1 = x << (sh - 64), a0 = 0;

    add_sssaaaaaa(*t2, *t1, *t0, *t2, *t1, *t0, a2, a1, a0);
}

#define SPLIT_LIMBS_FOLD(x, xh) \
    do { xh = V_ADD(xh, V_SRLI(x, 32)); x = V_AND(x, m32); } while (0)

#define SPLIT_LIMBS_H32_ADD(x, xh, p) \
    do { x = V_ADD(x, p); xh = V_ADD(xh, V_SRLI(p, 32)); } while (0)

/*
    LOAD2(k): the vector of entries k, ..., k + VL - 1 of vec2 (in the order
    of the variant); ELT2(k): entry k.
*/
#define SPLIT_LIMBS_BODY(LOAD2, ELT2)                                       \
    const int b = FLINT_BIT_COUNT(mod.n - 1);                               \
    const VT m32 = V_SET1(UWORD(0xFFFFFFFF));                               \
    VT L = V_ZERO(), M = V_ZERO(), H = V_ZERO();                            \
    VT Lh = V_ZERO(), Mh = V_ZERO(), Hh = V_ZERO();                         \
    ulong t2 = 0, t1 = 0, t0 = 0, res;                                      \
    slong i = 0, stop;                                                      \
    const slong vecs = len - len % VL;                                      \
                                                                            \
    if (b <= 61)                                                            \
    {                                                                       \
        const int s = (b + 1) / 2;                                          \
        const ulong l0 = (UWORD(1) << s) - 1, l1 = (UWORD(1) << (b - s)) - 1; \
        const ulong pmax = FLINT_MAX(l0 * l0, 2 * l0 * l1);                 \
        const slong F = (slong) FLINT_MIN((UWORD_MAX - UWORD(0xFFFFFFFF)) / pmax, \
                                          UWORD(1) << 40);                  \
        const VT mk = V_SET1(l0), sv = V_SET1(s);                           \
        slong nf = 0;                                                       \
                                                                            \
        while (i < vecs)                                                    \
        {                                                                   \
            stop = FLINT_MIN(vecs, i + F * VL);                             \
            for ( ; i < stop; i += VL)                                      \
            {                                                               \
                VT a = V_LOAD(vec1 + i), c = LOAD2(i);                      \
                VT a0 = V_AND(a, mk), a1 = V_SRLV(a, sv);                   \
                VT c0 = V_AND(c, mk), c1 = V_SRLV(c, sv);                   \
                L = V_ADD(L, V_MUL(a0, c0));                                \
                H = V_ADD(H, V_MUL(a1, c1));                                \
                M = V_ADD(M, V_ADD(V_MUL(a0, c1), V_MUL(a1, c0)));          \
            }                                                               \
            /* lanes below 2^32 at the start of each chunk and at the end */ \
            SPLIT_LIMBS_FOLD(L, Lh);                                        \
            SPLIT_LIMBS_FOLD(M, Mh);                                        \
            SPLIT_LIMBS_FOLD(H, Hh);                                        \
            if (++nf == SPLIT_DUMP)                                         \
            {                                                               \
                split_limbs_add_shifted(&t2, &t1, &t0, V_HSUM(Lh), 32);     \
                split_limbs_add_shifted(&t2, &t1, &t0, V_HSUM(Mh), s + 32); \
                split_limbs_add_shifted(&t2, &t1, &t0, V_HSUM(Hh), 2*s + 32); \
                Lh = Mh = Hh = V_ZERO();                                    \
                nf = 0;                                                     \
            }                                                               \
        }                                                                   \
                                                                            \
        split_limbs_add_shifted(&t2, &t1, &t0, V_HSUM(L), 0);               \
        split_limbs_add_shifted(&t2, &t1, &t0, V_HSUM(Lh), 32);             \
        split_limbs_add_shifted(&t2, &t1, &t0, V_HSUM(M), s);               \
        split_limbs_add_shifted(&t2, &t1, &t0, V_HSUM(Mh), s + 32);         \
        split_limbs_add_shifted(&t2, &t1, &t0, V_HSUM(H), 2 * s);           \
        split_limbs_add_shifted(&t2, &t1, &t0, V_HSUM(Hh), 2 * s + 32);     \
    }                                                                       \
    else                                                                    \
    {                                                                       \
        while (i < vecs)                                                    \
        {                                                                   \
            stop = FLINT_MIN(vecs, i + H32_CHUNK * VL);                     \
            for ( ; i < stop; i += VL)                                      \
            {                                                               \
                VT a = V_LOAD(vec1 + i), c = LOAD2(i);                      \
                VT a1 = V_SRLI(a, 32), c1 = V_SRLI(c, 32);                  \
                SPLIT_LIMBS_H32_ADD(L, Lh, V_MUL(a, c));                    \
                SPLIT_LIMBS_H32_ADD(H, Hh, V_MUL(a1, c1));                  \
                SPLIT_LIMBS_H32_ADD(M, Mh, V_MUL(a, c1));                   \
                SPLIT_LIMBS_H32_ADD(M, Mh, V_MUL(a1, c));                   \
            }                                                               \
            L = V_SUB(L, V_SLLI(Lh, 32));                                   \
            M = V_SUB(M, V_SLLI(Mh, 32));                                   \
            H = V_SUB(H, V_SLLI(Hh, 32));                                   \
            split_limbs_add_shifted(&t2, &t1, &t0, V_HSUM(L), 0);           \
            split_limbs_add_shifted(&t2, &t1, &t0, V_HSUM(Lh), 32);         \
            split_limbs_add_shifted(&t2, &t1, &t0, V_HSUM(M), 32);          \
            split_limbs_add_shifted(&t2, &t1, &t0, V_HSUM(Mh), 64);         \
            split_limbs_add_shifted(&t2, &t1, &t0, V_HSUM(H), 64);          \
            split_limbs_add_shifted(&t2, &t1, &t0, V_HSUM(Hh), 96);         \
            L = M = H = Lh = Mh = Hh = V_ZERO();                            \
        }                                                                   \
    }                                                                       \
                                                                            \
    for ( ; i < len; i++)                                                   \
    {                                                                       \
        ulong p1, p0;                                                       \
        umul_ppmm(p1, p0, vec1[i], ELT2(i));                                \
        add_sssaaaaaa(t2, t1, t0, t2, t1, t0, UWORD(0), p1, p0);            \
    }                                                                       \
                                                                            \
    if (t2 == 0)                                                            \
        NMOD2_RED2(res, t1, t0, mod);                                       \
    else                                                                    \
    {                                                                       \
        NMOD_RED(t2, t2, mod);                                              \
        NMOD_RED3(res, t2, t1, t0, mod);                                    \
    }                                                                       \
    return res;

ulong
_nmod_vec_dot_split_limbs(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
{
    /* the parameters may come from a longer length: short products are
       faster with the scalar code */
    if (len < NMOD_VEC_DOT_SPLIT_LIMBS_MIN_LEN)
        return _nmod_vec_dot(vec1, vec2, len, mod, _nmod_vec_dot_params(len, mod));

#define LOAD2(k) V_LOAD(vec2 + (k))
#define ELT2(k) vec2[k]
    SPLIT_LIMBS_BODY(LOAD2, ELT2)
#undef LOAD2
#undef ELT2
}

ulong
_nmod_vec_dot_split_limbs_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
{
    if (len < NMOD_VEC_DOT_SPLIT_LIMBS_MIN_LEN)
        return _nmod_vec_dot_rev(vec1, vec2, len, mod, _nmod_vec_dot_params(len, mod));

#if defined(__AVX512F__)
    const __m512i rev = _mm512_set_epi64(0, 1, 2, 3, 4, 5, 6, 7);
#define LOAD2(k) _mm512_permutexvar_epi64(rev, V_LOAD(vec2 + len - 8 - (k)))
#else
#define LOAD2(k) _mm256_permute4x64_epi64(V_LOAD(vec2 + len - 4 - (k)), 0x1B)
#endif
#define ELT2(k) vec2[len - 1 - (k)]
    SPLIT_LIMBS_BODY(LOAD2, ELT2)
#undef LOAD2
#undef ELT2
}

ulong
_nmod_vec_dot_split_limbs_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset,
                              slong len, nmod_t mod)
{
    if (len < NMOD_VEC_DOT_SPLIT_LIMBS_MIN_LEN)
        return _nmod_vec_dot_ptr(vec1, vec2, offset, len, mod, _nmod_vec_dot_params(len, mod));

#if defined(__AVX512F__)
    const void * base = (const void *) (offset * (slong) sizeof(ulong));
#define LOAD2(k) _mm512_i64gather_epi64(V_LOAD(vec2 + (k)), base, 1)
#else
    /* gathers are slow on some AVX2 machines */
#define LOAD2(k) _mm256_set_epi64x(vec2[(k) + 3][offset], vec2[(k) + 2][offset], \
                                   vec2[(k) + 1][offset], vec2[k][offset])
#endif
#define ELT2(k) vec2[k][offset]
    SPLIT_LIMBS_BODY(LOAD2, ELT2)
#undef LOAD2
#undef ELT2
}

#else

/* never selected by _nmod_vec_dot_params without NMOD_VEC_HAVE_DOT_SPLIT_LIMBS */
ulong
_nmod_vec_dot_split_limbs(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
{
    return _nmod_vec_dot3(vec1, vec2, len, mod);
}

ulong
_nmod_vec_dot_split_limbs_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
{
    return _nmod_vec_dot3_rev(vec1, vec2, len, mod);
}

ulong
_nmod_vec_dot_split_limbs_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset,
                              slong len, nmod_t mod)
{
    return _nmod_vec_dot3_ptr(vec1, vec2, offset, len, mod);
}

#endif
