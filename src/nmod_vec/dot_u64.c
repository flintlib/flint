/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Dot products for moduli above 2^52 (up to 2^64) with AVX512-IFMA. The
    entries are split into 32-bit halves, x = x1 2^32 + x0, and the four
    products x0 y0, x0 y1, x1 y0, x1 y1 (below 2^64) are accumulated as
    their low and high 52-bit halves (vpmadd52luq / huq): eight
    accumulators, whose independence hides the latency of the fused
    multiply-adds, eight IFMA, four shifts / masks and two loads per 8
    terms. Every U64_DOT_CHUNK iterations the lanes (below 2^60) are
    summed and the accumulators combined into a three-limb total,

      x y = x0 y0 + 2^32 (x0 y1 + x1 y0) + 2^64 x1 y1,

    reduced at the very end; the result is exact for any length. The tail
    of fewer than 8 terms goes through masked loads.

    The three variants (vec1[i] * vec2[i], vec1[i] * vec2[len-1-i] and
    vec1[i] * vec2[i][offset]) differ only in how the 8 entries of vec2
    are gathered.
*/

#include "nmod.h"
#include "nmod_vec.h"

#if NMOD_VEC_HAVE_DOT_U64

#include <immintrin.h>

/* iterations between two dumps: 250 + 1 (masked tail) halves per lane,
   so that the sums of two accumulators stay below 2^61 */
#define U64_DOT_CHUNK 250

typedef struct
{
    __m512i l00, h00, l01, h01, l10, h10, l11, h11;
} u64_dot_acc;

FLINT_FORCE_INLINE void
u64_dot_acc_zero(u64_dot_acc * a)
{
    a->l00 = a->h00 = a->l01 = a->h01 = _mm512_setzero_si512();
    a->l10 = a->h10 = a->l11 = a->h11 = _mm512_setzero_si512();
}

FLINT_FORCE_INLINE void
u64_dot_acc_step(u64_dot_acc * a, __m512i x, __m512i y)
{
    const __m512i m32 = _mm512_set1_epi64(UWORD(0xFFFFFFFF));
    const __m512i x0 = _mm512_and_si512(x, m32), x1 = _mm512_srli_epi64(x, 32);
    const __m512i y0 = _mm512_and_si512(y, m32), y1 = _mm512_srli_epi64(y, 32);

    a->l00 = _mm512_madd52lo_epu64(a->l00, x0, y0);
    a->h00 = _mm512_madd52hi_epu64(a->h00, x0, y0);
    a->l01 = _mm512_madd52lo_epu64(a->l01, x0, y1);
    a->h01 = _mm512_madd52hi_epu64(a->h01, x0, y1);
    a->l10 = _mm512_madd52lo_epu64(a->l10, x1, y0);
    a->h10 = _mm512_madd52hi_epu64(a->h10, x1, y0);
    a->l11 = _mm512_madd52lo_epu64(a->l11, x1, y1);
    a->h11 = _mm512_madd52hi_epu64(a->h11, x1, y1);
}

/* t2:t1:t0 += s * 2^shift, s < 2^64, 0 < shift < 64 */
#define U64_DOT_ADD_SHIFTED(t2, t1, t0, s, shift)                       \
    add_sssaaaaaa(t2, t1, t0, t2, t1, t0,                               \
                  UWORD(0), (s) >> (64 - (shift)), (s) << (shift))
#define U64_DOT_ADD_SHIFTED64(t2, t1, t0, s, shift)                     \
    add_sssaaaaaa(t2, t1, t0, t2, t1, t0,                               \
                  (s) >> (64 - (shift)), (s) << (shift), UWORD(0))

/* t2:t1:t0 += the total of the accumulators, whose lanes are below 2^60 */
FLINT_FORCE_INLINE void
u64_dot_acc_dump(ulong * t2, ulong * t1, ulong * t0, const u64_dot_acc * a)
{
    ulong s;

    s = _mm512_reduce_add_epi64(a->l00);
    add_sssaaaaaa(*t2, *t1, *t0, *t2, *t1, *t0, UWORD(0), UWORD(0), s);
    s = _mm512_reduce_add_epi64(a->h00);
    U64_DOT_ADD_SHIFTED(*t2, *t1, *t0, s, 52);
    s = _mm512_reduce_add_epi64(_mm512_add_epi64(a->l01, a->l10));
    U64_DOT_ADD_SHIFTED(*t2, *t1, *t0, s, 32);
    s = _mm512_reduce_add_epi64(_mm512_add_epi64(a->h01, a->h10));
    U64_DOT_ADD_SHIFTED64(*t2, *t1, *t0, s, 20);
    s = _mm512_reduce_add_epi64(a->l11);
    add_sssaaaaaa(*t2, *t1, *t0, *t2, *t1, *t0, UWORD(0), s, UWORD(0));
    s = _mm512_reduce_add_epi64(a->h11);
    U64_DOT_ADD_SHIFTED64(*t2, *t1, *t0, s, 52);
}

#define U64_DOT_BODY(LOAD2, LOAD2_MASKED)                                   \
    u64_dot_acc acc;                                                        \
    ulong t2 = 0, t1 = 0, t0 = 0, res;                                      \
    slong i = 0, stop;                                                      \
    const slong vecs = len - len % 8;                                       \
                                                                            \
    u64_dot_acc_zero(&acc);                                                 \
                                                                            \
    while (i < vecs)                                                        \
    {                                                                       \
        stop = FLINT_MIN(vecs, i + 8 * U64_DOT_CHUNK);                      \
        for ( ; i < stop; i += 8)                                           \
            u64_dot_acc_step(&acc,                                          \
                _mm512_loadu_si512((const void *) (vec1 + i)), LOAD2(i));   \
        if (i < vecs)                                                       \
        {                                                                   \
            u64_dot_acc_dump(&t2, &t1, &t0, &acc);                          \
            u64_dot_acc_zero(&acc);                                         \
        }                                                                   \
    }                                                                       \
                                                                            \
    if (i < len)                                                            \
    {                                                                       \
        const __mmask8 m = (__mmask8) ((1 << (len - i)) - 1);               \
        u64_dot_acc_step(&acc,                                              \
                _mm512_maskz_loadu_epi64(m, (const void *) (vec1 + i)),     \
                LOAD2_MASKED(i, m));                                        \
    }                                                                       \
                                                                            \
    u64_dot_acc_dump(&t2, &t1, &t0, &acc);                                  \
    NMOD_RED(t2, t2, mod);                                                  \
    NMOD_RED3(res, t2, t1, t0, mod);                                        \
    return res;

ulong
_nmod_vec_dot_u64(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
{
    /* the parameters may come from a longer length: short products are
       faster with the scalar code */
    if (len < NMOD_VEC_DOT_U64_MIN_LEN)
        return _nmod_vec_dot3(vec1, vec2, len, mod);

#define LOAD2(k) _mm512_loadu_si512((const void *) (vec2 + (k)))
#define LOAD2_MASKED(k, m) _mm512_maskz_loadu_epi64(m, (const void *) (vec2 + (k)))
    U64_DOT_BODY(LOAD2, LOAD2_MASKED)
#undef LOAD2
#undef LOAD2_MASKED
}

ulong
_nmod_vec_dot_u64_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
{
    if (len < NMOD_VEC_DOT_U64_MIN_LEN)
        return _nmod_vec_dot3_rev(vec1, vec2, len, mod);

    /* as in dot_u52.c */
    const __m512i rev = _mm512_set_epi64(0, 1, 2, 3, 4, 5, 6, 7);
#define LOAD2(k) _mm512_permutexvar_epi64(rev, \
        _mm512_loadu_si512((const void *) (vec2 + len - 8 - (k))))
#define LOAD2_MASKED(k, m) _mm512_permutexvar_epi64( \
        _mm512_sub_epi64(_mm512_set1_epi64(len - (k) - 1), \
                         _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0)), \
        _mm512_maskz_loadu_epi64(m, (const void *) vec2))
    U64_DOT_BODY(LOAD2, LOAD2_MASKED)
#undef LOAD2
#undef LOAD2_MASKED
}

ulong
_nmod_vec_dot_u64_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset,
                      slong len, nmod_t mod)
{
    if (len < NMOD_VEC_DOT_U64_MIN_LEN)
        return _nmod_vec_dot3_ptr(vec1, vec2, offset, len, mod);

    const void * base = (const void *) (offset * (slong) sizeof(ulong));
#define LOAD2(k) _mm512_i64gather_epi64( \
        _mm512_loadu_si512((const void *) (vec2 + (k))), base, 1)
#define LOAD2_MASKED(k, m) _mm512_mask_i64gather_epi64( \
        _mm512_setzero_si512(), m, \
        _mm512_maskz_loadu_epi64(m, (const void *) (vec2 + (k))), base, 1)
    U64_DOT_BODY(LOAD2, LOAD2_MASKED)
#undef LOAD2
#undef LOAD2_MASKED
}

#else

/* never selected by _nmod_vec_dot_params without AVX512-IFMA */
ulong
_nmod_vec_dot_u64(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
{
    return _nmod_vec_dot3(vec1, vec2, len, mod);
}

ulong
_nmod_vec_dot_u64_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
{
    return _nmod_vec_dot3_rev(vec1, vec2, len, mod);
}

ulong
_nmod_vec_dot_u64_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset,
                      slong len, nmod_t mod)
{
    return _nmod_vec_dot3_ptr(vec1, vec2, offset, len, mod);
}

#endif
