/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Dot products for moduli up to 2^52 with AVX512-IFMA, the strategy of
    nmod_mat_mul_u52: the 104-bit product of two entries below 2^52 is
    accumulated as its low and high 52-bit halves (vpmadd52luq / huq) in
    two 64-bit lanes.

    Four independent pairs of accumulators hide the latency of the fused
    multiply-adds: 32 terms per loop iteration, two IFMA per 8 terms, and
    the loop is bound by the two loads per 8 terms. Every U52_DOT_CHUNK
    iterations the four pairs are added (each lane then holds at most
    4 * U52_DOT_CHUNK + 4 halves, below 2^61) and summed across the lanes
    into a two-limb total, reduced at the very end: the total is at most
    len * 2^104, which fits two limbs for any practical length (len below
    2^24), and _nmod_vec_dot_params selects this method only when the
    unreduced dot product fits two limbs. The tail of fewer than 8 terms
    goes through masked loads.

    The three variants (vec1[i] * vec2[i], vec1[i] * vec2[len-1-i] and
    vec1[i] * vec2[i][offset]) differ only in how the 8 entries of vec2
    are gathered.
*/

#include "nmod.h"
#include "nmod_vec.h"

#if NMOD_VEC_HAVE_DOT_U52

#include <immintrin.h>

#define U52_DOT_NACC 4
#define U52_DOT_UNROLL (8 * U52_DOT_NACC)
/* 4 * 126 + 4 vectors of the remainder = 508 halves per lane < 2^61 */
#define U52_DOT_CHUNK 126

typedef struct
{
    __m512i lo[U52_DOT_NACC];
    __m512i hi[U52_DOT_NACC];
} u52_dot_acc;

FLINT_FORCE_INLINE void
u52_dot_acc_zero(u52_dot_acc * a)
{
    int j;
    for (j = 0; j < U52_DOT_NACC; j++)
    {
        a->lo[j] = _mm512_setzero_si512();
        a->hi[j] = _mm512_setzero_si512();
    }
}

FLINT_FORCE_INLINE void
u52_dot_acc_step(u52_dot_acc * a, int j, __m512i x, __m512i y)
{
    a->lo[j] = _mm512_madd52lo_epu64(a->lo[j], x, y);
    a->hi[j] = _mm512_madd52hi_epu64(a->hi[j], x, y);
}

/* t1:t0 += sum of the lanes of the lo's + 2^52 * sum of the lanes of the
   hi's, the lanes being below 2^61 */
FLINT_FORCE_INLINE void
u52_dot_acc_dump(ulong * t1, ulong * t0, const u52_dot_acc * a)
{
    __m512i lo = _mm512_add_epi64(_mm512_add_epi64(a->lo[0], a->lo[1]),
                                  _mm512_add_epi64(a->lo[2], a->lo[3]));
    __m512i hi = _mm512_add_epi64(_mm512_add_epi64(a->hi[0], a->hi[1]),
                                  _mm512_add_epi64(a->hi[2], a->hi[3]));
    ulong slo = _mm512_reduce_add_epi64(lo);
    ulong shi = _mm512_reduce_add_epi64(hi);

    add_ssaaaa(*t1, *t0, *t1, *t0, shi >> 12, shi << 52);
    add_ssaaaa(*t1, *t0, *t1, *t0, UWORD(0), slo);
}

/*
    Chunks of U52_DOT_CHUNK unrolled iterations; after the last chunk, the
    remaining vectors go to the accumulators in turn, the last one masked.
*/
#define U52_DOT_BODY(LOAD2, LOAD2_MASKED)                                   \
    u52_dot_acc acc;                                                        \
    ulong t1 = 0, t0 = 0, res;                                              \
    slong i = 0, stop;                                                      \
    const slong full = len - len % U52_DOT_UNROLL;                          \
    const slong vecs = len - len % 8;                                       \
    int j;                                                                  \
                                                                            \
    u52_dot_acc_zero(&acc);                                                 \
                                                                            \
    while (i < full)                                                        \
    {                                                                       \
        stop = FLINT_MIN(full, i + U52_DOT_CHUNK * U52_DOT_UNROLL);         \
        for ( ; i < stop; i += U52_DOT_UNROLL)                              \
            for (j = 0; j < U52_DOT_NACC; j++)                              \
                u52_dot_acc_step(&acc, j,                                   \
                    _mm512_loadu_si512((const void *) (vec1 + i + 8 * j)),  \
                    LOAD2(i + 8 * j));                                      \
        if (i < full)                                                       \
        {                                                                   \
            u52_dot_acc_dump(&t1, &t0, &acc);                               \
            u52_dot_acc_zero(&acc);                                         \
        }                                                                   \
    }                                                                       \
                                                                            \
    for (j = 0; i < vecs; i += 8, j++)                                      \
        u52_dot_acc_step(&acc, j,                                           \
                _mm512_loadu_si512((const void *) (vec1 + i)), LOAD2(i));   \
                                                                            \
    if (i < len)                                                            \
    {                                                                       \
        const __mmask8 m = (__mmask8) ((1 << (len - i)) - 1);               \
        u52_dot_acc_step(&acc, j % U52_DOT_NACC,                            \
                _mm512_maskz_loadu_epi64(m, (const void *) (vec1 + i)),     \
                LOAD2_MASKED(i, m));                                        \
    }                                                                       \
                                                                            \
    u52_dot_acc_dump(&t1, &t0, &acc);                                       \
    NMOD2_RED2(res, t1, t0, mod);                                           \
    return res;

ulong
_nmod_vec_dot_u52(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
{
    /* the parameters may come from a longer length: short products are
       faster with the scalar code */
    if (len < NMOD_VEC_DOT_U52_MIN_LEN)
        return _nmod_vec_dot2(vec1, vec2, len, mod);

#define LOAD2(k) _mm512_loadu_si512((const void *) (vec2 + (k)))
#define LOAD2_MASKED(k, m) _mm512_maskz_loadu_epi64(m, (const void *) (vec2 + (k)))
    U52_DOT_BODY(LOAD2, LOAD2_MASKED)
#undef LOAD2
#undef LOAD2_MASKED
}

ulong
_nmod_vec_dot_u52_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
{
    if (len < NMOD_VEC_DOT_U52_MIN_LEN)
        return _nmod_vec_dot2_rev(vec1, vec2, len, mod);

    /* lane j of LOAD2(k) is vec2[len - 1 - k - j]; for the masked tail of
       c = len - k terms, lane j < c is vec2[c - 1 - j] (the lanes j >= c
       meet zeros of vec1) */
    const __m512i rev = _mm512_set_epi64(0, 1, 2, 3, 4, 5, 6, 7);
#define LOAD2(k) _mm512_permutexvar_epi64(rev, \
        _mm512_loadu_si512((const void *) (vec2 + len - 8 - (k))))
#define LOAD2_MASKED(k, m) _mm512_permutexvar_epi64( \
        _mm512_sub_epi64(_mm512_set1_epi64(len - (k) - 1), \
                         _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0)), \
        _mm512_maskz_loadu_epi64(m, (const void *) vec2))
    U52_DOT_BODY(LOAD2, LOAD2_MASKED)
#undef LOAD2
#undef LOAD2_MASKED
}

ulong
_nmod_vec_dot_u52_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset,
                      slong len, nmod_t mod)
{
    if (len < NMOD_VEC_DOT_U52_MIN_LEN)
        return _nmod_vec_dot2_ptr(vec1, vec2, offset, len, mod);

    /* vec2[k][offset]: the 8 row pointers, each plus 8*offset bytes */
    const void * base = (const void *) (offset * (slong) sizeof(ulong));
#define LOAD2(k) _mm512_i64gather_epi64( \
        _mm512_loadu_si512((const void *) (vec2 + (k))), base, 1)
#define LOAD2_MASKED(k, m) _mm512_mask_i64gather_epi64( \
        _mm512_setzero_si512(), m, \
        _mm512_maskz_loadu_epi64(m, (const void *) (vec2 + (k))), base, 1)
    U52_DOT_BODY(LOAD2, LOAD2_MASKED)
#undef LOAD2
#undef LOAD2_MASKED
}

#else

/* never selected by _nmod_vec_dot_params without AVX512-IFMA */
ulong
_nmod_vec_dot_u52(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
{
    return _nmod_vec_dot2(vec1, vec2, len, mod);
}

ulong
_nmod_vec_dot_u52_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
{
    return _nmod_vec_dot2_rev(vec1, vec2, len, mod);
}

ulong
_nmod_vec_dot_u52_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset,
                      slong len, nmod_t mod)
{
    return _nmod_vec_dot2_ptr(vec1, vec2, offset, len, mod);
}

#endif
