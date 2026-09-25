/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Dot products for moduli below 2^50 in double precision, the strategy
    of nmod_mat_mul_fp50 (the mulmod of fft_small, see nmod_mat/mul_fp_vec.h
    for the primitives and their bounds): every product is reduced to
    (-9/8 n, 9/8 n) as it is formed, four such terms are accumulated in a
    lane and the accumulator is brought back to (-0.51 n, 0.51 n), so that
    it never exceeds 5.01 n < 2^53 in absolute value and every operation is
    exact. Four accumulator vectors hide the latency of the fused
    multiply-adds. Meant for machines without AVX512-IFMA (AVX2, AVX-512
    without IFMA, NEON), where the alternative is the scalar two-limb code.

    The three variants gather the entries of vec2 in order, reversed, or
    through row pointers; the last two go through a small buffer.
*/

#include "nmod.h"
#include "nmod_vec.h"

#if NMOD_VEC_HAVE_DOT_FP50

#include "nmod_mat/mul_fp_vec.h"

#define FP50_DOT_NACC 4
#define FP50_DOT_UNROLL (FPV_VL * FP50_DOT_NACC)

/* the horizontal sum of a vector in (-0.51 n, 0.51 n), as an integer
   offset by 8 n so that it is nonnegative (n < 2^50: no overflow) */
FLINT_FORCE_INLINE ulong
fp50_dot_hsum(fpv acc, ulong n)
{
    double buf[FPV_VL];
    double s = 0.0;
    int l;

    fpv_storeu(buf, acc);
    for (l = 0; l < FPV_VL; l++)
        s += buf[l];

    return (ulong) ((slong) s + 8 * (slong) n);
}

#define FP50_DOT_BODY(LOAD2)                                                \
    const ulong n = mod.n;                                                  \
    const fpv nv = fpv_set1((double) n);                                    \
    const fpv ninv = fpv_set1(1.0 / (double) n);                            \
    fpv acc[FP50_DOT_NACC], r;                                              \
    ulong t1, t0, res;                                                      \
    slong i = 0, full;                                                      \
    int j, cnt = 0;                                                         \
                                                                            \
    for (j = 0; j < FP50_DOT_NACC; j++)                                     \
        acc[j] = fpv_zero();                                                \
                                                                            \
    /* one term per accumulator and iteration, a reduction every four */    \
    full = len - len % FP50_DOT_UNROLL;                                     \
    for ( ; i < full; i += FP50_DOT_UNROLL)                                 \
    {                                                                       \
        for (j = 0; j < FP50_DOT_NACC; j++)                                 \
            acc[j] = fpv_add(acc[j], fpv_mulmod(fpv_load_u64(vec1 + i + j * FPV_VL), \
                                                LOAD2(i + j * FPV_VL), nv, ninv)); \
        if (++cnt == 4)                                                     \
        {                                                                   \
            for (j = 0; j < FP50_DOT_NACC; j++)                             \
                acc[j] = fpv_reduce_pm1n(acc[j], nv, ninv);                 \
            cnt = 0;                                                        \
        }                                                                   \
    }                                                                       \
    if (cnt != 0)                                                           \
        for (j = 0; j < FP50_DOT_NACC; j++)                                 \
            acc[j] = fpv_reduce_pm1n(acc[j], nv, ninv);                     \
                                                                            \
    /* remaining full vectors, one at a time, reduced every time */         \
    full = len - len % FPV_VL;                                              \
    for ( ; i < full; i += FPV_VL)                                          \
        acc[0] = fpv_reduce_pm1n(fpv_add(acc[0],                            \
            fpv_mulmod(fpv_load_u64(vec1 + i), LOAD2(i), nv, ninv)),        \
            nv, ninv);                                                      \
                                                                            \
    r = fpv_reduce_pm1n(fpv_add(acc[0], acc[1]), nv, ninv);                 \
    r = fpv_reduce_pm1n(fpv_add(r, acc[2]), nv, ninv);                      \
    r = fpv_reduce_pm1n(fpv_add(r, acc[3]), nv, ninv);                      \
                                                                            \
    t1 = 0;                                                                 \
    t0 = fp50_dot_hsum(r, n);                                               \
    for ( ; i < len; i++)                                                   \
    {                                                                       \
        ulong s1, s0;                                                       \
        umul_ppmm(s1, s0, vec1[i], TAIL2(i));                               \
        add_ssaaaa(t1, t0, t1, t0, s1, s0);                                 \
    }                                                                       \
                                                                            \
    NMOD2_RED2(res, t1, t0, mod);                                           \
    return res;

ulong
_nmod_vec_dot_fp50(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
{
    /* the parameters may come from a longer length: short products are
       faster with the scalar code */
    if (len < NMOD_VEC_DOT_FP50_MIN_LEN)
        return _nmod_vec_dot2(vec1, vec2, len, mod);

#define LOAD2(k) fpv_load_u64(vec2 + (k))
#define TAIL2(k) vec2[k]
    FP50_DOT_BODY(LOAD2)
#undef LOAD2
#undef TAIL2
}

FLINT_FORCE_INLINE fpv
fp50_dot_load_rev(nn_srcptr p)
{
    ulong buf[FPV_VL];
    int l;
    for (l = 0; l < FPV_VL; l++)
        buf[l] = p[-l];
    return fpv_load_u64(buf);
}

ulong
_nmod_vec_dot_fp50_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
{
    /* the parameters may come from a longer length: short products are
       faster with the scalar code */
    if (len < NMOD_VEC_DOT_FP50_MIN_LEN)
        return _nmod_vec_dot2_rev(vec1, vec2, len, mod);

#define LOAD2(k) fp50_dot_load_rev(vec2 + len - 1 - (k))
#define TAIL2(k) vec2[len - 1 - (k)]
    FP50_DOT_BODY(LOAD2)
#undef LOAD2
#undef TAIL2
}

FLINT_FORCE_INLINE fpv
fp50_dot_load_ptr(const nn_ptr * p, slong offset)
{
    ulong buf[FPV_VL];
    int l;
    for (l = 0; l < FPV_VL; l++)
        buf[l] = p[l][offset];
    return fpv_load_u64(buf);
}

ulong
_nmod_vec_dot_fp50_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset,
                       slong len, nmod_t mod)
{
    /* the parameters may come from a longer length: short products are
       faster with the scalar code */
    if (len < NMOD_VEC_DOT_FP50_MIN_LEN)
        return _nmod_vec_dot2_ptr(vec1, vec2, offset, len, mod);

#define LOAD2(k) fp50_dot_load_ptr(vec2 + (k), offset)
#define TAIL2(k) vec2[k][offset]
    FP50_DOT_BODY(LOAD2)
#undef LOAD2
#undef TAIL2
}

#else

/* never selected by _nmod_vec_dot_params without a vector backend */
ulong
_nmod_vec_dot_fp50(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
{
    return _nmod_vec_dot2(vec1, vec2, len, mod);
}

ulong
_nmod_vec_dot_fp50_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
{
    return _nmod_vec_dot2_rev(vec1, vec2, len, mod);
}

ulong
_nmod_vec_dot_fp50_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset,
                       slong len, nmod_t mod)
{
    return _nmod_vec_dot2_ptr(vec1, vec2, offset, len, mod);
}

#endif
