/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_MAT_MUL_U52_VEC_H
#define NMOD_MAT_MUL_U52_VEC_H

/*
    AVX512-IFMA primitives shared by nmod_mat_mul_u52 (mul_u52.c) and the
    vector-matrix products of nmod_vec_mul.c: the modulus-derived
    constants and the reduction of a pair of accumulators (lo, hi), which
    hold sums of low and high 52-bit halves of products, to the canonical
    residue of hi*2^52 + lo. See mul_u52.c for the error analysis, which
    holds for lo, hi <= 3n + KC n^2 with KC <= 511 (in particular for the
    sums of the halves of up to KC products of residues, plus a residue).
*/

#include <immintrin.h>
#include "longlong.h"
#include "nmod_mat/impl.h"

#if NMOD_MAT_HAVE_MUL_U52

typedef struct
{
    ulong n;
    ulong c52;      /* 2^52 mod n */
    ulong w;        /* floor(c52 * 2^52 / n), for the Shoup step */
    double ninv;
    int lo_only;    /* n <= 2^26: all high halves vanish */
}
u52_ctx_struct;

static inline void
u52_ctx_init(u52_ctx_struct * ctx, ulong n)
{
    ulong q, r;

    FLINT_ASSERT(n >= 1 && n <= (UWORD(1) << 52));

    ctx->n = n;
    ctx->ninv = 1.0 / (double) n;
    ctx->c52 = (UWORD(1) << 52) % n;

    /* c52 * 2^52 has its high limb c52 >> 12 < n, as udiv_qrnnd wants */
    udiv_qrnnd(q, r, ctx->c52 >> 12, ctx->c52 << 52, n);
    ctx->w = q;

    ctx->lo_only = (n <= (UWORD(1) << 26));
}

typedef struct
{
    __m512i nv;
    __m512i c52v;
    __m512i wv;
    __m512d ninvv;
}
u52_consts;

FLINT_FORCE_INLINE u52_consts
u52_consts_init(const u52_ctx_struct * ctx)
{
    u52_consts C;

    C.nv = _mm512_set1_epi64(ctx->n);
    C.c52v = _mm512_set1_epi64(ctx->c52);
    C.wv = _mm512_set1_epi64(ctx->w);
    C.ninvv = _mm512_set1_pd(ctx->ninv);

    return C;
}

#define u52_madd52lo(acc, a, b) _mm512_madd52lo_epu64(acc, a, b)
#define u52_madd52hi(acc, a, b) _mm512_madd52hi_epu64(acc, a, b)

/* canonical residue of x for 0 <= x <= 3n + KC n^2, KC <= 511 (x < 2^61) */
FLINT_FORCE_INLINE __m512i
u52_red(__m512i x, const u52_consts * C)
{
    __m512d d = _mm512_cvtepu64_pd(x);
    __m512d qd = _mm512_mul_pd(d, C->ninvv);
    __m512i q, r;
    __mmask8 m;

    q = _mm512_cvt_roundpd_epu64(qd, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    r = _mm512_sub_epi64(x, _mm512_mullo_epi64(q, C->nv));

    m = _mm512_cmplt_epi64_mask(r, _mm512_setzero_si512());
    r = _mm512_mask_add_epi64(r, m, r, C->nv);
    m = _mm512_cmpge_epi64_mask(r, C->nv);
    r = _mm512_mask_sub_epi64(r, m, r, C->nv);

    return r;
}

/* canonical residue of hi*2^52 + lo for lo, hi < 2^61 */
FLINT_FORCE_INLINE __m512i
u52_red2(__m512i lo, __m512i hi, const u52_consts * C)
{
    const __m512i zero = _mm512_setzero_si512();
    __m512i h, q, t1, t2, r;

    h = u52_red(hi, C);

    /* r = h*c52 - q*n in [0, 2n), each product as its exact 104 bits */
    q = u52_madd52hi(zero, h, C->wv);
    t1 = _mm512_add_epi64(u52_madd52lo(zero, h, C->c52v),
                          _mm512_slli_epi64(u52_madd52hi(zero, h, C->c52v), 52));
    t2 = _mm512_add_epi64(u52_madd52lo(zero, q, C->nv),
                          _mm512_slli_epi64(u52_madd52hi(zero, q, C->nv), 52));
    r = _mm512_sub_epi64(t1, t2);

    return u52_red(_mm512_add_epi64(r, lo), C);
}

#endif  /* NMOD_MAT_HAVE_MUL_U52 */

#endif
