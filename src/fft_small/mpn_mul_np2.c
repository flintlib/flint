/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "flint.h"
#include "ulong_extras.h"
#include "mpn_extras.h"
#include "longlong.h"
#include "fft_small.h"

#include "machine_vectors.h"

/*
    Two-prime multiplication for the Toom boundary (roughly 200-2000
    limbs), where the wide-chunk pipelines lose 40-50% of their time
    to input conversion and CRT. With 42-bit chunks:

    - packed chunks lie below both primes, so the input conversion is
      a prime-independent bit slicing with no spreading multiplies;
    - the CRT is a single Garner step, vectorized in 53-bit lanes,
      with the 2^-depth normalization folded into its constants;
    - the recomposition adds one < 2^100 value per chunk into the
      output at a 42-bit stride.

    Coefficients accumulate min(alen, blen) products of (2^42-1)^2,
    so exactness needs 84 + ceil(log2(min alen, blen)) below
    log2(p1 p2), giving about 20000 limbs of headroom for the
    smaller operand; callers stay far under this.
*/

/* smallest feasible depth, then the largest chunk size the
   exactness bound blen * 2^(2 bits) <= p1 p2 allows: larger chunks
   shorten the transforms and often drop a whole depth level near the
   staircase boundaries */
/* returns 0 when no chunk size satisfies the exactness bound
   len_bound * blen * 2^(2 bits) <= p1 p2 */
int
_fft_small_np2_choose(const mpn_ctx_struct * R, ulong an, ulong bn,
           ulong len_bound, ulong * bits, ulong * depth)
{
    ulong d, b, dmax = 30;

    for (d = LG_BLK_SZ; d <= dmax; d++)
    {
        for (b = 45; b >= 40; b--)
        {
            ulong alen = (64 * an + b - 1) / b;
            ulong blen = (64 * bn + b - 1) / b;
            ulong t, thi;
            if (alen + blen - 1 > (UWORD(1) << d))
                continue;
            umul_ppmm(thi, t, len_bound, blen);
            if (thi != 0 || 2 * b + FLINT_BIT_COUNT(t) > R->np2_prodbits)
                continue;
            *bits = b;
            *depth = d;
            return 1;
        }
    }
    return 0;
}

/*
    Packing variants, all writing both prime copies directly (8-22%
    faster than packing once and copying, the more so once the
    buffers leave L1):

    - byte-window (default): since (pos & 7) + bits <= 52, each chunk
      is one unaligned 64-bit load at byte offset pos >> 3. Fastest
      measured on x86-64 by 2x over the combine.
    - aligned combine (default on arm64, where unaligned loads can be
      slow): two aligned word loads with a shift-range-safe combine.
    - AVX2 gather (kept for comparison only): byte-scale vpgatherqq +
      vpsrlvq + magic-number conversion; 1.9x slower than the scalar
      byte loads on tested x86-64, since scattered unaligned loads
      are what scalar load ports do well and gathers do poorly there.
      May be worth revisiting on microarchitectures with fast gathers.
*/
#if defined(FLINT_NP2_PACK_GATHER) && defined(__AVX2__)

#include <immintrin.h>

void
_fft_small_np2_pack(double * d, double * d2, nn_srcptr x, ulong xn, ulong len,
         ulong itrunc, ulong bits)
{
    ulong j, jsafe, pos, w, sh, v;
    ulong mask = (UWORD(1) << bits) - 1;
    const unsigned char * xb = (const unsigned char *) x;
    __m256i vmask = _mm256_set1_epi64x(mask);
    __m256i magic = _mm256_set1_epi64x(0x4330000000000000UL);
    __m256d magicd = _mm256_castsi256_pd(magic);
    __m256i vbits4 = _mm256_set_epi64x(3 * bits, 2 * bits, bits, 0);

    jsafe = FLINT_MIN(64 * (xn - 1) / bits + 1, len);
    for (j = 0; j + 4 <= jsafe; j += 4)
    {
        __m256i p = _mm256_add_epi64(_mm256_set1_epi64x(j * bits), vbits4);
        __m256i off = _mm256_srli_epi64(p, 3);
        __m256i s = _mm256_and_si256(p, _mm256_set1_epi64x(7));
        __m256i vv = _mm256_i64gather_epi64((const long long *) xb, off, 1);
        __m256d e;
        vv = _mm256_and_si256(_mm256_srlv_epi64(vv, s), vmask);
        e = _mm256_sub_pd(_mm256_castsi256_pd(_mm256_or_si256(vv, magic)),
                          magicd);
        _mm256_storeu_pd(d + j, e);
        _mm256_storeu_pd(d2 + j, e);
    }
    for (; j < len; j++)
    {
        double e;
        pos = j * bits;
        w = pos >> 6;
        sh = pos & 63;
        v = x[w] >> sh;
        if (sh > 64 - bits && w + 1 < xn)
            v |= x[w + 1] << (64 - sh);
        e = (double) (v & mask);
        d[j] = e; d2[j] = e;
    }
    for (; j < itrunc; j++)
    {
        d[j] = 0.0; d2[j] = 0.0;
    }
}

#elif defined(__aarch64__) || defined(FLINT_NP2_PACK_ALIGNED)

void
_fft_small_np2_pack(double * d, double * d2, nn_srcptr x, ulong xn, ulong len,
         ulong itrunc, ulong bits)
{
    ulong j, jsafe, pos, w, sh, v;
    ulong mask = (UWORD(1) << bits) - 1;

    jsafe = (xn >= 2) ? FLINT_MIN((64 * (xn - 1) - 1) / bits + 1, len) : 0;
    for (j = 0; j < jsafe; j++)
    {
        double e;
        pos = j * bits;
        w = pos >> 6;
        sh = pos & 63;
        e = (double) ((((x[w] >> sh) |
                        ((x[w + 1] << (63 - sh)) << 1))) & mask);
        d[j] = e; d2[j] = e;
    }
    for (; j < len; j++)
    {
        double e;
        pos = j * bits;
        w = pos >> 6;
        sh = pos & 63;
        v = x[w] >> sh;
        if (sh > 64 - bits && w + 1 < xn)
            v |= x[w + 1] << (64 - sh);
        e = (double) (v & mask);
        d[j] = e; d2[j] = e;
    }
    for (; j < itrunc; j++)
    {
        d[j] = 0.0; d2[j] = 0.0;
    }
}

#else  /* byte-window default */

/* writes both prime copies directly: measured 8-22% faster than
   packing once and copying, the more so once the buffers leave L1 */
void
_fft_small_np2_pack(double * d, double * d2, nn_srcptr x, ulong xn, ulong len,
         ulong itrunc, ulong bits)
{
    ulong j, jsafe, pos, w, sh, v;
    ulong mask = (UWORD(1) << bits) - 1;
    const unsigned char * xb = (const unsigned char *) x;

    /* since (pos & 7) + bits <= 7 + 45 < 64, each chunk fits in one
       unaligned 64-bit load at byte offset pos >> 3, replacing the
       two-word combine entirely; measured twice as fast as the
       combine, while a gather-based vector version is twice as slow */
    jsafe = (xn >= 1) ? 64 * (xn - 1) / bits + 1 : 0;
    jsafe = FLINT_MIN(jsafe, len);

    for (j = 0; j + 4 <= jsafe; j += 4)
    {
        ulong v0, v1, v2, v3;
        ulong p0 = (j + 0) * bits, p1 = (j + 1) * bits;
        ulong p2 = (j + 2) * bits, p3 = (j + 3) * bits;
        memcpy(&v0, xb + (p0 >> 3), 8);
        memcpy(&v1, xb + (p1 >> 3), 8);
        memcpy(&v2, xb + (p2 >> 3), 8);
        memcpy(&v3, xb + (p3 >> 3), 8);
        double e0 = (double) ((v0 >> (p0 & 7)) & mask);
        double e1 = (double) ((v1 >> (p1 & 7)) & mask);
        double e2 = (double) ((v2 >> (p2 & 7)) & mask);
        double e3 = (double) ((v3 >> (p3 & 7)) & mask);
        d[j + 0] = e0; d2[j + 0] = e0;
        d[j + 1] = e1; d2[j + 1] = e1;
        d[j + 2] = e2; d2[j + 2] = e2;
        d[j + 3] = e3; d2[j + 3] = e3;
    }
    for (; j < jsafe; j++)
    {
        double e;
        pos = j * bits;
        memcpy(&v, xb + (pos >> 3), 8);
        e = (double) ((v >> (pos & 7)) & mask);
        d[j] = e; d2[j] = e;
    }
    for (; j < len; j++)
    {
        double e;
        pos = j * bits;
        w = pos >> 6;
        sh = pos & 63;
        v = x[w] >> sh;
        if (sh > 64 - bits && w + 1 < xn)
            v |= x[w + 1] << (64 - sh);
        e = (double) (v & mask);
        d[j] = e; d2[j] = e;
    }
    for (; j < itrunc; j++)
    {
        d[j] = 0.0; d2[j] = 0.0;
    }
}

#endif

/* z[0 .. zh - zl) = limbs [zl, zh) of a*b. Exact for zl = 0. For
   zl > 0 a lower approximation at bit granularity: chunks straddling
   bit 64 zl contribute their bits at or above it (the below-window
   part masked off), and only chunks lying entirely below are
   excluded, so the deficit is the carry of a sum of values each
   below 2^(64 zl) - under two ulps of limb zl, within the range
   driver's documented min(an, bn, lo) bound. */

/*
    The window of _mpn_ctx_mpn_mul_window_np2 as a transposed product
    (see _nmod_poly_mul_tr_mpn_ctx), when that is shorter. The
    recomposition reads the chunks [j0, jtop) of the product, with j0 the
    first chunk that can reach bit 64 zl (the straddlers) and jtop the
    first one above limb zh, and nothing else; these coefficients are
    computed exactly as

        c[c0 + j] = sum_i r[i] v[i + j],    0 <= j < m = jtop - c0,

    with r the sl chunks of b reversed (b's zero chunks above blen
    included) and v the chunks of a from s0 = c0 - (sl - 1) on, at the
    cost of a product of lengths sl and m instead of alen and blen. With
    c0 * bits and (sl - 1) * bits multiples of 64 (and c0 one of 8, for
    the vectorized Garner step), v starts at limb s0 * bits / 64 of a (or
    at zero chunks in front of it), and the recomposition runs on the
    shifted chunks as on the full product: same chunks, same chunk
    positions relative to the limbs, so the result is the same, bit for
    bit, as that of the window of the full product.

    Returns 0, doing nothing, if the transposed product is not shorter.
*/
static int
_mpn_ctx_mpn_mul_window_np2_tr(mpn_ctx_t R, nn_ptr z, nn_ptr zlow, ulong nlow,
                        ulong zl, ulong zh,
                        nn_srcptr a, ulong an, nn_srcptr b, ulong bn,
                        ulong bits, ulong alen, ulong blen, ulong zlen)
{
    ulong g, q, jstr, jtop, c0, sl, m, pn, ztr, btr, otr, depth, dsz, j;
    slong s0;
    sd_fft_ctx_struct * Q0 = R->ffts + 0, * Q1 = R->ffts + 1;
    double * d0, * d1, * e0, * e1;

    /* g * bits is the smallest multiple of bits divisible by 64 */
    g = 64;
    while (g > 1 && (g / 2) * bits % 64 == 0)
        g /= 2;
    q = FLINT_MAX(g, 8);

    jstr = (64 * zl > 99) ? (64 * zl - 99) / bits : 0;
    jtop = FLINT_MIN(zlen, 64 * zh / bits + 1);
    c0 = jstr - jstr % q;
    m = jtop - c0;
    sl = blen + (g - (blen - 1) % g) % g;
    pn = sl + m - 1;

    FLINT_ASSERT(jstr < jtop);

    if (m >= alen || n_round_up(pn, BLK_SZ) >= n_round_up(zlen, BLK_SZ))
        return 0;

    ztr = n_round_up(pn, BLK_SZ);
    btr = n_round_up(sl, BLK_SZ);
    otr = n_round_up(m, BLK_SZ);
    depth = FLINT_MAX(LG_BLK_SZ, n_clog2(ztr));
    dsz = sd_fft_ctx_data_size(depth);
    s0 = (slong) c0 - (slong) (sl - 1);

    FLINT_ASSERT(c0 * bits % 64 == 0 && (sl - 1) * bits % 64 == 0);
    FLINT_ASSERT(s0 < (slong) alen);

    d0 = mpn_ctx_fit_buffer(R, 4 * dsz * sizeof(double));
    d1 = d0 + dsz;
    e0 = d1 + dsz;
    e1 = e0 + dsz;

    /* v: zero chunks in front of a if s0 < 0, otherwise a from limb
       s0 * bits / 64; zeros up to ztr (the chunks beyond pn meet only
       zero rows of the transposed product) */
    if (s0 < 0)
    {
        ulong k = -s0;

        for (j = 0; j < k; j++)
        {
            d0[j] = 0.0;
            d1[j] = 0.0;
        }

        _fft_small_np2_pack(d0 + k, d1 + k, a, an,
                            FLINT_MIN(alen, pn - k), ztr - k, bits);
    }
    else
    {
        ulong off = (ulong) s0 * bits / 64;
        ulong xn = an - off;

        _fft_small_np2_pack(d0, d1, a + off, xn,
                            FLINT_MIN(pn, (64 * xn + bits - 1) / bits), ztr, bits);
    }

    sd_ifft_trunc_tr(Q0, d0, depth, ztr);
    sd_ifft_trunc_tr(Q1, d1, depth, ztr);

    /* r: the chunks of b, zero padded to sl, reversed */
    _fft_small_np2_pack(e0, e1, b, bn, blen, btr, bits);
    for (j = 0; j < sl / 2; j++)
    {
        FLINT_SWAP(double, e0[j], e0[sl - 1 - j]);
        FLINT_SWAP(double, e1[j], e1[sl - 1 - j]);
    }

    sd_fft_trunc(Q0, e0, depth, btr, ztr);
    sd_fft_trunc(Q1, e1, depth, btr, ztr);

    for (int k = 0; k < 2; k++)
    {
        sd_fft_ctx_struct * Q = k ? Q1 : Q0;
        vec8d n8 = vec8d_set_d(Q->p), ninv8 = vec8d_set_d(Q->pinv);
        double * X = k ? d1 : d0, * Y = k ? e1 : e0;
        for (j = 0; j < ztr; j += 8)
            vec8d_store(X + j, vec8d_mulmod(vec8d_load(X + j),
                                            vec8d_load(Y + j), n8, ninv8));
    }

    /* the transposes of the transforms are unnormalized like the
       transforms, so the recomposition normalizes by 2^-depth as for
       the window of the full product */
    sd_fft_trunc_tr(Q0, d0, depth, otr, ztr);
    sd_fft_trunc_tr(Q1, d1, depth, otr, ztr);

    _fft_small_np2_crt_recompose_split(R, z, zlow, nlow, zl - c0 * bits / 64,
                                 zh - c0 * bits / 64, d0, d1, m, otr,
                                 bits, depth);

    return 1;
}

void
_mpn_ctx_mpn_mul_window_np2(mpn_ctx_t R, nn_ptr z, ulong zl, ulong zh,
                        nn_srcptr a, ulong an, nn_srcptr b, ulong bn)
{
    _mpn_ctx_mpn_mul_window_np2_split(R, z, NULL, 0, zl, zh, a, an, b, bn);
}

/* As above, with the lowest nlow limbs of the window written to zlow and
   the remaining zh - zl - nlow limbs to z. */
void
_mpn_ctx_mpn_mul_window_np2_split(mpn_ctx_t R, nn_ptr z, nn_ptr zlow,
                        ulong nlow, ulong zl, ulong zh,
                        nn_srcptr a, ulong an, nn_srcptr b, ulong bn)
{
    ulong bits, depth;
    ulong alen, blen, zlen, zn = zh;

    FLINT_ASSERT(zl < zh && zh <= an + bn);

    {
        int ok = _fft_small_np2_choose(R, an, bn, 1, &bits, &depth);
        FLINT_ASSERT(ok);
        (void) ok;
    }
    alen = (64 * an + bits - 1) / bits;
    blen = (64 * bn + bits - 1) / bits;
    zlen = alen + blen - 1;
    ulong atr = n_round_up(alen, BLK_SZ);
    ulong btr = n_round_up(blen, BLK_SZ);
    ulong ztr = n_round_up(zlen, BLK_SZ);
    int sqr = (a == b && an == bn);
    sd_fft_ctx_struct * Q0 = R->ffts + 0, * Q1 = R->ffts + 1;
    ulong dsz = sd_fft_ctx_data_size(depth);
    double * d0, * d1, * e0, * e1;
    ulong j;

    FLINT_ASSERT(an >= bn && bn >= 1);

    if (!sqr && _mpn_ctx_mpn_mul_window_np2_tr(R, z, zlow, nlow, zl, zh,
                                              a, an, b, bn,
                                              bits, alen, blen, zlen))
        return;

    d0 = mpn_ctx_fit_buffer(R, (sqr ? 2 : 4) * dsz * sizeof(double));
    d1 = d0 + dsz;
    e0 = d1 + dsz;
    e1 = e0 + dsz;

    _fft_small_np2_pack(d0, d1, a, an, alen, atr, bits);
    sd_fft_trunc(Q0, d0, depth, atr, ztr);
    sd_fft_trunc(Q1, d1, depth, atr, ztr);

    if (!sqr)
    {
        _fft_small_np2_pack(e0, e1, b, bn, blen, btr, bits);
        sd_fft_trunc(Q0, e0, depth, btr, ztr);
        sd_fft_trunc(Q1, e1, depth, btr, ztr);
    }
    else
    {
        e0 = d0;
        e1 = d1;
    }

    for (int k = 0; k < 2; k++)
    {
        sd_fft_ctx_struct * Q = k ? Q1 : Q0;
        vec8d n8 = vec8d_set_d(Q->p), ninv8 = vec8d_set_d(Q->pinv);
        double * X = k ? d1 : d0, * Y = k ? e1 : e0;
        for (j = 0; j < ztr; j += 8)
            vec8d_store(X + j, vec8d_mulmod(vec8d_load(X + j),
                                            vec8d_load(Y + j), n8, ninv8));
    }

    sd_ifft_trunc(Q0, d0, depth, ztr);
    sd_ifft_trunc(Q1, d1, depth, ztr);

    _fft_small_np2_crt_recompose_split(R, z, zlow, nlow, zl, zn, d0, d1,
                                       zlen, (ulong) ztr, bits, depth);
}

/* limb i of a window whose lowest nlow limbs go to zlow and the rest to z */
static inline void
_np2_write_limb(nn_ptr z, nn_ptr zlow, ulong nlow, ulong i, ulong val)
{
    if (i < nlow)
        zlow[i] = val;
    else
        z[i - nlow] = val;
}

/* Garner and register-window recomposition of two-prime lanes
   d0, d1 (writable; consumed in place): z[0 .. zh - zl) receives
   limbs [zl, zh) of sum_j v_j 2^(j bits) over nchunks chunk values,
   with straddlers of bit 64 zl masked to bit granularity. For
   depth != UWORD_MAX the lanes carry raw inverse-transform values
   and are normalized by 2^-depth during Garner; for
   depth == UWORD_MAX they are already plain residues (the
   negacyclic path, whose untwist folds the normalization). Shared by
   the multiplication window, the op export and the negacyclic
   reconstruction; the signed export keeps its own centered variant. */
/* Garner reconstruction of the two-prime lanes into the limb window
   [zl, zh) of z. Destructive: the lane arrays d0, d1 are consumed as
   scratch by the vectorized sweep. Callers hold either the product
   buffers of a multiplication in progress or a transform already
   sacrificed to a destructive conversion; non-destructive conversions
   copy first and call the fast destructive code, the default
   semantics of output conversion throughout fft_small. */
void
_fft_small_np2_crt_recompose(const mpn_ctx_struct * R, nn_ptr z,
    ulong zl, ulong zh, double * d0, double * d1, ulong nchunks,
    ulong navail8, ulong bits, ulong depth)
{
    _fft_small_np2_crt_recompose_split(R, z, NULL, 0, zl, zh, d0, d1,
                                       nchunks, navail8, bits, depth);
}

/* As above, with the lowest nlow limbs of the window, [zl, zl + nlow),
   written to zlow instead of z, which receives [zl + nlow, zh). */
void
_fft_small_np2_crt_recompose_split(const mpn_ctx_struct * R, nn_ptr z,
    nn_ptr zlow, ulong nlow, ulong zl, ulong zh, double * d0, double * d1,
    ulong nchunks, ulong navail8, ulong bits, ulong depth)
{
    const sd_fft_ctx_struct * Q0 = R->ffts + 0, * Q1 = R->ffts + 1;
    ulong p1u = (ulong) Q0->p;
    ulong j;
    vec8d P1 = vec8d_set_d(Q0->p), P1i = vec8d_set_d(Q0->pinv);
    vec8d P2 = vec8d_set_d(Q1->p), P2i = vec8d_set_d(Q1->pinv);
    vec8d C = vec8d_set_d((double) R->np2_ip1);

    {
        ulong jgar = FLINT_MIN(navail8,
                               n_round_up(64 * zh / bits + 1, 8));
        ulong j0 = (64 * zl > 99)
                 ? (((64 * zl - 99) / bits) & ~(ulong) 7) : 0;

        if (depth != UWORD_MAX)
        {
            vec8d S1 = vec8d_set_d((double) R->np2_s1[depth]);
            vec8d S2 = vec8d_set_d((double) R->np2_s2[depth]);
            for (j = j0; j < jgar; j += 8)
            {
                vec8d r1 = vec8d_mulmod(vec8d_load(d0 + j), S1, P1, P1i);
                vec8d r2 = vec8d_mulmod(vec8d_load(d1 + j), S2, P2, P2i);
                vec8d t;
                r1 = vec8d_reduce_to_0n(r1, P1, P1i);
                r2 = vec8d_reduce_to_0n(r2, P2, P2i);
                t = vec8d_mulmod(vec8d_add(vec8d_sub(r2, r1), P2),
                                 C, P2, P2i);
                t = vec8d_reduce_to_0n(t, P2, P2i);
                vec8d_store(d0 + j, r1);
                vec8d_store(d1 + j, t);
            }
        }
        else
        {
            for (j = j0; j < jgar; j += 8)
            {
                vec8d r1 = vec8d_reduce_to_0n(vec8d_load(d0 + j), P1, P1i);
                vec8d r2 = vec8d_reduce_to_0n(vec8d_load(d1 + j), P2, P2i);
                vec8d t = vec8d_mulmod(vec8d_add(vec8d_sub(r2, r1), P2),
                                       C, P2, P2i);
                t = vec8d_reduce_to_0n(t, P2, P2i);
                vec8d_store(d0 + j, r1);
                vec8d_store(d1 + j, t);
            }
        }
    }

    {
        ulong acc0 = 0, acc1 = 0, acc2 = 0, acc3 = 0;
        ulong jtop = FLINT_MIN(nchunks, 64 * zh / bits + 1);
        ulong jlo = (zl == 0) ? 0 : (64 * zl + bits - 1) / bits;
        ulong jstr = (64 * zl > 99) ? (64 * zl - 99) / bits : 0;
        ulong w0 = zl;

        for (j = jstr; j < FLINT_MIN(jlo, jtop); j++)
        {
            ulong r1 = (ulong) d0[j], tt = (ulong) d1[j];
            ulong vhi, vlo, sh2;
            ulong pos = j * bits;

            if (pos + 100 <= 64 * zl)
                continue;

            umul_ppmm(vhi, vlo, p1u, tt);
            add_ssaaaa(vhi, vlo, vhi, vlo, 0, r1);
            sh2 = 64 * zl - pos;
            if (sh2 >= 64)
            {
                vlo = vhi >> (sh2 - 64);
                vhi = 0;
            }
            else if (sh2 > 0)
            {
                vlo = (vlo >> sh2) | (vhi << (64 - sh2));
                vhi >>= sh2;
            }
            add_ssaaaa(acc1, acc0, acc1, acc0, vhi, vlo);
        }

        /* the diverted limbs [zl, zl + nlow) first, with the branchy
           writer (a few iterations), then the rest through the shifted
           base zbase = z - nlow - zl without any test per limb */
        nn_ptr zbase = z - nlow - zl;

        for (j = jlo; j < jtop && w0 < zl + nlow; j++)
        {
            ulong r1 = (ulong) d0[j], tt = (ulong) d1[j];
            ulong vhi, vlo, x0, x1, x2;
            ulong pos = j * bits, w = pos >> 6, sh = pos & 63;

            while (w0 < w)
            {
                _np2_write_limb(z, zlow, nlow, w0 - zl, acc0);
                w0++;
                acc0 = acc1; acc1 = acc2; acc2 = acc3; acc3 = 0;
            }

            umul_ppmm(vhi, vlo, p1u, tt);
            add_ssaaaa(vhi, vlo, vhi, vlo, 0, r1);

            x0 = vlo << sh;
            x1 = ((vlo >> (63 - sh)) >> 1) | (vhi << sh);
            x2 = (vhi >> (63 - sh)) >> 1;

            add_ssssaaaaaaaa(acc3, acc2, acc1, acc0,
                             acc3, acc2, acc1, acc0, 0, x2, x1, x0);
        }
        while (w0 < zh && w0 < zl + nlow)
        {
            _np2_write_limb(z, zlow, nlow, w0 - zl, acc0);
            w0++;
            acc0 = acc1; acc1 = acc2; acc2 = acc3; acc3 = 0;
        }

        for ( ; j < jtop; j++)
        {
            ulong r1 = (ulong) d0[j], tt = (ulong) d1[j];
            ulong vhi, vlo, x0, x1, x2;
            ulong pos = j * bits, w = pos >> 6, sh = pos & 63;

            while (w0 < w)
            {
                zbase[w0++] = acc0;
                acc0 = acc1; acc1 = acc2; acc2 = acc3; acc3 = 0;
            }

            umul_ppmm(vhi, vlo, p1u, tt);
            add_ssaaaa(vhi, vlo, vhi, vlo, 0, r1);

            x0 = vlo << sh;
            x1 = ((vlo >> (63 - sh)) >> 1) | (vhi << sh);
            x2 = (vhi >> (63 - sh)) >> 1;

            add_ssssaaaaaaaa(acc3, acc2, acc1, acc0,
                             acc3, acc2, acc1, acc0, 0, x2, x1, x0);
        }
        while (w0 < zh)
        {
            zbase[w0++] = acc0;
            acc0 = acc1; acc1 = acc2; acc2 = acc3; acc3 = 0;
        }
    }
}

void
_mpn_ctx_mpn_mullow_np2(mpn_ctx_t R, nn_ptr z, ulong zh, nn_srcptr a,
                        ulong an, nn_srcptr b, ulong bn)
{
    _mpn_ctx_mpn_mul_window_np2(R, z, 0, zh, a, an, b, bn);
}

void
_mpn_ctx_mpn_mul_np2(mpn_ctx_t R, nn_ptr z, nn_srcptr a, ulong an,
                     nn_srcptr b, ulong bn)
{
    _mpn_ctx_mpn_mullow_np2(R, z, an + bn, a, an, b, bn);
}
