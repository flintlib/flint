/*
    Copyright (C) 2026 Brian Heckel

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Internal (not installed) helpers for batched AVX2 Montgomery arithmetic on
    4x4 matrices over Z/nZ (n odd, n < 2^246), using LAZY reduction: each output
    entry sums its 4 raw products into a wide accumulator and reduces once.

    A value is s = ceil((bits+2)/31) limbs of 31 bits (2 guard bits so the
    single-subtract lazy REDC is valid for n < 2^(31s-2)).  The matmul takes
    s as a COMPILE-TIME argument and is force-inlined; bm_matmul dispatches on
    ctx->s through a switch so each size gets a fully-unrolled, right-sized
    instantiation (a runtime-s loop does not unroll and was measurably slower).
    Arrays are sized to the s=8 (248-bit) maximum; only limbs [0,s) are touched.
    AVX2-only; guard the use site with #ifdef __AVX2__.
*/

#ifndef MPN_MOD_BATCHED_MONT_H
#define MPN_MOD_BATCHED_MONT_H

#include "mpn_mod.h"
#include "gr_mat.h"
#include "fmpz.h"

#ifdef __AVX2__

#include <immintrin.h>

#if defined(__GNUC__)
#define BM_FORCE_INLINE static inline __attribute__((always_inline))
#else
#define BM_FORCE_INLINE static inline
#endif

#define BM_BITS     31
#define BM_MAXLIMB  8                               /* max 31-bit limbs (248-bit) */
#define BM_DIM      4
#define BM_MASK     (((slong) 1 << BM_BITS) - 1)

typedef struct
{
    __m256i n[BM_MAXLIMB];
    __m256i r2[BM_MAXLIMB];
    __m256i n0inv;
    __m256i mask;
    slong   s;
}
bm_ctx_t;

typedef ulong bm_mat[BM_DIM][BM_DIM][BM_MAXLIMB];

static inline void
bm_pack31(ulong * out, nn_srcptr a, slong nlimbs, slong s)
{
    slong o;
    for (o = 0; o < s; o++)
    {
        ulong p = (ulong) o * BM_BITS;
        slong w = (slong) (p >> 6);
        int off = (int) (p & 63);
        ulong bits = 0;
        if (w < nlimbs)
            bits = a[w] >> off;
        if (off + BM_BITS > 64 && w + 1 < nlimbs)
            bits |= a[w + 1] << (64 - off);
        out[o] = bits & BM_MASK;
    }
}

static inline void
bm_unpack31(nn_ptr a, const ulong * in, slong nlimbs, slong s)
{
    slong o, w;
    for (w = 0; w < nlimbs; w++)
        a[w] = 0;
    for (o = 0; o < s; o++)
    {
        ulong p = (ulong) o * BM_BITS;
        ulong v = in[o];
        w = (slong) (p >> 6);
        int off = (int) (p & 63);
        if (w < nlimbs)
            a[w] |= v << off;
        if (off + BM_BITS > 64 && w + 1 < nlimbs)
            a[w + 1] |= v >> (64 - off);
    }
}

static inline ulong
bm_n0inv(ulong n0)
{
    ulong inv = n0;
    inv *= 2 - n0 * inv;
    inv *= 2 - n0 * inv;
    inv *= 2 - n0 * inv;
    inv *= 2 - n0 * inv;
    return ((ulong) (-inv)) & BM_MASK;
}

static inline void
bm_ctx_init(bm_ctx_t * C, gr_ctx_t ctx)
{
    slong nlimbs = MPN_MOD_CTX_NLIMBS(ctx);
    nn_srcptr nd = MPN_MOD_CTX_MODULUS(ctx);
    slong bits = MPN_MOD_CTX_MODULUS_BITS(ctx);
    slong s = (bits + 2 + BM_BITS - 1) / BM_BITS;
    ulong n31[BM_MAXLIMB], r2nn[BM_MAXLIMB], r231[BM_MAXLIMB];
    fmpz_t N, R2;
    slong l;

    C->s = s;
    bm_pack31(n31, nd, nlimbs, s);
    for (l = 0; l < s; l++)
        C->n[l] = _mm256_set1_epi64x(n31[l]);
    C->n0inv = _mm256_set1_epi64x(bm_n0inv(n31[0]));
    C->mask  = _mm256_set1_epi64x(BM_MASK);

    fmpz_init(N);
    fmpz_init(R2);
    fmpz_set_ui_array(N, nd, nlimbs);
    fmpz_one(R2);
    fmpz_mul_2exp(R2, R2, 2 * BM_BITS * s);
    fmpz_mod(R2, R2, N);
    for (l = 0; l < nlimbs; l++)
        r2nn[l] = 0;
    fmpz_get_ui_array(r2nn, nlimbs, R2);
    bm_pack31(r231, r2nn, nlimbs, s);
    for (l = 0; l < s; l++)
        C->r2[l] = _mm256_set1_epi64x(r231[l]);
    fmpz_clear(N);
    fmpz_clear(R2);
}

/*
   Shared building blocks (all force-inlined; s is a compile-time constant at the
   hot call sites, so every loop below unrolls to a size-s specialization).
   Each lane holds an independent number; carries ripple *within* a lane across
   the limb array, so there are no cross-lane permutes.
*/

/*
   T[off .. off+s] += a[0..s-1] * b, propagating 31-bit carries.  The top limb
   ACCUMULATES (T[off+s] += carry) rather than being overwritten, so several
   calls at different offsets sum into one wide accumulator (lazy reduction).
*/
BM_FORCE_INLINE void
bm_fma_column(__m256i * T, const __m256i * a, __m256i b, __m256i mask,
              int off, const int s)
{
    __m256i carry = _mm256_setzero_si256();
    int j;
    for (j = 0; j < s; j++)
    {
        __m256i p = _mm256_add_epi64(_mm256_add_epi64(T[off + j],
                          _mm256_mul_epu32(a[j], b)), carry);
        T[off + j] = _mm256_and_si256(p, mask);
        carry = _mm256_srli_epi64(p, BM_BITS);
    }
    T[off + s] = _mm256_add_epi64(T[off + s], carry);
}

/* Ripple 31-bit carries left across T[0..len], leaving T[0..len-1] < 2^31. */
BM_FORCE_INLINE void
bm_normalize(__m256i * T, __m256i mask, const int len)
{
    int j;
    for (j = 0; j < len; j++)
    {
        __m256i carry = _mm256_srli_epi64(T[j], BM_BITS);
        T[j] = _mm256_and_si256(T[j], mask);
        T[j + 1] = _mm256_add_epi64(T[j + 1], carry);
    }
}

/*
   One Montgomery reduction step: cancel limb T[pos] by adding m*n, where
   m = T[pos] * n0inv (mod 2^31), then ripple the carry up through T[.. top].
*/
BM_FORCE_INLINE void
bm_redc_step(__m256i * T, int pos, const bm_ctx_t * bc, __m256i mask,
             const int s, int top)
{
    __m256i m = _mm256_and_si256(_mm256_mul_epu32(T[pos], bc->n0inv), mask);
    __m256i carry = _mm256_setzero_si256();
    int j;
    for (j = 0; j < s; j++)
    {
        __m256i p = _mm256_add_epi64(_mm256_add_epi64(T[pos + j],
                          _mm256_mul_epu32(m, bc->n[j])), carry);
        T[pos + j] = _mm256_and_si256(p, mask);
        carry = _mm256_srli_epi64(p, BM_BITS);
    }
    for (j = pos + s; j <= top; j++)
    {
        __m256i p = _mm256_add_epi64(T[j], carry);
        T[j] = _mm256_and_si256(p, mask);
        carry = _mm256_srli_epi64(p, BM_BITS);
    }
}

/*
   Final conditional subtraction.  hi[0..s-1] is the tentative result and
   overflow its (s+1)-th limb; subtract n once when hi >= n (i.e. no borrow, or
   the overflow limb is set), writing the canonical residue in [0,n) to out.
*/
BM_FORCE_INLINE void
bm_final_reduce(__m256i * out, const __m256i * hi, __m256i overflow,
                const bm_ctx_t * bc, __m256i mask, const int s)
{
    __m256i d[BM_MAXLIMB];
    __m256i borrow = _mm256_setzero_si256();
    __m256i take_d;
    int j;
    for (j = 0; j < s; j++)
    {
        __m256i sub = _mm256_sub_epi64(_mm256_sub_epi64(hi[j], bc->n[j]), borrow);
        borrow = _mm256_srli_epi64(sub, 63);
        d[j] = _mm256_and_si256(sub, mask);
    }
    take_d = _mm256_or_si256(_mm256_cmpeq_epi64(borrow, _mm256_setzero_si256()),
                             _mm256_cmpgt_epi64(overflow, _mm256_setzero_si256()));
    for (j = 0; j < s; j++)
        out[j] = _mm256_blendv_epi8(hi[j], d[j], take_d);
}

/*
   One CIOS iteration for the conversion path: t += a[]*xi, then cancel t[0] and
   shift the array down one limb (classic interleaved multiply-and-reduce).
*/
BM_FORCE_INLINE void
bm_cios_step(__m256i * t, const __m256i * a, __m256i xi, const bm_ctx_t * C,
             __m256i mask, const int s)
{
    __m256i carry, m;
    int j;

    carry = _mm256_setzero_si256();
    for (j = 0; j < s; j++)
    {
        __m256i p = _mm256_add_epi64(_mm256_add_epi64(t[j],
                          _mm256_mul_epu32(a[j], xi)), carry);
        t[j]  = _mm256_and_si256(p, mask);
        carry = _mm256_srli_epi64(p, BM_BITS);
    }
    {
        __m256i p = _mm256_add_epi64(t[s], carry);
        t[s]     = _mm256_and_si256(p, mask);
        t[s + 1] = _mm256_srli_epi64(p, BM_BITS);
    }

    m = _mm256_and_si256(_mm256_mul_epu32(t[0], C->n0inv), mask);
    carry = _mm256_srli_epi64(
                _mm256_add_epi64(t[0], _mm256_mul_epu32(m, C->n[0])), BM_BITS);
    for (j = 1; j < s; j++)
    {
        __m256i p = _mm256_add_epi64(_mm256_add_epi64(t[j],
                          _mm256_mul_epu32(m, C->n[j])), carry);
        t[j - 1] = _mm256_and_si256(p, mask);
        carry    = _mm256_srli_epi64(p, BM_BITS);
    }
    {
        __m256i p = _mm256_add_epi64(t[s], carry);
        t[s - 1] = _mm256_and_si256(p, mask);
        carry = _mm256_srli_epi64(p, BM_BITS);
        t[s] = _mm256_add_epi64(t[s + 1], carry);
    }
}

/*
   Reducing CIOS multiply r = a*x*R^-1 mod n, used only for the to/from-
   Montgomery conversions (runtime s; not the hot path).
*/
static inline void
bm_mul(__m256i * r, const __m256i * a, const __m256i * x, const bm_ctx_t * C)
{
    const __m256i mask = C->mask;
    const slong s = C->s;
    __m256i t[BM_MAXLIMB + 2];
    slong i;

    for (i = 0; i < s + 2; i++)
        t[i] = _mm256_setzero_si256();
    for (i = 0; i < s; i++)
        bm_cios_step(t, a, x[i], C, mask, s);

    bm_final_reduce(r, t, t[s], C, mask, s);
}

static inline void
bm_from(__m256i * r, const __m256i * a, const bm_ctx_t * C)
{
    __m256i one[BM_MAXLIMB];
    slong j;
    one[0] = _mm256_set1_epi64x(1);
    for (j = 1; j < C->s; j++)
        one[j] = _mm256_setzero_si256();
    bm_mul(r, a, one, C);
}

static inline void
bm_load(bm_mat Mm, const gr_mat_t A, gr_ctx_t ctx, const bm_ctx_t * bc)
{
    slong nlimbs = MPN_MOD_CTX_NLIMBS(ctx);
    const slong s = bc->s;
    ulong tmp[4];
    slong i, j, l;
    for (i = 0; i < BM_DIM; i++)
    {
        ulong pa[BM_DIM][BM_MAXLIMB];
        __m256i soa[BM_MAXLIMB], mo[BM_MAXLIMB];
        for (j = 0; j < BM_DIM; j++)
            bm_pack31(pa[j], (nn_srcptr) gr_mat_entry_srcptr(A, i, j, ctx), nlimbs, s);
        for (l = 0; l < s; l++)
            soa[l] = _mm256_setr_epi64x(pa[0][l], pa[1][l], pa[2][l], pa[3][l]);
        bm_mul(mo, soa, bc->r2, bc);
        for (l = 0; l < s; l++)
        {
            _mm256_storeu_si256((__m256i *) tmp, mo[l]);
            Mm[i][0][l] = tmp[0]; Mm[i][1][l] = tmp[1];
            Mm[i][2][l] = tmp[2]; Mm[i][3][l] = tmp[3];
        }
    }
}

static inline void
bm_store(gr_mat_t C, const bm_mat Mm, gr_ctx_t ctx, const bm_ctx_t * bc)
{
    slong nlimbs = MPN_MOD_CTX_NLIMBS(ctx);
    const slong s = bc->s;
    slong i, j, l;
    for (i = 0; i < BM_DIM; i++)
    {
        __m256i soa[BM_MAXLIMB], mo[BM_MAXLIMB];
        ulong lanes[BM_MAXLIMB][4];
        for (l = 0; l < s; l++)
            soa[l] = _mm256_setr_epi64x(Mm[i][0][l], Mm[i][1][l],
                                        Mm[i][2][l], Mm[i][3][l]);
        bm_from(mo, soa, bc);
        for (l = 0; l < s; l++)
            _mm256_storeu_si256((__m256i *) lanes[l], mo[l]);
        for (j = 0; j < BM_DIM; j++)
        {
            ulong out31[BM_MAXLIMB];
            for (l = 0; l < s; l++)
                out31[l] = lanes[l][j];
            bm_unpack31((nn_ptr) gr_mat_entry_ptr(C, i, j, ctx), out31, nlimbs, s);
        }
    }
}

/*
   Cm = Am * Bm, lazy reduction.  s is a COMPILE-TIME constant at each call site
   (see bm_matmul), so all loops unroll to a size-s specialization.
*/
BM_FORCE_INLINE void
bm_matmul_lazy(bm_mat Cm, const bm_mat Am, const bm_mat Bm, const bm_ctx_t * bc,
               const int s)
{
    const __m256i mask = bc->mask;
    __m256i Brow[BM_DIM][BM_MAXLIMB];
    ulong tmp[4];
    int i, d, bi, j;

    // Pack each row of B into SoA lanes: lane j carries output column j.
    for (d = 0; d < BM_DIM; d++)
        for (j = 0; j < s; j++)
            Brow[d][j] = _mm256_setr_epi64x(Bm[d][0][j], Bm[d][1][j],
                                            Bm[d][2][j], Bm[d][3][j]);

    for (i = 0; i < BM_DIM; i++)
    {
        __m256i T[2 * BM_MAXLIMB + 1];
        __m256i av[BM_MAXLIMB];
        __m256i out[BM_MAXLIMB];

        for (j = 0; j < 2 * s + 1; j++)
            T[j] = _mm256_setzero_si256();

        // Lazy: sum all four raw products A[i][d]*B[d] into the accumulator T
        // before reducing (one REDC per output entry instead of four).
        for (d = 0; d < BM_DIM; d++)
        {
            for (j = 0; j < s; j++)
                av[j] = _mm256_set1_epi64x(Am[i][d][j]);
            for (bi = 0; bi < s; bi++)
                bm_fma_column(T, av, Brow[d][bi], mask, bi, s);
        }

        // Reduce the summed 2s-limb accumulator once: normalize carries, run s
        // Montgomery reduction steps, then the final conditional subtract.
        bm_normalize(T, mask, 2 * s);
        for (bi = 0; bi < s; bi++)
            bm_redc_step(T, bi, bc, mask, s, 2 * s);
        bm_final_reduce(out, T + s, T[2 * s], bc, mask, s);

        for (j = 0; j < s; j++)
        {
            _mm256_storeu_si256((__m256i *) tmp, out[j]);
            Cm[i][0][j] = tmp[0]; Cm[i][1][j] = tmp[1];
            Cm[i][2][j] = tmp[2]; Cm[i][3][j] = tmp[3];
        }
    }
}

/* dispatch to the size-specialized (unrolled) instantiation */
static inline void
bm_matmul(bm_mat Cm, const bm_mat Am, const bm_mat Bm, const bm_ctx_t * bc)
{
    switch (bc->s)
    {
        case 2: bm_matmul_lazy(Cm, Am, Bm, bc, 2); break;
        case 3: bm_matmul_lazy(Cm, Am, Bm, bc, 3); break;
        case 4: bm_matmul_lazy(Cm, Am, Bm, bc, 4); break;
        case 5: bm_matmul_lazy(Cm, Am, Bm, bc, 5); break;
        case 6: bm_matmul_lazy(Cm, Am, Bm, bc, 6); break;
        case 7: bm_matmul_lazy(Cm, Am, Bm, bc, 7); break;
        default: bm_matmul_lazy(Cm, Am, Bm, bc, 8); break;
    }
}

#endif /* __AVX2__ */

#endif /* MPN_MOD_BATCHED_MONT_H */
