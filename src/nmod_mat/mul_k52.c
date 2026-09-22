/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Matrix multiplication over Z/nZ for n <= 2^52 with integer SIMD: a
    two-limb Karatsuba product on 32x32 -> 64 bit multiplications.

    nmod_mat_mul_u32 stops at 2^32 because its products must fit 64 bits.
    Here an entry a, lifted to the symmetric range |a| <= n/2 < 2^51, is
    split into two balanced limbs

        a = a1 * 2^27 + a0,    -2^26 <= a0 < 2^26,    |a1| <= 2^24,

    so that a*b = a1 b1 2^54 + (a0 b1 + a1 b0) 2^27 + a0 b0. The three
    coefficients are accumulated separately, the middle one as
    (a0 + a1)(b0 + b1) - a0 b0 - a1 b1 (Karatsuba), which makes three
    widening multiplications per product instead of four:

        acc0 += a0 * b0,    accM += (a0 + a1) * (b0 + b1),    acc2 += a1 * b1.

    Every partial product is below 2^52.7 in absolute value, so
    a k-block of KC = 256 products accumulates in signed 64-bit lanes with
    no reduction at all. The operands of a "k step" are one 64-bit lane per
    entry holding the two int32 limbs (a0 low, a1 high): on x86 vpmuldq reads
    the low halves, so the first limb is the lane itself, the second a shift,
    and their sum an add; on NEON the narrowing vmovn / vshrn do the same. One
    splat, one shift and one add per A entry, the same per B vector, and three
    multiply-adds per accumulator triple.

    At the end of a block the three accumulators are combined modulo n in
    double precision with the mulmod of fft_small (mul_fp_vec.h). Each
    accumulator x is first taken to (-0.51 n, 0.51 n): writing
    x = hi 2^32 + lo with hi signed, hi and lo convert exactly to doubles
    and x = mulmod(hi, 2^32 mod n) + lo followed by a reduce_pm1n, which
    is valid for any n <= 2^52 since |hi| < 2^31. Then

        acc0 + 2^27 accM' + 2^54 acc2  (accM' = accM - acc0 - acc2)

    is r0 + mulmod(rM, 2^27 mod n) + mulmod(r2, 2^54 mod n), with rM and r2
    brought to [-n/2, n/2] before the two mulmods and their results brought
    back to [-n/2, n/2] after (the mulmod inputs are then at most n^2/2 in
    absolute value, for which fft_small's analysis gives an output below
    1.26 n, see mul_fp_vec.h): the sum is below 1.51 n < 2^53 in absolute
    value, hence exact, and a last reduce_pm1n and a sign fix give the
    canonical residue, which is below 2^52 and stores as an integer.

    NOTE The floating point modular reduction is what limits the modulus to
    2^52: the integer part alone would accept 54 bits.

    On machines without AVX512-IFMA (where nmod_mat_mul_u52 should be used),
    this is a good candidate for 33 to 52 bit moduli. It does 3 multiply-adds
    plus overhead per product, against 4-5 dgemm passes and a CRT for
    nmod_mat_mul_blas, and 6-7 floating point operations for the
    all-floating-point kernel of mul_fp50.c. The profile
    p-mul_tune.c can be used to compare complete multiplications).

    The packing, the register-tile microkernel, the blocked core and the
    thread split are those of mul_blocked_templ.h; the backends (AVX-512,
    AVX2, NEON, plain C) follow the selection of mul_fp_vec.h.
*/

#include <string.h>
#include "nmod_mat.h"

#if FLINT_BITS == 64

#include "nmod.h"
#include "thread_pool.h"
#include "thread_support.h"
#include "nmod_mat/mul_fp_vec.h"

/* tile geometry and blocking ************************************************/

/*
    A tile holds 3*MR*NACC accumulators, 3*NACC B operands and the A
    operands of MR rows (3 on x86, 2 on NEON); the register budget picks
    MR and NACC (32 vectors on AVX-512 and NEON, 16 on AVX2). The defaults
    are the widest spill-free tiles, which on AVX2 and AVX-512 are also
    the best geometries measured by a standalone microkernel probe; on
    AVX-512 that probe preferred 4 x 2 (one vector over the register file,
    a few spills) and on NEON 4 x 2 as well, both by about 1.5%.
*/
#define K52_VL FPV_VL
#if defined(FPV_AVX512)
# ifndef K52_MR
#  define K52_MR 3          /* 18 + 6 + 3 = 27 zmm */
# endif
# ifndef K52_NACC
#  define K52_NACC 2
# endif
#elif defined(FPV_AVX2)
# ifndef K52_MR
#  define K52_MR 2          /* 6 + 3 + 3 = 12 ymm */
# endif
# ifndef K52_NACC
#  define K52_NACC 1
# endif
#elif defined(FPV_NEON)
# ifndef K52_MR
#  define K52_MR 3          /* 18 + 6 + 6 = 30 of 32 (4 x 2 spills) */
# endif
# ifndef K52_NACC
#  define K52_NACC 2
# endif
#else
# ifndef K52_MR
#  define K52_MR 2
# endif
# ifndef K52_NACC
#  define K52_NACC 2
# endif
#endif

#ifndef K52_KC
# define K52_KC 256
#endif
#ifndef K52_MC
# define K52_MC 96
#endif
#ifndef K52_NC
# define K52_NC 2048
#endif
/* products per thread, as U32_MT_MIN_WORK in mul_u32.c */
#ifndef K52_MT_MIN_WORK
# define K52_MT_MIN_WORK 500000.0
#endif

/* limb size: a0 in [-2^26, 2^26), |a1| <= 2^24 for |a| <= 2^51 */
#define K52_L 27

/* the middle term accM - acc0 - acc2 of a block, below
   KC * (2^52.7 + 2^52 + 2^48.1) + 2^52, must fit 63 bits */
#if K52_KC > 512
# error "K52_KC too large for the 64-bit accumulators"
#endif

/* modulus-derived parameters ************************************************/

typedef struct
{
    ulong n;
    double nd;
    double ninv;
    double c32;     /* 2^32 mod n */
    double cL;      /* 2^L mod n */
    double c2L;     /* 2^(2L) mod n */
}
k52_ctx_struct;

static void
k52_ctx_init(k52_ctx_struct * ctx, ulong n)
{
    FLINT_ASSERT(n >= 1 && n <= (UWORD(1) << 52));

    ctx->n = n;
    ctx->nd = (double) n;
    ctx->ninv = 1.0 / (double) n;
    ctx->c32 = (double) ((UWORD(1) << 32) % n);
    ctx->cL = (double) ((UWORD(1) << K52_L) % n);
    ctx->c2L = (double) ((UWORD(1) << (2 * K52_L)) % n);
}

/* symmetric lift then balanced split; the two int32 limbs share a word,
   a0 in the low half */
FLINT_FORCE_INLINE ulong
k52_lift(ulong a, ulong n)
{
    slong s = (a > n / 2) ? (slong) a - (slong) n : (slong) a;
    slong a0 = ((s + (WORD(1) << (K52_L - 1))) & ((WORD(1) << K52_L) - 1))
                   - (WORD(1) << (K52_L - 1));
    slong a1 = (s - a0) >> K52_L;

    return (ulong) (uint32_t) a0 | ((ulong) (uint32_t) a1 << 32);
}

/* primitives ****************************************************************/

/*
    Each backend provides

      k52_vi           vector of K52_VL signed 64-bit lanes (fpv_i)
      k52_acc          the accumulator triple acc0, accM, acc2
      k52_bv           the three B operands b0, b1, b0 + b1 of a column group
      k52_av           the A operands of an entry
      k52_load_bstep   the K52_NACC B operand triples of a k step
      k52_load_a       the A operands of one packed entry
      k52_acc_mul_add  the three widening multiply-adds of a product
      k52_sub          lane-wise difference (for accM - acc0 - acc2)

    and shares the rest (zero, load of C, reduction) below.
*/

#if defined(FPV_AVX512)

/* AVX-512 (F + DQ) *********************************************************/

typedef __m512i k52_vi;
typedef struct { k52_vi a0, aM, a2; } k52_acc;
typedef struct { __m512i b0, b1, bs; } k52_bv;
typedef struct { __m512i a0, a1, as; } k52_av;

FLINT_FORCE_INLINE k52_vi k52_vi_zero(void) { return _mm512_setzero_si512(); }
FLINT_FORCE_INLINE k52_vi k52_load_c(const ulong * p) { return _mm512_loadu_si512((const void *) p); }
FLINT_FORCE_INLINE k52_vi k52_sub(k52_vi a, k52_vi b) { return _mm512_sub_epi64(a, b); }

/* vpmuldq reads the low 32 bits of each lane: the a0 limb as loaded, the
   a1 limb after a shift, and the sum of the two in the low half */
FLINT_FORCE_INLINE void
k52_load_bstep(k52_bv * bv, const ulong * p)
{
    slong v;

    for (v = 0; v < K52_NACC; v++)
    {
        __m512i t = _mm512_load_si512((const void *) (p + v * K52_VL));
        bv[v].b0 = t;
        bv[v].b1 = _mm512_srli_epi64(t, 32);
        bv[v].bs = _mm512_add_epi32(bv[v].b0, bv[v].b1);
    }
}

FLINT_FORCE_INLINE k52_av
k52_load_a(const ulong * p)
{
    k52_av a;
    a.a0 = _mm512_set1_epi64((long long) *p);
    a.a1 = _mm512_srli_epi64(a.a0, 32);
    a.as = _mm512_add_epi32(a.a0, a.a1);
    return a;
}

FLINT_FORCE_INLINE k52_acc
k52_acc_mul_add(k52_acc acc, k52_av a, k52_bv b)
{
    acc.a0 = _mm512_add_epi64(acc.a0, _mm512_mul_epi32(a.a0, b.b0));
    acc.aM = _mm512_add_epi64(acc.aM, _mm512_mul_epi32(a.as, b.bs));
    acc.a2 = _mm512_add_epi64(acc.a2, _mm512_mul_epi32(a.a1, b.b1));
    return acc;
}

#elif defined(FPV_AVX2)

/* AVX2 *********************************************************************/

typedef __m256i k52_vi;
typedef struct { k52_vi a0, aM, a2; } k52_acc;
typedef struct { __m256i b0, b1, bs; } k52_bv;
typedef struct { __m256i a0, a1, as; } k52_av;

FLINT_FORCE_INLINE k52_vi k52_vi_zero(void) { return _mm256_setzero_si256(); }
FLINT_FORCE_INLINE k52_vi k52_load_c(const ulong * p) { return _mm256_loadu_si256((const __m256i *) p); }
FLINT_FORCE_INLINE k52_vi k52_sub(k52_vi a, k52_vi b) { return _mm256_sub_epi64(a, b); }

FLINT_FORCE_INLINE void
k52_load_bstep(k52_bv * bv, const ulong * p)
{
    slong v;

    for (v = 0; v < K52_NACC; v++)
    {
        __m256i t = _mm256_load_si256((const __m256i *) (p + v * K52_VL));
        bv[v].b0 = t;
        bv[v].b1 = _mm256_srli_epi64(t, 32);
        bv[v].bs = _mm256_add_epi32(bv[v].b0, bv[v].b1);
    }
}

FLINT_FORCE_INLINE k52_av
k52_load_a(const ulong * p)
{
    k52_av a;
    a.a0 = _mm256_set1_epi64x((long long) *p);
    a.a1 = _mm256_srli_epi64(a.a0, 32);
    a.as = _mm256_add_epi32(a.a0, a.a1);
    return a;
}

FLINT_FORCE_INLINE k52_acc
k52_acc_mul_add(k52_acc acc, k52_av a, k52_bv b)
{
    acc.a0 = _mm256_add_epi64(acc.a0, _mm256_mul_epi32(a.a0, b.b0));
    acc.aM = _mm256_add_epi64(acc.aM, _mm256_mul_epi32(a.as, b.bs));
    acc.a2 = _mm256_add_epi64(acc.a2, _mm256_mul_epi32(a.a1, b.b1));
    return acc;
}

#elif defined(FPV_NEON)

/* AArch64 NEON *************************************************************/

/*
    smlal takes narrow int32x2 operands, or one such operand and a lane of
    another: the limbs of B are narrowed from the 64-bit lanes, and the two
    limbs of an A entry are the two lanes of one 64-bit register, their sum
    a pairwise add of it.
*/
typedef int64x2_t k52_vi;
typedef struct { k52_vi a0, aM, a2; } k52_acc;
typedef struct { int32x2_t b0, b1, bs; } k52_bv;
typedef struct { int32x2_t a01, as; } k52_av;

FLINT_FORCE_INLINE k52_vi k52_vi_zero(void) { return vdupq_n_s64(0); }
FLINT_FORCE_INLINE k52_vi k52_load_c(const ulong * p) { return vreinterpretq_s64_u64(vld1q_u64((const uint64_t *) p)); }
FLINT_FORCE_INLINE k52_vi k52_sub(k52_vi a, k52_vi b) { return vsubq_s64(a, b); }

FLINT_FORCE_INLINE void
k52_load_bstep(k52_bv * bv, const ulong * p)
{
    slong v;

    for (v = 0; v < K52_NACC; v++)
    {
        int64x2_t t = vreinterpretq_s64_u64(vld1q_u64((const uint64_t *) (p + v * K52_VL)));
        bv[v].b0 = vmovn_s64(t);
        bv[v].b1 = vshrn_n_s64(t, 32);
        bv[v].bs = vadd_s32(bv[v].b0, bv[v].b1);
    }
}

FLINT_FORCE_INLINE k52_av
k52_load_a(const ulong * p)
{
    k52_av a;
    a.a01 = vcreate_s32(*p);            /* lane 0 = a0, lane 1 = a1 */
    a.as = vpadd_s32(a.a01, a.a01);
    return a;
}

FLINT_FORCE_INLINE k52_acc
k52_acc_mul_add(k52_acc acc, k52_av a, k52_bv b)
{
    acc.a0 = vmlal_lane_s32(acc.a0, b.b0, a.a01, 0);
    acc.aM = vmlal_lane_s32(acc.aM, b.bs, a.as, 0);
    acc.a2 = vmlal_lane_s32(acc.a2, b.b1, a.a01, 1);
    return acc;
}

#else

/* plain C ******************************************************************/

typedef fpv_i k52_vi;
typedef struct { k52_vi a0, aM, a2; } k52_acc;
typedef struct { slong b0[K52_VL], b1[K52_VL], bs[K52_VL]; } k52_bv;
typedef struct { slong a0, a1, as; } k52_av;

FLINT_FORCE_INLINE k52_vi
k52_vi_zero(void)
{
    k52_vi a;
    slong i;
    for (i = 0; i < K52_VL; i++)
        a.v[i] = 0;
    return a;
}

FLINT_FORCE_INLINE k52_vi
k52_load_c(const ulong * p)
{
    k52_vi a;
    slong i;
    for (i = 0; i < K52_VL; i++)
        a.v[i] = (slong) p[i];
    return a;
}

FLINT_FORCE_INLINE k52_vi
k52_sub(k52_vi a, k52_vi b)
{
    slong i;
    for (i = 0; i < K52_VL; i++)
        a.v[i] -= b.v[i];
    return a;
}

FLINT_FORCE_INLINE void
k52_load_bstep(k52_bv * bv, const ulong * p)
{
    slong v, i;

    for (v = 0; v < K52_NACC; v++)
        for (i = 0; i < K52_VL; i++)
        {
            ulong t = p[v * K52_VL + i];
            bv[v].b0[i] = (int32_t) (uint32_t) t;
            bv[v].b1[i] = (int32_t) (uint32_t) (t >> 32);
            bv[v].bs[i] = bv[v].b0[i] + bv[v].b1[i];
        }
}

FLINT_FORCE_INLINE k52_av
k52_load_a(const ulong * p)
{
    k52_av a;
    a.a0 = (int32_t) (uint32_t) *p;
    a.a1 = (int32_t) (uint32_t) (*p >> 32);
    a.as = a.a0 + a.a1;
    return a;
}

FLINT_FORCE_INLINE k52_acc
k52_acc_mul_add(k52_acc acc, k52_av a, k52_bv b)
{
    slong i;
    for (i = 0; i < K52_VL; i++)
    {
        acc.a0.v[i] += a.a0 * b.b0[i];
        acc.aM.v[i] += a.as * b.bs[i];
        acc.a2.v[i] += a.a1 * b.b1[i];
    }
    return acc;
}

#endif

/* shared across backends ****************************************************/

FLINT_FORCE_INLINE k52_acc
k52_acc_zero(void)
{
    k52_acc a;
    a.a0 = k52_vi_zero();
    a.aM = k52_vi_zero();
    a.a2 = k52_vi_zero();
    return a;
}

/* the entry of C, canonical, goes to the weight-1 column; it goes to the
   middle column too, so that it cancels in accM - acc0 - acc2 */
FLINT_FORCE_INLINE k52_acc
k52_acc_load_c(const ulong * p)
{
    k52_acc a;
    a.a0 = k52_load_c(p);
    a.aM = a.a0;
    a.a2 = k52_vi_zero();
    return a;
}

typedef struct
{
    fpv nv;
    fpv ninvv;
    fpv c32v;
    fpv cLv;
    fpv c2Lv;
}
k52_consts;

FLINT_FORCE_INLINE k52_consts
k52_consts_init(const k52_ctx_struct * ctx)
{
    k52_consts C;

    C.nv = fpv_set1(ctx->nd);
    C.ninvv = fpv_set1(ctx->ninv);
    C.c32v = fpv_set1(ctx->c32);
    C.cLv = fpv_set1(ctx->cL);
    C.c2Lv = fpv_set1(ctx->c2L);

    return C;
}

/* signed 64-bit lanes x -> x mod n in (-0.51 n, 0.51 n), as doubles */
FLINT_FORCE_INLINE fpv
k52_red(k52_vi x, const k52_consts * C)
{
    fpv hd, ld, t;

    fpv_split_i64(x, &hd, &ld);
    t = fpv_add(fpv_mulmod(hd, C->c32v, C->nv, C->ninvv), ld);

    return fpv_reduce_pm1n(t, C->nv, C->ninvv);
}

/*
    acc0 + 2^L (accM - acc0 - acc2) + 2^(2L) acc2 mod n, canonical.

    Kept out of line: inlined, the many vector constants of the reduction
    stay live across the k loop, and on AVX2 (16 registers) that can make the
    compiler spill accumulators of the microkernel; the call runs once per
    accumulator triple per k-block, which is negligible.
*/
FLINT_STATIC_NOINLINE fpv
k52_finish(k52_acc acc, const k52_consts * C)
{
    fpv r0, rM, r2, v;

    r0 = k52_red(acc.a0, C);
    rM = k52_red(k52_sub(k52_sub(acc.aM, acc.a0), acc.a2), C);
    r2 = k52_red(acc.a2, C);

    /* to [-n/2, n/2] first: the mulmod then rounds a quotient of at most
       (n-1)/2 < 2^51, which is what fpv_rint_mul wants (mul_fp_vec.h) */
    rM = fpv_pm1n_to_pmhn(rM, C->nv);
    r2 = fpv_pm1n_to_pmhn(r2, C->nv);

    rM = fpv_pm1n_to_pmhn(fpv_mulmod(rM, C->cLv, C->nv, C->ninvv), C->nv);
    r2 = fpv_pm1n_to_pmhn(fpv_mulmod(r2, C->c2Lv, C->nv, C->ninvv), C->nv);

    v = fpv_add(fpv_add(r0, rM), r2);

    return fpv_reduce_0n(fpv_reduce_pm1n(v, C->nv, C->ninvv), C->nv);
}

/* template instantiation ****************************************************/

#define BT_NAME(x) k52_##x
#define BT_ENTRY ulong
#define BT_PACKED ulong
#define BT_CTX k52_ctx_struct
#define BT_MR K52_MR
#define BT_VL K52_VL
#define BT_NACC K52_NACC
#define BT_KC K52_KC
#define BT_MC K52_MC
#define BT_NC K52_NC
#define BT_MT_MIN_WORK K52_MT_MIN_WORK
#define BT_LIFT(x, ctx) k52_lift((x), (ctx)->n)
#define BT_BSLOT(j) (j)
#define BT_ACC k52_acc
#define BT_BV k52_bv
#define BT_AV k52_av
#define BT_CONSTS k52_consts
#define BT_CONSTS_INIT(ctx) k52_consts_init(ctx)
#define BT_ACC_ZERO() k52_acc_zero()
#define BT_LOAD_C(p) k52_acc_load_c(p)
#define BT_STORE_C(p, acc) fpv_store_u64(p, acc)
#define BT_LOAD_BSTEP(bv, p) k52_load_bstep(bv, p)
#define BT_LOAD_A(p) k52_load_a(p)
#define BT_MUL_ADD(acc, a, b, C) k52_acc_mul_add(acc, a, b)
/* no in-loop reduction: the accumulators hold a whole k-block */
#define BT_CADENCE(ctx) K52_KC
#define BT_FOLD(acc, C, ctx) (acc)
#define BT_FINISH(acc, C, ctx) k52_finish(acc, C)
#include "mul_blocked_templ.h"

/* public entry *************************************************************/

int
nmod_mat_mul_k52(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    slong m = A->r;
    slong k = A->c;
    slong n = B->c;
    ulong modn = C->mod.n;
    k52_ctx_struct ctx;

    FLINT_ASSERT(C->r == A->r);
    FLINT_ASSERT(C->c == B->c);
    FLINT_ASSERT(A->c == B->r);

    if (modn > (UWORD(1) << 52))
        return 0;

    if (m <= 0 || n <= 0)
        return 1;

    if (k <= 0 || modn == 1)
    {
        nmod_mat_zero(C);
        return 1;
    }

    if (C == A || C == B)
    {
        nmod_mat_t T;
        nmod_mat_init(T, m, n, modn);
        nmod_mat_mul_k52(T, A, B);
        nmod_mat_swap_entrywise(C, T);
        nmod_mat_clear(T);
        return 1;
    }

    k52_ctx_init(&ctx, modn);

    k52_core_mt(C->entries, C->stride, A->entries, A->stride,
                B->entries, B->stride, m, k, n, &ctx,
                flint_get_num_threads());

    return 1;
}

#else

/* the packing and the reduction here assume 64-bit words */
int
nmod_mat_mul_k52(nmod_mat_t FLINT_UNUSED(C),
                 const nmod_mat_t FLINT_UNUSED(A),
                 const nmod_mat_t FLINT_UNUSED(B))
{
    return 0;
}

#endif
