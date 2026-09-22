/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Matrix multiplication over Z/nZ for n <= 2^52 with AVX512-IFMA.

    The IFMA instructions vpmadd52luq / vpmadd52huq multiply the low 52 bits
    of two unsigned 64-bit lanes and add the low, respectively the high, 52
    bits of the 104-bit product to a 64-bit accumulator. For canonical
    entries a, b < n <= 2^52 the pair (lo, hi) with ab = hi*2^52 + lo is
    therefore exact, and the two accumulators have 12 bits of headroom each:
    up to 4095 products can be summed with no reduction at all, so a k-block
    of KC <= 2048 needs none. Compared with nmod_mat_mul_u32 this covers
    moduli up to 2^52 in one pass (where nmod_mat_mul_blas needs 4-5 dgemm
    passes and a CRT, and nmod_mat_mul_u32 stops at 2^32), at a cost of two
    IFMA per product-vector against one vpmuldq plus one vpaddq, i.e. at the
    same number of vector operations. When (n-1)^2 < 2^52, that is n <= 2^26,
    the high halves are all zero, and a second microkernel keeps only the
    low accumulator: one instruction per product-vector, twice the rows per
    tile. Measured against nmod_mat_mul_u32 (single thread, 32 to 2048):
    the single-IFMA mode is 1.3-1.35x faster on Ice Lake (two IFMA ports,
    where u32 is issue bound) and 1.1-1.2x on Zen 4 (one 512-bit multiply
    per cycle, where u32 still pays for its separate add and shift); the
    two-IFMA mode is 1.2-1.65x slower than u32 up to 30 bits, on par at 31
    bits on Ice Lake, and 1.3-2.7x faster from 32 bits on where u32 folds
    every 1-2 products. Hence the U52_* parameters of flint-mparam.h, see
    the dispatch in mul.c.

    The end-of-block reduction takes the two accumulators lo, hi < 2^61 to
    the canonical residue of hi*2^52 + lo:

      h  = hi mod n                     (Barrett through double precision)
      r  = h*c - q*n, c = 2^52 mod n,   (Shoup: q = floor(h*w / 2^52) with
           r in [0, 2n)                  w = floor(c*2^52 / n); the two
                                         104-bit products come from IFMA)
      C  = (r + lo) mod n               (Barrett again).

    The Barrett quotients are off by less than one: the argument t is
    either exact in double precision (t <= 3n + 256 n^2 is below 2^53 when
    n < 2^9) or below 2^61 with n >= 2^9, where the absolute error of its
    double approximation, at most 2^8, is below n/2; with the rounding to
    the nearest quotient the remainder lies in [-n, n] and one correction on
    each side suffices. This runs once per output entry per k-block.

    The packing, the register-tile microkernel, the blocked core and the
    thread split are those of mul_blocked_templ.h, shared with mul_u32.c.
    Panels hold the canonical entries as 64-bit lanes (IFMA reads 64-bit
    lanes), so no lift is needed, and B panels keep their columns in
    order.
*/

#include <string.h>
#include "nmod_mat.h"
#include "nmod_mat/impl.h"

#if NMOD_MAT_HAVE_MUL_U52

#include "nmod.h"
#include "longlong.h"
#include "machine_vectors.h"
#include "thread_pool.h"
#include "thread_support.h"

/* tile geometry and blocking ************************************************/

#define U52_VL 8
#define U52_NACC 2
/* lo+hi: 2*MR*NACC accumulators + NACC B operands + 1 A broadcast <= 32 */
#ifndef U52_MR_HI
# define U52_MR_HI 7
#endif
/* lo only: MR*NACC accumulators + NACC + 1 <= 32 */
#ifndef U52_MR_LO
# define U52_MR_LO 14
#endif

#ifndef U52_KC
# define U52_KC 256
#endif
#ifndef U52_MC
# define U52_MC 96
#endif
#ifndef U52_NC
# define U52_NC 2048
#endif
/* products per thread, as U32_MT_MIN_WORK in mul_u32.c */
#ifndef U52_MT_MIN_WORK
# define U52_MT_MIN_WORK 500000.0
#endif

/* KC + 1 products of 52 bits (the +1 for the entry of C) must fit 64 bits */
#if U52_KC > 2048
# error "U52_KC too large for the 64-bit accumulators"
#endif

/* modulus-derived parameters ************************************************/

typedef struct
{
    ulong n;
    ulong c52;      /* 2^52 mod n */
    ulong w;        /* floor(c52 * 2^52 / n), for the Shoup step */
    double ninv;
    int lo_only;    /* n <= 2^26: all high halves vanish */
}
u52_ctx_struct;

static void
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

/* primitives ****************************************************************/

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

/* canonical residue of x for 0 <= x < 2^61, x <= 3n + 256 n^2 */
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

/* the two accumulator flavours */
typedef struct { __m512i lo, hi; } u52_acc2;

FLINT_FORCE_INLINE u52_acc2
u52_acc2_zero(void)
{
    u52_acc2 a;
    a.lo = _mm512_setzero_si512();
    a.hi = _mm512_setzero_si512();
    return a;
}

FLINT_FORCE_INLINE u52_acc2
u52_acc2_load_c(const ulong * p)
{
    u52_acc2 a;
    a.lo = _mm512_loadu_si512((const void *) p);
    a.hi = _mm512_setzero_si512();
    return a;
}

FLINT_FORCE_INLINE u52_acc2
u52_acc2_mul_add(u52_acc2 acc, __m512i a, __m512i b)
{
    acc.lo = u52_madd52lo(acc.lo, a, b);
    acc.hi = u52_madd52hi(acc.hi, a, b);
    return acc;
}

FLINT_FORCE_INLINE __m512i
u52_acc2_finish(u52_acc2 acc, const u52_consts * C)
{
    return u52_red2(acc.lo, acc.hi, C);
}

FLINT_FORCE_INLINE __m512i
u52_acc1_load_c(const ulong * p) { return _mm512_loadu_si512((const void *) p); }

FLINT_FORCE_INLINE __m512i
u52_acc1_mul_add(__m512i acc, __m512i a, __m512i b)
{
    return u52_madd52lo(acc, a, b);
}

FLINT_FORCE_INLINE void
u52_store_c(ulong * p, __m512i a) { _mm512_storeu_si512((void *) p, a); }

FLINT_FORCE_INLINE void
u52_load_bstep(__m512i * bv, const ulong * p)
{
    bv[0] = _mm512_load_si512((const void *) p);
    bv[1] = _mm512_load_si512((const void *) (p + U52_VL));
}

FLINT_FORCE_INLINE __m512i
u52_load_a(const ulong * p) { return _mm512_set1_epi64(*p); }

/* template instantiations ***************************************************/

#define BT_ENTRY ulong
#define BT_PACKED ulong
#define BT_CTX u52_ctx_struct
#define BT_VL U52_VL
#define BT_NACC U52_NACC
#define BT_KC U52_KC
#define BT_MC U52_MC
#define BT_NC U52_NC
#define BT_MT_MIN_WORK U52_MT_MIN_WORK
#define BT_LIFT(x, ctx) (x)
#define BT_BSLOT(j) (j)
#define BT_BV __m512i
#define BT_AV __m512i
#define BT_CONSTS u52_consts
#define BT_CONSTS_INIT(ctx) u52_consts_init(ctx)
#define BT_STORE_C(p, acc) u52_store_c(p, acc)
#define BT_LOAD_BSTEP(bv, p) u52_load_bstep(bv, p)
#define BT_LOAD_A(p) u52_load_a(p)
/* no in-loop reduction: the accumulators hold a whole k-block */
#define BT_CADENCE(ctx) U52_KC
#define BT_FOLD(acc, C, ctx) (acc)

/* lo + hi accumulators, any n <= 2^52 */
#define BT_NAME(x) u52_hi_##x
#define BT_MR U52_MR_HI
#define BT_ACC u52_acc2
#define BT_ACC_ZERO() u52_acc2_zero()
#define BT_LOAD_C(p) u52_acc2_load_c(p)
#define BT_MUL_ADD(acc, a, b, C) u52_acc2_mul_add(acc, a, b)
#define BT_FINISH(acc, C, ctx) u52_acc2_finish(acc, C)
#include "mul_blocked_templ.h"
#undef BT_NAME
#undef BT_MR
#undef BT_ACC
#undef BT_ACC_ZERO
#undef BT_LOAD_C
#undef BT_MUL_ADD
#undef BT_FINISH

/* lo accumulator only, n <= 2^26 */
#define BT_NAME(x) u52_lo_##x
#define BT_MR U52_MR_LO
#define BT_ACC __m512i
#define BT_ACC_ZERO() _mm512_setzero_si512()
#define BT_LOAD_C(p) u52_acc1_load_c(p)
#define BT_MUL_ADD(acc, a, b, C) u52_acc1_mul_add(acc, a, b)
#define BT_FINISH(acc, C, ctx) u52_red(acc, C)
#include "mul_blocked_templ.h"
#undef BT_NAME
#undef BT_MR
#undef BT_ACC
#undef BT_ACC_ZERO
#undef BT_LOAD_C
#undef BT_MUL_ADD
#undef BT_FINISH

/* public entry *************************************************************/

int
nmod_mat_mul_u52(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    slong m = A->r;
    slong k = A->c;
    slong n = B->c;
    ulong modn = C->mod.n;
    u52_ctx_struct ctx;

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
        nmod_mat_mul_u52(T, A, B);
        nmod_mat_swap_entrywise(C, T);
        nmod_mat_clear(T);
        return 1;
    }

    u52_ctx_init(&ctx, modn);

    if (ctx.lo_only)
        u52_lo_core_mt(C->entries, C->stride, A->entries, A->stride,
                       B->entries, B->stride, m, k, n, &ctx,
                       flint_get_num_threads());
    else
        u52_hi_core_mt(C->entries, C->stride, A->entries, A->stride,
                       B->entries, B->stride, m, k, n, &ctx,
                       flint_get_num_threads());

    return 1;
}

#else

/* no AVX512-IFMA at compile time (or a 32-bit build) */
int
nmod_mat_mul_u52(nmod_mat_t FLINT_UNUSED(C),
                 const nmod_mat_t FLINT_UNUSED(A),
                 const nmod_mat_t FLINT_UNUSED(B))
{
    return 0;
}

#endif
