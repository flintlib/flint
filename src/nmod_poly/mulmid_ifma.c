/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "impl.h"
#if FLINT_HAVE_FFT_SMALL
# include "fft_small.h"
# include "fft_small/impl.h"
#endif

int
_nmod_poly_is_fft_small_prime(nmod_t mod)
{
#if FLINT_HAVE_FFT_SMALL
    return _fft_small_mpn_ctx_prime_index(get_default_mpn_ctx(), mod.n) >= 0;
#else
    (void) mod;
    return 0;
#endif
}

/*
    Coefficients nlo to nhi - 1 of the product, for n <= 2^52 and a shorter
    operand of at most 511 coefficients (see impl.h). Assumes len1, len2 >= 1
    and 0 <= nlo < nhi <= len1 + len2 - 1, does not allow aliasing.

    With a the shorter operand, of length la, and b the longer one, of
    length lb, c[k] = sum_i a[i] b[k - i]. The outputs are computed by
    blocks of R vectors of 8 coefficients (W = 8R outputs from k0), whose
    accumulators stay in registers for all their terms. The terms are taken
    in diagonal order: for each offset d, one load of the window
    b[k0 + d], ..., b[k0 + d + 7] serves the R vectors, vector r with
    a[8r - d] broadcast. The range of d is cut into runs on which the vectors
    with 0 <= 8r - d < la form a fixed interval [rlo, rhi], each run going
    to unrolled code: no work outside the product apart from the lanes of
    the vectors at both ends. The blocks whose windows leave b read them
    from zero-padded copies of its ends.

    The products of 52-bit entries are accumulated as their low and high
    52-bit halves (vpmadd52luq, vpmadd52huq); for n <= 2^26 they are below
    2^52 and the low halves suffice. With at most 511 terms, every
    accumulator stays below 511 * 2^52 < 2^61. The sum V = lo + 2^52 hi
    of each output is reduced, for odd n > 2^26, by two steps of REDC in
    radix 2^52, a having been multiplied by 2^104 mod n (1.20-1.30x faster
    than u52_red2; with the low halves alone, u52_red is faster than one
    step of REDC with the multiplication of a); for n a power of 2, by a
    mask; otherwise by u52_red2 (u52_red for V = lo) of
    nmod_mat/mul_u52_vec.h, valid for lo, hi <= 3n + 511 n^2.
*/

#if NMOD_POLY_HAVE_IFMA

#include <immintrin.h>
#include "nmod_mat/mul_u52_vec.h"

#define IFMA_LO(acc, x, y) _mm512_madd52lo_epu64(acc, x, y)
#define IFMA_HI(acc, x, y) _mm512_madd52hi_epu64(acc, x, y)

/*
    REDC in radix 2^52, n odd: a step maps X = x0 + 2^52 x1 (x0 < 2^52) to
    (X + m n) / 2^52 with m = x0 n' mod 2^52 and n' = -1/n mod 2^52, that is
    to x1 + hi52(m n) + [x0 != 0], since x0 + lo52(m n) is 0 or 2^52.
    For X = lo + 2^52 hi, lo, hi < 2^61, a first step gives Y < 2^62, and a
    second one a result below 2^10 + n <= 2n (n >= 2^10).
*/
typedef struct
{
    __m512i n;
    __m512i nprime;     /* -1/n mod 2^52 */
    __m512i mask52;
    __m512i one;
}
redc52_consts;

static void
redc52_consts_init(redc52_consts * C, ulong n)
{
    C->n = _mm512_set1_epi64(n);
    C->nprime = _mm512_set1_epi64((-n_binvert(n)) & ((UWORD(1) << 52) - 1));
    C->mask52 = _mm512_set1_epi64((UWORD(1) << 52) - 1);
    C->one = _mm512_set1_epi64(1);
}

FLINT_FORCE_INLINE __m512i
redc52_step(__m512i x, __m512i x1, const redc52_consts * C)
{
    const __m512i zero = _mm512_setzero_si512();
    __m512i x0 = _mm512_and_si512(x, C->mask52);
    __m512i m = IFMA_LO(zero, x0, C->nprime);
    __m512i r = IFMA_HI(x1, m, C->n);

    return _mm512_mask_add_epi64(r, _mm512_test_epi64_mask(x0, x0), r, C->one);
}

/* (lo + 2^52 hi) 2^-104 mod n, canonical */
FLINT_FORCE_INLINE __m512i
redc52_2(__m512i lo, __m512i hi, const redc52_consts * C)
{
    __m512i r = redc52_step(lo, _mm512_add_epi64(hi, _mm512_srli_epi64(lo, 52)), C);

    r = redc52_step(r, _mm512_srli_epi64(r, 52), C);

    return _mm512_mask_sub_epi64(r, _mm512_cmpge_epu64_mask(r, C->n), r, C->n);
}

typedef struct
{
    nn_srcptr b;        /* the longer operand, length lb */
    slong lb;
    nn_srcptr head;     /* head[8 + j] = b[j] for j < hend, zero around */
    slong hend;
    nn_srcptr tail;     /* tail[j - tbeg] = b[j] for tbeg <= j < lb, zero after */
    slong tbeg;
    int redc;           /* 1: REDC (lo / hi), 2: n a power of 2, 0: u52_red */
    __m512i mask;       /* n - 1 for n a power of 2 */
    redc52_consts rc;
    u52_consts uc;
}
ifma_ctx_struct;

/* bk with bk[d] = b[k0 + d] for da <= d < dend + 7: b itself where the
   windows stay inside, otherwise one of its copies */
FLINT_FORCE_INLINE const ulong *
ifma_bk(const ifma_ctx_struct * C, slong k0, slong da, slong dend)
{
    if (k0 + da >= 0 && k0 + dend + 6 < C->lb)
        return C->b + k0;
    else if (k0 + da < 0 || C->tail == NULL)
        return C->head + 8 + k0;
    else
        return C->tail + (k0 - C->tbeg);
}

/* the switch on the run of vectors rlo to rhi of a block, IFMA_RUN being
   defined by the block */
#define IFMA_CASE(RLO, RHI) \
    case (RLO) * 8 + (RHI): IFMA_RUN(RLO, RHI); break;

#define IFMA_CASES_R6 \
    IFMA_CASE(0, 0) IFMA_CASE(0, 1) IFMA_CASE(0, 2) IFMA_CASE(0, 3) \
    IFMA_CASE(0, 4) IFMA_CASE(0, 5) IFMA_CASE(1, 1) IFMA_CASE(1, 2) \
    IFMA_CASE(1, 3) IFMA_CASE(1, 4) IFMA_CASE(1, 5) IFMA_CASE(2, 2) \
    IFMA_CASE(2, 3) IFMA_CASE(2, 4) IFMA_CASE(2, 5) IFMA_CASE(3, 3) \
    IFMA_CASE(3, 4) IFMA_CASE(3, 5) IFMA_CASE(4, 4) IFMA_CASE(4, 5) \
    IFMA_CASE(5, 5)

#define IFMA_CASES_R8 IFMA_CASES_R6 \
    IFMA_CASE(0, 6) IFMA_CASE(1, 6) IFMA_CASE(2, 6) IFMA_CASE(3, 6) \
    IFMA_CASE(4, 6) IFMA_CASE(5, 6) IFMA_CASE(6, 6) IFMA_CASE(0, 7) \
    IFMA_CASE(1, 7) IFMA_CASE(2, 7) IFMA_CASE(3, 7) IFMA_CASE(4, 7) \
    IFMA_CASE(5, 7) IFMA_CASE(6, 7) IFMA_CASE(7, 7)

#define IFMA_SWITCH(R, key) \
    if ((R) == 8) { switch (key) { IFMA_CASES_R8 default: break; } } \
    else { switch (key) { IFMA_CASES_R6 default: break; } }

/* the offsets of a block, from da to dend - 1: those with a term,
   0 <= 8r - d < la for some r <= rmax, 0 <= k0 + d + l < lb for some lane
   l; then, by runs, rlo, rhi and the end db of the run */
#define IFMA_BLOCK_RANGE \
    rmax = FLINT_MIN(R, (nout + 7) / 8) - 1; \
    da = FLINT_MAX(-k0 - 7, -(la - 1)); \
    dend = FLINT_MIN(C->lb - 1 - k0, 8 * rmax) + 1;

#define IFMA_RUN_RANGE \
    rlo = (da <= 0) ? 0 : (da + 7) / 8; \
    rhi = FLINT_MIN((la - 1 + da) / 8, rmax); \
    /* the run ends where rlo or rhi grows */ \
    db = FLINT_MIN(dend, 8 * rlo + 1); \
    if (rhi < rmax) \
    { \
        t = 8 * (rhi + 1) - (la - 1); \
        db = FLINT_MIN(db, t); \
    }

FLINT_FORCE_INLINE __mmask8
ifma_store_mask(slong nout, int r)
{
    return (8 * r + 8 <= nout) ? 0xff : (__mmask8) ((UWORD(1) << (nout - 8 * r)) - 1);
}

/* the vectors rlo to rhi, offsets d from da to db - 1 */
FLINT_FORCE_INLINE void
ifma_run(__m512i * lo, __m512i * hi, const ulong * bk, nn_srcptr a,
         slong da, slong db, int rlo, int rhi, int with_hi)
{
    slong d;
    int r;

    for (d = da; d < db; d++)
    {
        __m512i w = _mm512_loadu_si512((const void *) (bk + d));

        for (r = rlo; r <= rhi; r++)
        {
            __m512i x = _mm512_set1_epi64(a[8 * r - d]);

            lo[r] = IFMA_LO(lo[r], x, w);
            if (with_hi)
                hi[r] = IFMA_HI(hi[r], x, w);
        }
    }
}

/* reduce and store the vector r of a block */
FLINT_FORCE_INLINE void
ifma_store(nn_ptr res, slong nout, int r, __m512i lo, __m512i hi,
           const ifma_ctx_struct * C, int with_hi)
{
    __m512i v;

    if (C->redc == 1)
        v = redc52_2(lo, hi, &C->rc);
    else if (C->redc == 2)
        v = _mm512_and_si512(with_hi ? _mm512_add_epi64(lo, _mm512_slli_epi64(hi, 52)) : lo, C->mask);
    else
        v = with_hi ? u52_red2(lo, hi, &C->uc) : u52_red(lo, &C->uc);

    _mm512_mask_storeu_epi64(res + 8 * r, ifma_store_mask(nout, r), v);
}

#define IFMA_STORE(r) \
    if ((r) < R && (r) <= rmax) \
        ifma_store(res, nout, r, lo[r], hi[r], C, with_hi);

/* the nout <= 8R outputs from k0 */
FLINT_FORCE_INLINE void
ifma_block(nn_ptr res, slong nout, nn_srcptr a, slong la, slong k0,
           const ifma_ctx_struct * C, int with_hi, int R)
{
    __m512i lo[8], hi[8];
    const ulong * bk;
    slong da, db, dend, rlo, rhi, rmax, t;
    int r;

    for (r = 0; r < 8; r++)
        lo[r] = hi[r] = _mm512_setzero_si512();

    IFMA_BLOCK_RANGE
    bk = ifma_bk(C, k0, da, dend);

    for ( ; da < dend; da = db)
    {
        IFMA_RUN_RANGE
#define IFMA_RUN(RLO, RHI) ifma_run(lo, hi, bk, a, da, db, RLO, RHI, with_hi)
        IFMA_SWITCH(R, rlo * 8 + rhi)
#undef IFMA_RUN
    }

    /* constant indices: clang keeps the accumulators in memory otherwise */
    IFMA_STORE(0) IFMA_STORE(1) IFMA_STORE(2) IFMA_STORE(3)
    IFMA_STORE(4) IFMA_STORE(5) IFMA_STORE(6) IFMA_STORE(7)
}

#define IFMA_BLOCKS(NAME, WITH_HI, R) \
static void \
NAME(nn_ptr res, nn_srcptr a, slong la, slong nlo, slong nhi, \
     const ifma_ctx_struct * C) \
{ \
    slong k0; \
    for (k0 = nlo; k0 < nhi; k0 += 8 * (R)) \
        ifma_block(res + k0 - nlo, FLINT_MIN(8 * (R), nhi - k0), a, la, k0, \
                   C, WITH_HI, R); \
}

IFMA_BLOCKS(ifma_blocks_lohi_r8, 1, 8)
IFMA_BLOCKS(ifma_blocks_lohi_r6, 1, 6)
IFMA_BLOCKS(ifma_blocks_lo_r8, 0, 8)
IFMA_BLOCKS(ifma_blocks_lo_r6, 0, 6)

/* as: room for la coefficients; inlined: 1-3 ns less at 8 x 8 */
FLINT_FORCE_INLINE void
_ifma52(nn_ptr res, nn_srcptr a, slong la, slong nlo, slong nhi,
        ifma_ctx_struct * C, nn_ptr as, nmod_t mod)
{
    int with_hi = (mod.n > (UWORD(1) << 26));

    C->redc = with_hi && (mod.n & 1);
    if ((mod.n & (mod.n - 1)) == 0)
    {
        C->redc = 2;
        C->mask = _mm512_set1_epi64(mod.n - 1);
    }
    else if (C->redc)
    {
        ulong c = nmod_set_ui(UWORD(1) << 52, mod);
        _nmod_vec_scalar_mul_nmod(as, a, la, nmod_mul(c, c, mod), mod);
        a = as;
        redc52_consts_init(&C->rc, mod.n);
    }
    else
    {
        u52_ctx_struct ctx;
        u52_ctx_init(&ctx, mod.n);
        C->uc = u52_consts_init(&ctx);
    }

    /* 6 vectors per block for short operands (prototype, Oct 2026) */
    if (with_hi)
    {
        if (la >= 20)
            ifma_blocks_lohi_r8(res, a, la, nlo, nhi, C);
        else
            ifma_blocks_lohi_r6(res, a, la, nlo, nhi, C);
    }
    else
    {
        if (la >= 24)
            ifma_blocks_lo_r8(res, a, la, nlo, nhi, C);
        else
            ifma_blocks_lo_r6(res, a, la, nlo, nhi, C);
    }
}

/* the margin of the copies of the ends of b: a block reads b from
   k0 - la + 1 to k0 + 8R - 1 */
#define IFMA_MARGIN 72

void
_nmod_poly_mulmid_ifma(nn_ptr res, nn_srcptr poly1, slong len1,
                       nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod)
{
    ifma_ctx_struct C;
    nn_srcptr a;
    nn_ptr head, tail, as;
    slong la, lb, tlen;
    TMP_INIT;

    if (len1 >= len2)
    {
        a = poly2; la = len2; C.b = poly1; lb = len1;
    }
    else
    {
        a = poly1; la = len1; C.b = poly2; lb = len2;
    }

    la = FLINT_MIN(la, nhi);
    lb = FLINT_MIN(lb, nhi);
    C.lb = lb;

    FLINT_ASSERT(la <= 511);
    FLINT_ASSERT(mod.n <= (UWORD(1) << 52));
    FLINT_ASSERT(0 <= nlo && nlo < nhi && nhi <= la + lb - 1);

    TMP_START;

    /* the head of b, with 8 zeros before and after; all of b if short */
    C.hend = FLINT_MIN(lb, la + IFMA_MARGIN);
    tlen = (C.hend < lb) ? lb - FLINT_MAX(0, lb - la - IFMA_MARGIN) : 0;
    head = TMP_ALLOC((C.hend + 16 + (tlen ? tlen + 8 : 0) + la) * sizeof(ulong));
    _nmod_vec_zero(head, 8);
    memcpy(head + 8, C.b, C.hend * sizeof(ulong));
    _nmod_vec_zero(head + 8 + C.hend, 8);
    C.head = head;

    /* the tail of b, with 8 zeros after */
    tail = head + C.hend + 16;
    if (tlen)
    {
        C.tbeg = lb - tlen;
        memcpy(tail, C.b + C.tbeg, tlen * sizeof(ulong));
        _nmod_vec_zero(tail + tlen, 8);
        C.tail = tail;
        as = tail + tlen + 8;
    }
    else
    {
        C.tbeg = 0;
        C.tail = NULL;
        as = tail;
    }

    _ifma52(res, a, la, nlo, nhi, &C, as, mod);

    TMP_END;
}

#else

/* no AVX512-IFMA at compile time (or a 32-bit build) */
void
_nmod_poly_mulmid_ifma(nn_ptr res, nn_srcptr poly1, slong len1,
                       nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod)
{
    _nmod_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi, mod);
}

#endif
