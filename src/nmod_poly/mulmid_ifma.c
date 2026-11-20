/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "longlong.h"
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
    Coefficients nlo to nhi - 1 of the product, for any n and a shorter
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

    For n > 2^52, the sum is accumulated in three limbs of weights 1, 2^52
    and 2^104, in one of two ways.
    - Seven IFMA per product of vectors: with x = x0 + 2^52 x1 and
      y = y0 + 2^52 y1 (x1, y1 < 2^12), the limbs receive lo(x0 y0);
      hi(x0 y0) + lo(x0 y1) + lo(x1 y0); hi(x0 y1) + hi(x1 y0) + lo(x1 y1),
      where lo and hi are the 52-bit halves. With 511 terms they stay below
      2^61, 1533 * 2^52 < 2^63 and 2^34.
    - CRT, where vpmullq is fast (Zen 4), from
      FLINT_NMOD_POLY_IFMA_CRT_MIN_LEN terms (always in
      _nmod_poly_mulmid_ifma_crt): the exact sum c < 511 n^2 < 2^137 is
      computed modulo 2^64 (vpmullq, vpaddq), modulo p1 = 2^52 - 47 (two
      IFMA, low and high halves) and modulo p2 = 2^26 - 5 (one IFMA: the
      products are below 2^52), and 2^64 p1 p2 > 2^141. See
      ifma_crt_limbs.
    The three limbs V of each output are reduced, for odd n, by a two-limb
    REDC in radix 2^52, a having been multiplied by 2^104 mod n; for n a
    power of 2, by a mask; otherwise by NMOD_RED3.
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
    For X = lo + 2^52 hi, lo, hi < 2^62, a first step gives Y < 2^63, and a
    second one a result below 2^11 + n <= 2n (n >= 2^11).
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

/*
    Two-limb REDC in radix 2^52, n = n0 + 2^52 n1 odd, 2^52 < n < 2^64:
    V = v0 + 2^52 v1 + 2^104 v2 with v0, v1, v2 < 2^63 and
    V < (n - 2^12) 2^104 is mapped to V 2^-104 mod n, canonical. Each of
    the two steps adds m n to the partial sum, m = t0 n' mod 2^52 for the
    low limb t0, and drops the zero low limb.
*/
typedef struct
{
    __m512i n0, n1, nprime, mask52, one, n;
}
redc104_consts;

static void
redc104_consts_init(redc104_consts * C, ulong n)
{
    C->n0 = _mm512_set1_epi64(n & ((UWORD(1) << 52) - 1));
    C->n1 = _mm512_set1_epi64(n >> 52);
    C->nprime = _mm512_set1_epi64((-n_binvert(n)) & ((UWORD(1) << 52) - 1));
    C->mask52 = _mm512_set1_epi64((UWORD(1) << 52) - 1);
    C->one = _mm512_set1_epi64(1);
    C->n = _mm512_set1_epi64(n);
}

FLINT_FORCE_INLINE __m512i
redc104(__m512i v0, __m512i v1, __m512i v2, const redc104_consts * C)
{
    const __m512i zero = _mm512_setzero_si512();
    __m512i m, t0, t1;
    __mmask8 nz, ge;

    v1 = _mm512_add_epi64(v1, _mm512_srli_epi64(v0, 52));
    v0 = _mm512_and_si512(v0, C->mask52);
    v2 = _mm512_add_epi64(v2, _mm512_srli_epi64(v1, 52));
    v1 = _mm512_and_si512(v1, C->mask52);

    m = IFMA_LO(zero, v0, C->nprime);
    nz = _mm512_test_epi64_mask(v0, v0);
    t0 = IFMA_HI(v1, m, C->n0);
    t0 = IFMA_LO(t0, m, C->n1);
    t0 = _mm512_mask_add_epi64(t0, nz, t0, C->one);
    t1 = IFMA_HI(v2, m, C->n1);
    t1 = _mm512_add_epi64(t1, _mm512_srli_epi64(t0, 52));
    t0 = _mm512_and_si512(t0, C->mask52);

    m = IFMA_LO(zero, t0, C->nprime);
    nz = _mm512_test_epi64_mask(t0, t0);
    v0 = IFMA_HI(t1, m, C->n0);
    v0 = IFMA_LO(v0, m, C->n1);
    v0 = _mm512_mask_add_epi64(v0, nz, v0, C->one);
    v1 = IFMA_HI(zero, m, C->n1);
    v1 = _mm512_add_epi64(v1, _mm512_srli_epi64(v0, 52));
    v0 = _mm512_and_si512(v0, C->mask52);

    ge = _mm512_cmpgt_epu64_mask(v1, C->n1)
         | (_mm512_cmpeq_epu64_mask(v1, C->n1) & _mm512_cmpge_epu64_mask(v0, C->n0));
    t0 = _mm512_add_epi64(v0, _mm512_slli_epi64(v1, 52));

    return _mm512_mask_sub_epi64(t0, ge, t0, C->n);
}

/* the primes of the CRT for n > 2^52, and constants modulo p2 */
#define IFMA_CRT_P1 UWORD(4503599627370449)     /* 2^52 - 47 */
#define IFMA_CRT_P2 UWORD(67108859)             /* 2^26 - 5 */
#define IFMA_CRT_P2_2E40 UWORD(81920)           /* 2^40 mod p2 */
#define IFMA_CRT_P2_CS UWORD(15251985)          /* 2^104 / p1 mod p2 */

typedef struct
{
    nn_srcptr b;        /* the longer operand, length lb */
    slong lb;
    nn_srcptr head;     /* head[8 + j] = b[j] for j < hend, zero around */
    slong hend;
    nn_srcptr tail;     /* tail[j - tbeg] = b[j] for tbeg <= j < lb, zero after */
    slong tbeg;
    int redc;           /* 1: REDC, 2: n a power of 2, 0: u52_red (n <= 2^52)
                           or NMOD_RED3 (n > 2^52) */
    __m512i mask;       /* n - 1 for n a power of 2 */
    redc52_consts rc;
    u52_consts uc;
    /* n > 2^52 */
    nmod_t mod;
    redc104_consts rn;
    redc52_consts r1, r2;       /* the CRT primes p1, p2 */
    nn_srcptr pl[3];            /* CRT: b, b mod p1, b mod p2 from index g0 */
    slong g0;
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

/* n <= 2^52 ****************************************************************/

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

/* n > 2^52 *****************************************************************/

/*
    The vectors rlo to rhi, offsets d from da to db - 1, with bw[i][d] the
    plane i of b[k0 + d] and aw[i][j] that of a[j]:
    7 IFMA: bw[0] = b, aw[0] = a, aw[1] = a >> 52 (b >> 52 on the fly);
    CRT: bw[0] = b, bw[1] = b mod p1, bw[2] = b mod p2, aw[0] = a,
    aw[1] = 2^40 a mod p1, aw[2] = 2^40 a mod p2.
*/
FLINT_FORCE_INLINE void
ifma64_run(__m512i * c0, __m512i * c1, __m512i * c2, __m512i * c3,
           const ulong * const * bw, const ulong * const * aw,
           slong da, slong db, int rlo, int rhi, int crt)
{
    slong d;
    int r;

    for (d = da; d < db; d++)
    {
        __m512i w0 = _mm512_loadu_si512((const void *) (bw[0] + d));
        __m512i w1, w2;

        if (crt)
        {
            w1 = _mm512_loadu_si512((const void *) (bw[1] + d));
            w2 = _mm512_loadu_si512((const void *) (bw[2] + d));
        }
        else
        {
            w1 = _mm512_srli_epi64(w0, 52);
            w2 = w1;
        }

        for (r = rlo; r <= rhi; r++)
        {
            __m512i x0 = _mm512_set1_epi64(aw[0][8 * r - d]);
            __m512i x1 = _mm512_set1_epi64(aw[1][8 * r - d]);

            if (crt)
            {
                __m512i x2 = _mm512_set1_epi64(aw[2][8 * r - d]);

                c0[r] = _mm512_add_epi64(c0[r], _mm512_mullo_epi64(x0, w0));
                c1[r] = IFMA_LO(c1[r], x1, w1);
                c2[r] = IFMA_HI(c2[r], x1, w1);
                c3[r] = IFMA_LO(c3[r], x2, w2);
            }
            else
            {
                c0[r] = IFMA_LO(c0[r], x0, w0);
                c1[r] = IFMA_HI(c1[r], x0, w0);
                c1[r] = IFMA_LO(c1[r], x0, w1);
                c2[r] = IFMA_HI(c2[r], x0, w1);
                c1[r] = IFMA_LO(c1[r], x1, w0);
                c2[r] = IFMA_HI(c2[r], x1, w0);
                c2[r] = IFMA_LO(c2[r], x1, w1);
            }
        }
    }
}

/*
    The limbs P0 + 2^52 P1 + 2^104 P2 of the exact sum c from its residues:
    r0 = c mod 2^64; x1 + 2^52 x2 = 2^40 c mod p1 and x3 = 2^40 c mod p2
    (unreduced). With h = (c - r0) / 2^64 < 2^73 < p1 p2:
    h1 = h mod p1 = (x1 + 2^52 x2 - 2^40 r0) 2^-104 mod p1,
    h2 = h mod p2 = (x3 - 2^40 r0) 2^-104 mod p2 (one REDC each, after
    adding multiples of p1, p2 to keep the limbs nonnegative),
    s = (h2 - h1) / p1 mod p2, h = h1 + p1 s = u + 2^52 hH with
    u = h1 + lo52(p1 s) < 2^53, and c = r0 + 2^64 u + 2^116 hH.
*/
FLINT_FORCE_INLINE void
ifma_crt_limbs(__m512i * P0, __m512i * P1, __m512i * P2, __m512i r0,
               __m512i x1, __m512i x2, __m512i x3, const ifma_ctx_struct * C)
{
    const __m512i zero = _mm512_setzero_si512();
    const __m512i one = _mm512_set1_epi64(1);
    const __m512i two52 = _mm512_set1_epi64(UWORD(1) << 52);
    const __m512i mask52 = _mm512_set1_epi64((UWORD(1) << 52) - 1);
    const __m512i c40 = _mm512_set1_epi64(IFMA_CRT_P2_2E40);
    const __m512i p1 = _mm512_set1_epi64(IFMA_CRT_P1);
    __m512i lo, hi, h1, h2, s, u;

    /* 2^40 r0 = ((r0 mod 2^12) 2^40) + 2^52 (r0 >> 12), r0 >> 12 < 2 p1 */
    lo = _mm512_sub_epi64(_mm512_add_epi64(x1, two52),
                          _mm512_slli_epi64(_mm512_and_si512(r0, _mm512_set1_epi64(0xfff)), 40));
    hi = _mm512_sub_epi64(_mm512_add_epi64(x2, _mm512_set1_epi64(2 * IFMA_CRT_P1)),
                          _mm512_add_epi64(_mm512_srli_epi64(r0, 12), one));
    h1 = redc52_2(lo, hi, &C->r1);

    /* 2^40 r0 = c40 r0 mod p2, and c40 r0 = lo52(c40 r0') + 2^52 t with
       r0' = r0 mod 2^52 and t = hi52(c40 r0') + c40 (r0 >> 52) < 2^30
       < 2^14 p2 (c40 < 2^17) */
    lo = _mm512_sub_epi64(_mm512_add_epi64(x3, two52), IFMA_LO(zero, r0, c40));
    hi = _mm512_sub_epi64(_mm512_set1_epi64(IFMA_CRT_P2 << 14),
                          _mm512_add_epi64(IFMA_LO(IFMA_HI(zero, r0, c40),
                                           _mm512_srli_epi64(r0, 52), c40), one));
    h2 = redc52_2(lo, hi, &C->r2);

    /* s = (h2 cs + h1 (p2 - cs)) 2^-104 mod p2, cs = 2^104 / p1 mod p2 */
    lo = IFMA_LO(IFMA_LO(zero, h2, _mm512_set1_epi64(IFMA_CRT_P2_CS)),
                 h1, _mm512_set1_epi64(IFMA_CRT_P2 - IFMA_CRT_P2_CS));
    hi = IFMA_HI(zero, h1, _mm512_set1_epi64(IFMA_CRT_P2 - IFMA_CRT_P2_CS));
    s = redc52_2(lo, hi, &C->r2);

    u = _mm512_add_epi64(h1, IFMA_LO(zero, p1, s));
    *P0 = _mm512_and_si512(r0, mask52);
    *P1 = _mm512_add_epi64(_mm512_srli_epi64(r0, 52),
                           _mm512_slli_epi64(_mm512_and_si512(u, _mm512_set1_epi64((UWORD(1) << 40) - 1)), 12));
    *P2 = _mm512_add_epi64(_mm512_srli_epi64(u, 40),
                           _mm512_slli_epi64(IFMA_HI(zero, p1, s), 12));
}

/* reduce and store the vector r of a block */
FLINT_FORCE_INLINE void
ifma64_store(nn_ptr res, slong nout, int r, __m512i x0, __m512i x1,
             __m512i x2, __m512i x3, const ifma_ctx_struct * C, int crt)
{
    __m512i P0, P1, P2;
    __mmask8 mask = ifma_store_mask(nout, r);

    if (crt)
        ifma_crt_limbs(&P0, &P1, &P2, x0, x1, x2, x3, C);
    else
    {
        P0 = x0;
        P1 = x1;
        P2 = x2;
    }

    if (C->redc == 1)
    {
        _mm512_mask_storeu_epi64(res + 8 * r, mask, redc104(P0, P1, P2, &C->rn));
    }
    else if (C->redc == 2)
    {
        /* V mod 2^64 */
        P0 = _mm512_add_epi64(P0, _mm512_slli_epi64(P1, 52));
        _mm512_mask_storeu_epi64(res + 8 * r, mask, _mm512_and_si512(P0, C->mask));
    }
    else
    {
        ulong y0[8], y1[8], y2[8];
        slong j;

        _mm512_storeu_si512((void *) y0, P0);
        _mm512_storeu_si512((void *) y1, P1);
        _mm512_storeu_si512((void *) y2, P2);

        for (j = 0; j < 8 && 8 * r + j < nout; j++)
        {
            ulong w0, w1, w2;

            /* V = y0 + 2^52 y1 + 2^104 y2 in three words, w2 < 2^10 */
            add_ssaaaa(w1, w0, UWORD(0), y0[j], y1[j] >> 12, y1[j] << 52);
            add_ssaaaa(w2, w1, y2[j] >> 24, y2[j] << 40, UWORD(0), w1);
            NMOD_RED3(res[8 * r + j], w2, w1, w0, C->mod);
        }
    }
}

#define IFMA64_STORE(r) \
    if ((r) < R && (r) <= rmax) \
        ifma64_store(res, nout, r, c0[r], c1[r], c2[r], c3[r], C, crt);

/* the nout <= 8R outputs from k0 */
FLINT_FORCE_INLINE void
ifma64_block(nn_ptr res, slong nout, const ulong * const * aw, slong la,
             slong k0, const ifma_ctx_struct * C, int crt, int R)
{
    __m512i c0[8], c1[8], c2[8], c3[8];
    const ulong * bw[3];
    slong da, db, dend, rlo, rhi, rmax, t;
    int r;

    for (r = 0; r < 8; r++)
        c0[r] = c1[r] = c2[r] = c3[r] = _mm512_setzero_si512();

    IFMA_BLOCK_RANGE

    if (crt)
    {
        bw[0] = C->pl[0] + (k0 - C->g0);
        bw[1] = C->pl[1] + (k0 - C->g0);
        bw[2] = C->pl[2] + (k0 - C->g0);
    }
    else
    {
        bw[0] = ifma_bk(C, k0, da, dend);
        bw[1] = bw[2] = bw[0];
    }

    for ( ; da < dend; da = db)
    {
        IFMA_RUN_RANGE
#define IFMA_RUN(RLO, RHI) ifma64_run(c0, c1, c2, c3, bw, aw, da, db, RLO, RHI, crt)
        IFMA_SWITCH(R, rlo * 8 + rhi)
#undef IFMA_RUN
    }

    IFMA64_STORE(0) IFMA64_STORE(1) IFMA64_STORE(2) IFMA64_STORE(3)
    IFMA64_STORE(4) IFMA64_STORE(5) IFMA64_STORE(6) IFMA64_STORE(7)
}

#define IFMA64_BLOCKS(NAME, CRT, R) \
static void \
NAME(nn_ptr res, const ulong * const * aw, slong la, slong nlo, slong nhi, \
     const ifma_ctx_struct * C) \
{ \
    slong k0; \
    for (k0 = nlo; k0 < nhi; k0 += 8 * (R)) \
        ifma64_block(res + k0 - nlo, FLINT_MIN(8 * (R), nhi - k0), aw, la, \
                     k0, C, CRT, R); \
}

IFMA64_BLOCKS(ifma64_blocks_7_r8, 0, 8)
IFMA64_BLOCKS(ifma64_blocks_7_r6, 0, 6)
IFMA64_BLOCKS(ifma64_blocks_crt_r6, 1, 6)

/* the outputs of the CRT blocks by chunks of IFMA_CRT_CHUNK, with the
   planes of the part of b they read */
#define IFMA_CRT_CHUNK (8 * 6 * 32)

/* as: room for 3 la coefficients; the CRT for la >= crt_min_len > 0 */
static void
_ifma64(nn_ptr res, nn_srcptr a, slong la, slong nlo, slong nhi,
        ifma_ctx_struct * C, nn_ptr as, nmod_t mod, slong crt_min_len)
{
    const ulong * aw[3];
    nn_ptr a0 = as, a1 = as + la, a2 = as + 2 * la;
    int crt = (crt_min_len > 0 && la >= crt_min_len);
    slong i;

    C->mod = mod;

    /* a multiplied by 2^104 mod n for the REDC */
    if ((mod.n & (mod.n - 1)) == 0)
    {
        C->redc = 2;
        C->mask = _mm512_set1_epi64(mod.n - 1);
        _nmod_vec_set(a0, a, la);
    }
    else if (mod.n & 1)
    {
        ulong c = nmod_set_ui(UWORD(1) << 52, mod);
        C->redc = 1;
        redc104_consts_init(&C->rn, mod.n);
        _nmod_vec_scalar_mul_nmod(a0, a, la, nmod_mul(c, c, mod), mod);
    }
    else
    {
        C->redc = 0;
        _nmod_vec_set(a0, a, la);
    }

    aw[0] = a0;
    aw[1] = a1;
    aw[2] = a2;

    if (!crt)
    {
        for (i = 0; i < la; i++)
            a1[i] = a0[i] >> 52;

        if (la >= 20)
            ifma64_blocks_7_r8(res, aw, la, nlo, nhi, C);
        else
            ifma64_blocks_7_r6(res, aw, la, nlo, nhi, C);
    }
    else
    {
        nmod_t m1, m2;
        nn_ptr pl;
        slong k, kend, K, g0, plen, j0, j1, lb = C->lb;
        TMP_INIT;

        nmod_init(&m1, IFMA_CRT_P1);
        nmod_init(&m2, IFMA_CRT_P2);
        redc52_consts_init(&C->r1, IFMA_CRT_P1);
        redc52_consts_init(&C->r2, IFMA_CRT_P2);

        _nmod_vec_reduce(a1, a0, la, m1);
        _nmod_vec_scalar_mul_nmod(a1, a1, la, UWORD(1) << 40, m1);
        _nmod_vec_reduce(a2, a0, la, m2);
        _nmod_vec_scalar_mul_nmod(a2, a2, la, IFMA_CRT_P2_2E40, m2);

        /* chunks of K outputs, a multiple of the 48 of a block: the chunk
           of outputs k to kend - 1 reads b from k - la + 1 to at most
           k + K - 1, kept at offsets j - g0 from g0 = k - la + 1 - 8, with
           zeros outside b */
        K = FLINT_MIN(IFMA_CRT_CHUNK, (nhi - nlo + 47) / 48 * 48);
        plen = K + la + 16;
        TMP_START;
        pl = TMP_ALLOC(3 * plen * sizeof(ulong));

        for (k = nlo; k < nhi; k = kend)
        {
            kend = FLINT_MIN(nhi, k + K);
            g0 = k - la + 1 - 8;
            j0 = FLINT_MIN(FLINT_MAX(g0, 0), g0 + plen);
            j1 = FLINT_MAX(FLINT_MIN(g0 + plen, lb), j0);

            for (i = 0; i < 3; i++)
            {
                _nmod_vec_zero(pl + i * plen, j0 - g0);
                _nmod_vec_zero(pl + i * plen + (j1 - g0), plen - (j1 - g0));
            }

            _nmod_vec_set(pl + (j0 - g0), C->b + j0, j1 - j0);
            _nmod_vec_reduce(pl + plen + (j0 - g0), C->b + j0, j1 - j0, m1);
            _nmod_vec_reduce(pl + 2 * plen + (j0 - g0), C->b + j0, j1 - j0, m2);

            C->pl[0] = pl;
            C->pl[1] = pl + plen;
            C->pl[2] = pl + 2 * plen;
            C->g0 = g0;

            ifma64_blocks_crt_r6(res + (k - nlo), aw, la, k, kend, C);
        }

        TMP_END;
    }
}

/* the margin of the copies of the ends of b: a block reads b from
   k0 - la + 1 to k0 + 8R - 1 */
#define IFMA_MARGIN 72

/* inlined in the two functions below: 2-4 ns less per call at lengths 8
   to 16 (Zen 4, Ice Lake); mod by pointer, as gcc 14 and 15 read a copy by
   value with one 16-byte load from the two 8-byte stores of the
   dispatchers, which waits for these stores (about 10 ns per call) */
FLINT_FORCE_INLINE void
_mulmid_ifma(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2,
             slong len2, slong nlo, slong nhi, const nmod_t * mod,
             slong crt_min_len)
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
    FLINT_ASSERT(0 <= nlo && nlo < nhi && nhi <= la + lb - 1);

    TMP_START;

    /* the head of b, with 8 zeros before and after; all of b if short */
    C.hend = FLINT_MIN(lb, la + IFMA_MARGIN);
    tlen = (C.hend < lb) ? lb - FLINT_MAX(0, lb - la - IFMA_MARGIN) : 0;
    head = TMP_ALLOC((C.hend + 16 + (tlen ? tlen + 8 : 0) + 3 * la) * sizeof(ulong));
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

    if (mod->n <= (UWORD(1) << 52))
        _ifma52(res, a, la, nlo, nhi, &C, as, *mod);
    else
        _ifma64(res, a, la, nlo, nhi, &C, as, *mod, crt_min_len);

    TMP_END;
}

void
_nmod_poly_mulmid_ifma(nn_ptr res, nn_srcptr poly1, slong len1,
                       nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod)
{
    _mulmid_ifma(res, poly1, len1, poly2, len2, nlo, nhi, &mod,
                 FLINT_NMOD_POLY_IFMA_CRT_MIN_LEN);
}

void
_nmod_poly_mulmid_ifma_crt(nn_ptr res, nn_srcptr poly1, slong len1,
                           nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod)
{
    _mulmid_ifma(res, poly1, len1, poly2, len2, nlo, nhi, &mod, 1);
}

#else

/* no AVX512-IFMA at compile time (or a 32-bit build) */
void
_nmod_poly_mulmid_ifma(nn_ptr res, nn_srcptr poly1, slong len1,
                       nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod)
{
    _nmod_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi, mod);
}

void
_nmod_poly_mulmid_ifma_crt(nn_ptr res, nn_srcptr poly1, slong len1,
                           nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod)
{
    _nmod_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi, mod);
}

#endif
