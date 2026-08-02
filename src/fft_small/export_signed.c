/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <gmp.h>
#include "mpn_extras.h"
#include "fft_small.h"
#include "machine_vectors.h"
#include "crt_helpers.h"
#include "thread_pool.h"
#include "thread_support.h"

/* defined in mpn_helpers.c */
void _convert_block(ulong* Xs, sd_fft_ctx_struct* Rffts, double* d,
                    ulong dstride, ulong np, ulong I);

/*
    One block of slots, reconstructed and glued into the window.

    Everything in the inner loop is specialized on the geometry, exactly
    as the unsigned reconstruction specializes its own: the Chinese
    remainder dot product uses the same tuned even/odd-split kernels, and
    the centered lift, the shift into place and the signed accumulation
    into the window are unrolled at the coefficient length with carry
    intrinsics rather than calls into mpn with runtime lengths. The
    latter was what made this loop several times slower than the
    templated unsigned one.

    N is the coefficient length and M the number of significant cofactor
    limbs; M is a function of np, not of N alone (np = 4 and np = 5 share
    N = 4 but need M = 3 and M = 4), so the table is keyed on np.
*/
typedef void (* sig_block_func)(ulong * W, ulong wlen, ulong A, ulong bits,
        const crt_data_struct * Cd, const ulong * Xs, nn_srcptr prod,
        nn_srcptr prodh, ulong slot0, ulong jlo, ulong jhi);

#define DEFINE_SIG_BLOCK(NP, N, M) \
static void CAT(_sig_block, NP)(ulong * W, ulong wlen, ulong A, ulong bits, \
        const crt_data_struct * Cd, const ulong * Xs, nn_srcptr prod, \
        nn_srcptr prodh, ulong slot0, ulong jlo, ulong jhi) \
{ \
    for (ulong j = jlo; j < jhi; j++) \
    { \
        ulong r[N + 2], t[N + 2]; \
        ulong slot = slot0 + j; \
        ulong off = slot * bits; \
        ulong k = off / FLINT_BITS - A; \
        ulong sh = off % FLINT_BITS; \
        ulong l = 0, acc = 0; \
        slong q; \
        unsigned char cf; \
        int neg = 0; \
 \
        CAT3(_big_mul, N, M)(r, t, _crt_data_co_prime( \
                (crt_data_struct *) Cd, l, N), Xs[l * BLK_SZ + j]); \
        for (l++; l < NP; l++) \
            CAT3(_big_addmul, N, M)(r, t, _crt_data_co_prime( \
                    (crt_data_struct *) Cd, l, N), Xs[l * BLK_SZ + j]); \
        CAT(_reduce_big_sum, N)(r, t, (ulong *) prod); \
 \
        for (q = 0; q < N; q++) \
            acc |= r[q]; \
        if (acc == 0) \
            continue; \
 \
        /* centered lift */ \
        for (q = N - 1; q >= 0; q--) \
        { \
            if (r[q] != prodh[q]) \
            { \
                neg = (r[q] > prodh[q]); \
                break; \
            } \
        } \
        if (neg) \
        { \
            cf = 0; \
            for (q = 0; q < N; q++) \
                cf = _subborrow_ulong(cf, prod[q], r[q], &r[q]); \
        } \
 \
        /* shift into position within the limb */ \
        if (sh != 0) \
        { \
            r[N] = r[N - 1] >> (FLINT_BITS - sh); \
            for (q = N; q >= 2; q--) \
                r[q - 1] = (r[q - 1] << sh) | (r[q - 2] >> (FLINT_BITS - sh)); \
            r[0] = r[0] << sh; \
        } \
        else \
            r[N] = 0; \
 \
        /* signed accumulation into the window */ \
        if (!neg) \
        { \
            cf = 0; \
            for (q = 0; q <= N; q++) \
                cf = _addcarry_ulong(cf, W[k + q], r[q], &W[k + q]); \
            for (q = k + N + 1; cf && q < (slong) wlen; q++) \
                cf = _addcarry_ulong(0, W[q], UWORD(1), &W[q]); \
        } \
        else \
        { \
            cf = 0; \
            for (q = 0; q <= N; q++) \
                cf = _subborrow_ulong(cf, W[k + q], r[q], &W[k + q]); \
            for (q = k + N + 1; cf && q < (slong) wlen; q++) \
                cf = _subborrow_ulong(0, W[q], UWORD(1), &W[q]); \
        } \
    } \
}

DEFINE_SIG_BLOCK(4, 4, 3)
DEFINE_SIG_BLOCK(5, 4, 4)
DEFINE_SIG_BLOCK(6, 5, 4)
DEFINE_SIG_BLOCK(7, 6, 5)
DEFINE_SIG_BLOCK(8, 7, 6)

#undef DEFINE_SIG_BLOCK

static const sig_block_func _sig_block_tab[8 - 4 + 1] = {
    _sig_block_4, _sig_block_5, _sig_block_6, _sig_block_7, _sig_block_8
};

/*
    Signed counterpart of fft_small_export_mpn, for chunk values that may
    be negative (bounded in magnitude by less than half the prime
    product). The reconstruction is the composition

        prod_p (Z/pZ)[x] -> (CRT) -> Z[x] -> (evaluate at 2^bits) -> Z

    fused into a single pass: per slot, the Chinese remainder residue
    r in [0, prod) is lifted to the centered remainder

        r,          r <= prod/2
        r - prod,   r >  prod/2

    -- exactly the chunk's signed value, read where the wrap information
    still exists -- and the evaluation glues the signed chunks with
    borrow propagation. (Recovering the signs after an unsigned export is
    impossible even in principle: overlapping slot contributions destroy
    the wrap bits, and distinct in-range signed chunk vectors with
    different values can produce the same unsigned output integer.)

    The gluing uses a sliding two's-complement window: consecutive slots
    overlap by 64 coeff_len - bits bits, and the running prefix sum
    shifted past the current slot is bounded by twice the chunk magnitude
    bound, so all carry and borrow propagation is confined to the window;
    limbs below it are final and are flushed as they complete. The window
    carries slack so the slide is amortized over many slots.

    lo_limbs > 0 requests a truncated export of the limbs starting at
    limb lo_limbs of the value: only slots that can influence the
    requested window are processed, and the discarded low tail perturbs
    the result by at most one unit in the lowest returned limb (the
    standard mulhigh contract). With lo_limbs = 0 the export is exact.

    Threading follows the multiplication driver's crt stage: the slot
    range is split at BLK_SZ boundaries, each segment reconstructs into
    its own window and writes its own band of output limbs, and the
    signed spill each segment leaves above its band is added back in a
    serial stitch afterwards. The sign is resolved once, from the top of
    the assembled two's-complement result.

    z receives zn limbs of |value| starting at limb lo_limbs, zero padded
    at the top; *sign receives the sign. The first nslots chunk slots are
    read; the caller must provide 64 (lo_limbs + zn) >=
    nslots * bits + 64 coeff_len + 2.
*/

/* The window holds one whole block of slots plus a coefficient, so the
   slide runs once per block: a block advances BLK_SZ * bits / 64 limbs,
   which is a few thousand bytes -- comfortably L1-resident, and the
   flush it produces is exactly the output rate. */

/* one slot range: window reconstruction into the output band
   [Lbase, Lstop), leaving the spill above Lstop in w->spill */
typedef struct
{
    nn_ptr z;
    ulong zn;
    ulong lo_limbs;
    const fft_small_op_struct * X;
    const fft_small_plan_struct * P;
    ulong slot_start;
    ulong slot_stop;
    ulong Lbase;
    ulong Lstop;
    int is_last;
    nn_ptr spill;
    ulong spill_alloc;
    ulong spill_len;
    ulong spill_base;
}
sig_seg_struct;

static void
_sig_segment(sig_seg_struct * S)
{
    const fft_small_plan_struct * P = S->P;
    ulong bits = P->bits;
    ulong np = P->np;
    const crt_data_struct * Cd = P->crts + np - 1;
    ulong n = Cd->coeff_len;
    nn_srcptr prod = crt_data_prod_primes((crt_data_struct *) Cd);
    ulong prodh[MPN_CTX_NCRTS + 2];
    nn_ptr W;
    ulong Xs[MPN_CTX_NCRTS * BLK_SZ];
    ulong wlive = n + 3;
    ulong wlen = wlive + n_cdiv((ulong) BLK_SZ * bits, FLINT_BITS) + 8;
    ulong A = S->Lbase;
    ulong slot, k, idx;
    sig_block_func sig_block = _sig_block_tab[np - 4];

    mpn_rshift(prodh, prod, n, 1);
    W = flint_malloc(wlen * sizeof(ulong));
    flint_mpn_zero(W, wlen);

#define SIG_PUT(absidx, val) \
    do { \
        ulong _a = (absidx); \
        if (_a >= S->lo_limbs && _a - S->lo_limbs < S->zn) \
            S->z[_a - S->lo_limbs] = (val); \
    } while (0)

    /* one block at a time: the slide is checked per block (a block
       advances BLK_SZ * bits / 64 limbs, well inside the slack), and the
       block itself runs in the templated kernel */
    for (slot = S->slot_start; slot < S->slot_stop; )
    {
        ulong blk = slot / BLK_SZ;
        ulong bstart = blk * BLK_SZ;
        ulong jlo = slot - bstart;
        ulong jhi = n_min(BLK_SZ, S->slot_stop - bstart);
        ulong klast;

        _convert_block(Xs, P->ffts, S->X->data, P->stride, np, blk);

        /* room for every slot of this block plus a coefficient; the
           slide is by exactly this slot's offset, since only limbs
           strictly below it are final -- a fixed stride could jump past
           the first slot of the block and wrap its window index */
        klast = ((bstart + jhi - 1) * bits) / FLINT_BITS - A;
        if (klast + wlive > wlen)
        {
            ulong e = (W[wlen - 1] >> (FLINT_BITS - 1)) ? ~UWORD(0) : 0;
            ulong i2;
            ulong s0 = (slot * bits) / FLINT_BITS - A;


            FLINT_ASSERT(s0 <= wlen);
            FLINT_ASSERT(klast - s0 + wlive <= wlen);
            for (i2 = 0; i2 < s0; i2++)
                SIG_PUT(A + i2, W[i2]);
            memmove(W, W + s0, (wlen - s0) * sizeof(ulong));
            for (i2 = wlen - s0; i2 < wlen; i2++)
                W[i2] = e;
            A += s0;
        }

        sig_block(W, wlen, A, bits, Cd, Xs, prod, prodh, bstart, jlo, jhi);
        slot = bstart + jhi;
    }

    /* flush the part of the window below this segment's end, and hand
       the rest on as a signed spill */
    {
        ulong e = (W[wlen - 1] >> (FLINT_BITS - 1)) ? ~UWORD(0) : 0;
        ulong stop = S->is_last ? (S->lo_limbs + S->zn) : S->Lstop;

        for (idx = 0; idx < wlen && A + idx < stop; idx++)
            SIG_PUT(A + idx, W[idx]);

        if (S->is_last)
        {
            /* sign extension fills the remainder of the output */
            for (; A + idx < stop; idx++)
                SIG_PUT(A + idx, e);
            S->spill_len = 0;
        }
        else
        {
            ulong off = (S->Lstop > A) ? S->Lstop - A : 0;
            ulong m = (off < wlen) ? wlen - off : 0;

            S->spill = flint_malloc(wlen * sizeof(ulong));
            S->spill_alloc = wlen;

            for (k = 0; k < m; k++)
                S->spill[k] = W[off + k];
            for (; k < wlen; k++)
                S->spill[k] = e;
            S->spill_len = wlen;
            S->spill_base = S->Lstop;
        }
    }

    flint_free(W);
#undef SIG_PUT
}

/* add a two's-complement spill (slen limbs, sign extended above) into
   the output at absolute limb Labs */
static void
_sig_add_spill(nn_ptr z, ulong zn, ulong lo_limbs, ulong Labs,
               nn_srcptr spill, ulong slen)
{
    int neg;
    slong off;

    if (slen == 0)
        return;

    neg = (spill[slen - 1] >> (FLINT_BITS - 1)) != 0;
    off = (slong) Labs - (slong) lo_limbs;

    if (off >= (slong) zn)
        return;

    if (off < 0)
    {
        ulong skip = (ulong) (-off);
        if (skip >= slen)
            return;
        spill += skip;
        slen -= skip;
        off = 0;
    }

    if (slen > zn - (ulong) off)
        slen = zn - (ulong) off;

    mpn_add(z + off, z + off, zn - (ulong) off, spill, slen);

    /* sign extension above the spill: all ones, i.e. a borrow one limb
       past its top */
    if (neg && (ulong) off + slen < zn)
        mpn_sub_1(z + off + slen, z + off + slen,
                  zn - (ulong) off - slen, 1);
}

typedef struct
{
    sig_seg_struct * S;
}
sig_worker_arg;

static void
_sig_worker(void * varg)
{
    _sig_segment(((sig_worker_arg *) varg)->S);
}

ulong _sig_max_threads_used = 0;   /* TEMPORARY */

void
fft_small_export_mpn_signed_trunc(ulong* z, ulong zn, int * sign,
        const fft_small_op_t X, ulong nslots, ulong lo_limbs,
        const fft_small_plan_t P)
{
    ulong bits = P->bits;
    ulong n = (P->crts + P->np - 1)->coeff_len;
    ulong j0;
    thread_pool_handle * handles = NULL;
    slong nworkers = 0;
    ulong nthreads, l;

    FLINT_ASSERT(4 <= P->np && P->np <= 8);
    FLINT_ASSERT(P->offset == 0);
    FLINT_ASSERT(X->domain == FFT_SMALL_OP_PRODUCT);
    FLINT_ASSERT(FLINT_BITS * (lo_limbs + zn)
                 >= nslots * bits + FLINT_BITS * n + 2);
    {
        static const ulong _clen_of_np[8 - 4 + 1] = {4, 4, 5, 6, 7};
        FLINT_ASSERT(n == _clen_of_np[P->np - 4]);
        (void) _clen_of_np;
    }

    /* first slot that can influence limbs >= lo_limbs */
    j0 = 0;
    if (FLINT_BITS * lo_limbs > FLINT_BITS * n + 1)
        j0 = (FLINT_BITS * lo_limbs - FLINT_BITS * n - 1 + bits - 1) / bits;
    j0 = n_min(j0, nslots);
    j0 = (j0 / BLK_SZ) * BLK_SZ;

    if (j0 >= nslots)
    {
        flint_mpn_zero(z, zn);
        *sign = 0;
        return;
    }

    /* the segments cover every limb from the first segment's base to the
       top of the window, each written exactly once; only the prefix
       below that base (truncated exports) needs clearing */
    {
        ulong Lfirst = (j0 * bits) / FLINT_BITS;
        if (Lfirst > lo_limbs)
            flint_mpn_zero(z, n_min(Lfirst - lo_limbs, zn));
    }

    /* the reconstruction is O(output); thread it from roughly the same
       size the multiplication driver does */
    if (nslots - j0 >= 4 * BLK_SZ &&
            (nslots - j0) * bits >= 2048 * FLINT_BITS)
        nworkers = flint_request_threads(&handles, 8);
    nthreads = nworkers + 1;
    if (nthreads > _sig_max_threads_used) _sig_max_threads_used = nthreads;

    {
        sig_seg_struct * segs = FLINT_ARRAY_ALLOC(nthreads, sig_seg_struct);
        sig_worker_arg * args = FLINT_ARRAY_ALLOC(nthreads, sig_worker_arg);
        ulong span = nslots - j0;

        for (l = 0; l < nthreads; l++)
        {
            sig_seg_struct * S = segs + l;
            ulong a = j0 + n_round_up((l + 0) * span / nthreads, BLK_SZ);
            ulong b = j0 + n_round_up((l + 1) * span / nthreads, BLK_SZ);

            a = n_min(a, nslots);
            b = n_min(b, nslots);

            S->z = z;
            S->zn = zn;
            S->lo_limbs = lo_limbs;
            S->X = X;
            S->P = P;
            S->slot_start = a;
            S->slot_stop = b;
            S->Lbase = (a * bits) / FLINT_BITS;
            S->Lstop = (b * bits) / FLINT_BITS;
            S->is_last = (l + 1 == nthreads);
            S->spill = NULL;
            S->spill_alloc = 0;
            S->spill_len = 0;
            S->spill_base = S->Lstop;
            args[l].S = S;
        }

        for (l = nworkers; l > 0; l--)
            thread_pool_wake(global_thread_pool, handles[l - 1], 0,
                             _sig_worker, args + l);
        _sig_segment(segs + 0);
        for (l = nworkers; l > 0; l--)
            thread_pool_wait(global_thread_pool, handles[l - 1]);

        /* serial stitch of the signed spills across segment boundaries */
        for (l = 0; l + 1 < nthreads; l++)
            _sig_add_spill(z, zn, lo_limbs, segs[l].spill_base,
                           segs[l].spill, segs[l].spill_len);

        for (l = 0; l < nthreads; l++)
            flint_free(segs[l].spill);
        flint_free(segs);
        flint_free(args);
    }

    flint_give_back_threads(handles, nworkers);

    /* one sign resolution for the assembled two's-complement result */
    if (zn > 0 && (z[zn - 1] >> (FLINT_BITS - 1)))
    {
        *sign = 1;
        mpn_com(z, z, zn);
        mpn_add_1(z, z, zn, 1);
    }
    else
        *sign = 0;

    if (*sign && flint_mpn_zero_p(z, zn))
        *sign = 0;
}

void
fft_small_export_mpn_signed(ulong* z, ulong zn, int * sign,
        const fft_small_op_t X, ulong nslots, const fft_small_plan_t P)
{
    fft_small_export_mpn_signed_trunc(z, zn, sign, X, nslots, 0, P);
}
