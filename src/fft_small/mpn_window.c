/*
    Copyright (C) 2026 Fredrik Johansson

    Developed using Claude Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "mpn_extras.h"
#include "thread_pool.h"
#include "thread_support.h"
#include "fft_small.h"
#include "machine_vectors.h"

/*
    Window plans. The integers are split into chunks of P->bits bits and
    multiplied as chunk sequences. Only the coefficients c_lo <= k < c_hi
    can reach the limbs [lo, hi), where c_lo is the first coefficient
    whose reconstruction (coeff_len limbs) reaches limb lo and c_hi the
    first lying entirely at or above limb hi. When a power of two 2^d
    with max(operand chunks, c_hi) <= 2^d <= zlen and zlen - 2^d <= c_lo
    exists, the transforms are cyclic of length 2^d: the coefficients
    k >= 2^d alias onto k - 2^d < c_lo, which are never read as exact
    values.

    Error of the window exports. They read whole blocks of slots, so some
    slots below c_lo -- dropped coefficients, or aliased sums of two --
    are read. Every slot value is lifted to the centered range of the
    prime product, below 2^(FLINT_BITS coeff_len - 1) in absolute value,
    and a slot k < c_lo starts below bit FLINT_BITS (lo - coeff_len), so
    whatever it holds lies below limb lo; such slots hold less than
    B^lo / 2 in total either way, both as computed and in the exact
    value, and the window is within one unit of floor(x / B^lo).
*/
int
fft_small_plan_init_mpn_window(fft_small_plan_t P, mpn_ctx_t R,
                    ulong an_max, ulong bn_max, ulong len_bound,
                    int is_signed, ulong lo, ulong hi)
{
    ulong bits, nn, alen, blen, zlen, atrunc, btrunc, c_lo, c_hi, d, w;

    if (!fft_small_plan_init_mpn(P, R, an_max, bn_max, len_bound, is_signed))
        return 0;

    bits = P->bits;
    nn = (P->crts + P->np - 1)->coeff_len;
    alen = n_cdiv(FLINT_BITS * an_max, bits);
    blen = n_cdiv(FLINT_BITS * bn_max, bits);
    zlen = alen + blen - 1;
    atrunc = n_round_up(alen, BLK_SZ);
    btrunc = n_round_up(blen, BLK_SZ);

    c_lo = (lo > nn) ? n_min(zlen, n_cdiv((lo - nn) * FLINT_BITS, bits)) : 0;
    c_hi = n_min(zlen, n_cdiv(hi * FLINT_BITS, bits));

    d = n_flog2(zlen);
    w = n_pow2(d);

    if (d >= LG_BLK_SZ && d < P->depth &&
        atrunc <= w && btrunc <= w && c_hi <= w && w <= zlen && zlen <= c_lo + w)
    {
        P->depth = d;
        P->ztrunc = w;
        P->zl = c_lo;
        P->zh = c_hi;
        P->zn = zlen;

        /* the two-prime plans defer their normalization to the export,
           which indexes its scale factors by depth */
        if (P->np == 2)
            P->stride = sd_fft_ctx_data_size(d);
        else
            _fft_small_plan_set_normalizers(P);
    }

    return 1;
}

ulong
fft_small_plan_mpn_nslots(const fft_small_plan_t P, ulong an_max,
                          ulong bn_max)
{
    ulong alen = n_cdiv(FLINT_BITS * an_max, P->bits);
    ulong blen = n_cdiv(FLINT_BITS * bn_max, P->bits);
    return n_min(P->ztrunc, alen + blen - 1);
}

int
fft_small_plan_is_cyclic(const fft_small_plan_t P)
{
    return P->zn > P->ztrunc;
}

/* truncation granularity follows the plan's depth: whole blocks at
   block depths, the full (short) transform length below one block */
static ulong
_op_trunc(ulong x, const fft_small_plan_struct * P)
{
    return (P->depth < LG_BLK_SZ) ? n_pow2(P->depth) : n_round_up(x, BLK_SZ);
}

ulong
fft_small_export_mpn_signed_limbs(const fft_small_plan_t P, ulong nslots)
{
    /* the signed export reads round_up(nslots) slots and requires
       64 zn >= nslots bits + 64 coeff_len + 2 */
    ulong clen = (P->crts + P->np - 1)->coeff_len;
    return n_cdiv(nslots * P->bits + FLINT_BITS * clen + 2 + (FLINT_BITS - 1),
                  FLINT_BITS);
}

ulong
fft_small_export_mpn_limbs(const fft_small_plan_t P, ulong len_bound,
                           ulong nslots)
{
    ulong m = n_min(_op_trunc(nslots, P), P->zn);
    return n_cdiv(m * P->bits + 2 * P->bits + FLINT_BIT_COUNT(len_bound)
                  + FLINT_BIT_COUNT(P->zn) + 4, FLINT_BITS);
}

/* d = m*d over npts points (a multiple of 16) */
static void
_scale(const sd_fft_ctx_struct * Q, double * d, ulong m_, ulong npts)
{
    vec8d m = vec8d_set_d(vec1d_reduce_0n_to_pmhn((slong) m_, Q->p));
    vec8d n = vec8d_set_d(Q->p);
    vec8d ninv = vec8d_set_d(Q->pinv);
    ulong j = 0;

    do {
        vec8d x0, x1;
        x0 = vec8d_load(d + j + 0);
        x1 = vec8d_load(d + j + 8);
        x0 = vec8d_mulmod(x0, m, n, ninv);
        x1 = vec8d_mulmod(x1, m, n, ninv);
        vec8d_store(d + j + 0, x0);
        vec8d_store(d + j + 8, x1);
    } while (j += 16, j < npts);
}

typedef struct
{
    fft_small_op_struct * X;
    const fft_small_plan_struct * P;
    ulong itrunc;
    ulong clear_to;
    const ulong * scale;
    ulong start_pi;
    ulong stop_pi;
}
_ifft_trunc_arg;

static void
_ifft_trunc_worker(void * varg)
{
    _ifft_trunc_arg * W = (_ifft_trunc_arg *) varg;
    const fft_small_plan_struct * P = W->P;
    ulong i;

    for (i = W->start_pi; i < W->stop_pi; i++)
    {
        sd_fft_ctx_struct * Q = P->ffts + P->offset + i;
        double * d = W->X->data + P->stride * i;

        sd_ifft_trunc(Q, d, P->depth, W->itrunc);
        if (W->scale != NULL && W->scale[i] != 1)
            _scale(Q, d, W->scale[i], W->itrunc);
        if (W->clear_to > W->itrunc)
            memset(d + W->itrunc, 0, (W->clear_to - W->itrunc) * sizeof(double));
    }
}

void
fft_small_ifft_trunc(fft_small_op_t X, ulong itrunc, ulong clear_to,
                     const ulong * scale, const fft_small_plan_t P)
{
    thread_pool_handle * handles = NULL;
    slong nworkers = 0;
    ulong nthreads, i;
    _ifft_trunc_arg args[MPN_CTX_NCRTS];

    itrunc = n_min(_op_trunc(itrunc, P), n_pow2(P->depth));
    clear_to = n_min(clear_to, n_pow2(P->depth));

    /* thread over the primes from the size fft_small_ifft does */
    if (P->np >= 2 && itrunc >= 3072)
        nworkers = flint_request_threads(&handles, P->np);
    nthreads = FLINT_MIN((ulong) (nworkers + 1), P->np);

    for (i = 0; i < nthreads; i++)
    {
        args[i].X = X;
        args[i].P = P;
        args[i].itrunc = itrunc;
        args[i].clear_to = clear_to;
        args[i].scale = scale;
        args[i].start_pi = (i + 0) * P->np / nthreads;
        args[i].stop_pi = (i + 1) * P->np / nthreads;
    }

    for (i = nthreads - 1; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0,
                         _ifft_trunc_worker, args + i);
    _ifft_trunc_worker(args + 0);
    for (i = nthreads - 1; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    flint_give_back_threads(handles, nworkers);

    X->domain = FFT_SMALL_OP_PRODUCT;
}

ulong *
fft_small_stage_begin(ulong * z, ulong zn, ulong need,
                      fft_small_op_pool_struct * pool, int use_slabs)
{
    if (zn >= need)
        return z;

    if (pool != NULL)
        return fft_small_op_pool_scratch_take(pool, need, use_slabs);

    return flint_malloc(need * sizeof(ulong));
}

int
fft_small_stage_end(ulong * z, ulong zn, ulong * buf, ulong len,
                    fft_small_op_pool_struct * pool)
{
    int fits = (len <= zn);

    if (buf != z)
    {
        if (fits && len > 0)
            flint_mpn_copyi(z, buf, len);

        if (pool != NULL)
            fft_small_op_pool_scratch_give(pool, buf);
        else
            flint_free(buf);
    }

    if (fits && len < zn)
        flint_mpn_zero(z + len, zn - len);

    return fits;
}
