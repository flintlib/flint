/*
    Copyright (C) 2026 Fredrik Johansson

    Developed using Claude Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "ulong_extras.h"
#include "mpn_extras.h"
#include "fft_small.h"
#include "machine_vectors.h"
#include "thread_pool.h"
#include "thread_support.h"
#include <string.h>

/*
    Two products with a common operand, z1 = a x and z2 = b x, with the
    transform of x computed once: three forward and two inverse
    transforms against four and two for two independent products.

    Two-prime geometry over the range where the multiplication driver
    uses it (operands of at most FLINT_MPN_MUL_NP2_MAX_AN limbs), the
    wide geometry beyond (below). As in _mpn_ctx_mpn_mul_window_np2,
    each product is computed in one pass while its data is hot -- pack,
    forward transforms, pointwise product, inverse transforms, Garner
    recomposition -- over four lanes of scratch (two holding the
    transform of x, two for the current operand). The transforms of x
    run at the depth and truncation of the longer product; the other
    operand of each product, its pointwise product and its inverse
    transform are truncated to that product's own length, so a shorter
    b x costs less than a x in those stages.

    In the two-prime range a version on the transformed-operand API
    (fft_small_fft_mpn, fft_small_op_mul, fft_small_ifft_trunc), making
    the three forward transforms, the pointwise products and the
    inverse transforms separate passes over all three operands, was as
    slow as two independent products: the extra traffic costs what the
    shared transform saves.

    Like the multiplication driver, these functions have a plan for all
    lengths: two-prime when its exactness bound holds (always within
    its range), wide otherwise, which needs only that eight primes hold
    2 bits + log2(length) bits for chunks of just over a limb (lengths
    beyond 2^200 limbs).

    Exactness: at most min(max(alen, blen), xlen) chunk products
    accumulate on one coefficient, the bound the chunk size is chosen
    for, and each product is below B^(an + xn) (resp. B^(bn + xn)), so
    the recomposition of that many limbs is exact.
*/

static void
_pointwise(const sd_fft_ctx_struct * Q, double * X, const double * Y,
           ulong n)
{
    vec8d n8 = vec8d_set_d(Q->p), ninv8 = vec8d_set_d(Q->pinv);
    ulong j;

    for (j = 0; j < n; j += 8)
        vec8d_store(X + j, vec8d_mulmod(vec8d_load(X + j),
                                        vec8d_load(Y + j), n8, ninv8));
}

/* the plan of the products of x (xn limbs) by operands of up to yn_max
   limbs: the depth for the longest product, the exactness bound for
   min(yn_max, xn) accumulated chunk products; 0 above the two-prime
   range */
static int
_np2_plan(mpn_ctx_struct * R, ulong * bits, ulong * depth, ulong xn,
          ulong yn_max)
{
    if (FLINT_MAX(yn_max, xn) > FLINT_MPN_MUL_NP2_MAX_AN)
        return 0;

    return _fft_small_np2_choose(R, FLINT_MAX(yn_max, xn),
                                 FLINT_MIN(yn_max, xn), 1, bits, depth);
}

/* the two lanes x0, x1 of the transform of x, truncated to the
   longest product (by an operand of yn_max limbs) */
static void
_np2_transform(mpn_ctx_struct * R, double * x0, double * x1, nn_srcptr x,
               ulong xn, ulong yn_max, ulong bits, ulong depth)
{
    ulong xlen = n_cdiv(FLINT_BITS * xn, bits);
    ulong xtr = n_round_up(xlen, BLK_SZ);
    ulong ztr = n_round_up(n_cdiv(FLINT_BITS * yn_max, bits) + xlen - 1,
                           BLK_SZ);

    _fft_small_np2_pack(x0, x1, x, xn, xlen, xtr, bits);
    sd_fft_trunc(R->ffts + 0, x0, depth, xtr, ztr);
    sd_fft_trunc(R->ffts + 1, x1, depth, xtr, ztr);
}

/* z = x y (xn + yn limbs) from the lanes of x, in one pass while the
   data is hot: pack, forward transforms, pointwise products, inverse
   transforms, Garner recomposition, over the scratch lanes d0, d1; z
   may overlap y */
static void
_np2_mul_lanes(mpn_ctx_struct * R, nn_ptr z, nn_srcptr y, ulong yn,
               const double * x0, const double * x1, ulong xn, ulong bits,
               ulong depth, double * d0, double * d1)
{
    sd_fft_ctx_struct * Q0 = R->ffts + 0, * Q1 = R->ffts + 1;
    ulong xlen = n_cdiv(FLINT_BITS * xn, bits);
    ulong ylen = n_cdiv(FLINT_BITS * yn, bits);
    ulong ytr = n_round_up(ylen, BLK_SZ);
    ulong zlen = ylen + xlen - 1;
    ulong ztr = n_round_up(zlen, BLK_SZ);

    _fft_small_np2_pack(d0, d1, y, yn, ylen, ytr, bits);
    sd_fft_trunc(Q0, d0, depth, ytr, ztr);
    _pointwise(Q0, d0, x0, ztr);
    sd_ifft_trunc(Q0, d0, depth, ztr);
    sd_fft_trunc(Q1, d1, depth, ytr, ztr);
    _pointwise(Q1, d1, x1, ztr);
    sd_ifft_trunc(Q1, d1, depth, ztr);

    _fft_small_np2_crt_recompose(R, z, 0, yn + xn, d0, d1, zlen,
                                 ztr, bits, depth);
}

/*
    Products by a fixed operand with its transform computed once, for
    an operand used in many products (the power tables of binary
    splitting): the plan is fixed at init for partners of up to yn_max
    limbs, and each product costs the conversion and forward
    transforms of the partner, the pointwise products, the inverse
    transforms and the reconstruction -- one conversion and forward
    transform less than a product.

    Two-prime geometry up to FLINT_MPN_MUL_NP2_MAX_AN limbs, as the
    multiplication driver, with the per-product pipeline above. Above
    it, a wide plan: the partner is converted and transformed on the
    transformed-operand API (threaded over the primes), its pointwise
    product and inverse transform run per prime while the prime's data
    is hot (threaded over the primes), and it is reconstructed by the
    threaded export; the partner's transform lives in the calling
    thread's context buffer.

    Measured against flint_mpn_mul for partners of 1.75 times the
    length of x: two-prime 26-31% faster from 300 to 28000 limbs; wide
    16-22% faster from 43000 to 920000 limbs single-threaded.

    The transform of x is read only, so the products may run on several
    threads, each with its own context's scratch: the default contexts
    have the same primes, and the transforms extend their own tables.
*/

/* the points of the transform of x that the products read: those of
   the longest product */
static ulong
_np2_ztr(ulong xn, ulong yn_max, ulong bits)
{
    return n_round_up(n_cdiv(FLINT_BITS * yn_max, bits)
                      + n_cdiv(FLINT_BITS * xn, bits) - 1, BLK_SZ);
}

/* the wide plan of the products of x by operands of up to yn_max
   limbs: the automatic wide selection (np >= 3, as the multiplication
   driver beyond the two-prime range; forced below it, where the
   two-prime plan normally serves) */
static void
_wide_plan(fft_small_plan_t P, ulong xn, ulong yn_max)
{
    if (!fft_small_plan_init_mpn_np(P, get_default_mpn_ctx(), yn_max, xn,
                                    1, 0, 3))
        flint_throw(FLINT_ERROR, "fft_small: no plan for products of "
                    "%wu by %wu limbs\n", xn, yn_max);

    FLINT_ASSERT(P->np >= 3);
}

ulong
_fft_small_mpn_mul_precomp_size(ulong xn, ulong yn_max)
{
    ulong bits, depth, size;
    fft_small_plan_t P;

    FLINT_ASSERT(xn > 0 && yn_max > 0);

    if (_np2_plan(get_default_mpn_ctx(), &bits, &depth, xn, yn_max))
        return 2 * _np2_ztr(xn, yn_max, bits);

    _wide_plan(P, xn, yn_max);
    size = P->np * P->ztrunc;
    fft_small_plan_clear(P);
    return size;
}

void
_fft_small_mpn_mul_precomp_init(fft_small_mpn_mul_precomp_t P,
                                nn_srcptr x, ulong xn, ulong yn_max,
                                double * storage)
{
    mpn_ctx_struct * R = get_default_mpn_ctx();

    FLINT_ASSERT(xn > 0 && yn_max > 0);

    P->wide = 0;
    P->xn = xn;
    P->yn_max = yn_max;
    P->owns = (storage == NULL);

    /* the transform is made in the context's scratch at its full
       depth and only the points the products read are kept: the
       storage is up to half that of the full transform */
    if (_np2_plan(R, &P->bits, &P->depth, xn, yn_max))
    {
        ulong dsz = sd_fft_ctx_data_size(P->depth);
        double * t = mpn_ctx_fit_buffer(R, 2 * dsz * sizeof(double));

        P->xstride = _np2_ztr(xn, yn_max, P->bits);
        P->d = (storage != NULL) ? storage
                 : flint_aligned_alloc(64, 2 * P->xstride * sizeof(double));
        _np2_transform(R, t, t + dsz, x, xn, yn_max, P->bits, P->depth);
        memcpy(P->d, t, P->xstride * sizeof(double));
        memcpy(P->d + P->xstride, t + dsz, P->xstride * sizeof(double));
    }
    else
    {
        fft_small_op_t X;
        ulong i;

        _wide_plan(P->P, xn, yn_max);
        P->wide = 1;
        P->xstride = P->P->ztrunc;
        P->d = (storage != NULL) ? storage
                 : flint_aligned_alloc(64, P->P->np * P->xstride
                                           * sizeof(double));
        fft_small_op_init_borrowed(X, P->P, mpn_ctx_fit_buffer(R,
                                       fft_small_op_sizeof_data(P->P)));
        fft_small_fft_mpn(X, x, xn, P->P);
        for (i = 0; i < P->P->np; i++)
            memcpy(P->d + i * P->xstride, X->data + i * P->P->stride,
                   P->xstride * sizeof(double));
    }
}

void
_fft_small_mpn_mul_precomp_clear(fft_small_mpn_mul_precomp_t P)
{
    if (P->d != NULL)
    {
        if (P->owns)
            flint_aligned_free(P->d);
        if (P->wide)
            fft_small_plan_clear(P->P);
    }
    P->d = NULL;
}

/* the pointwise product by the transform of x, with the normalization
   factor folded in, then the inverse transform, per prime */
typedef struct
{
    const fft_small_plan_struct * P;
    double * y;
    const double * x;
    ulong xstride;
    ulong itrunc;
    ulong clear_to;
    ulong start_pi;
    ulong stop_pi;
}
_wide_worker_struct;

static void
_wide_worker(void * varg)
{
    _wide_worker_struct * W = (_wide_worker_struct *) varg;
    const fft_small_plan_struct * P = W->P;
    ulong i, j;

    for (i = W->start_pi; i < W->stop_pi; i++)
    {
        sd_fft_ctx_struct * Q = P->ffts + P->offset + i;
        double * d = W->y + P->stride * i;
        const double * e = W->x + W->xstride * i;
        vec8d p = vec8d_set_d(Q->p), pinv = vec8d_set_d(Q->pinv);
        vec8d m = vec8d_set_d(vec1d_reduce_0n_to_pmhn((slong) P->m[i], Q->p));

        for (j = 0; j < W->itrunc; j += 8)
            vec8d_store(d + j, vec8d_mulmod(vec8d_mulmod(vec8d_load(d + j),
                        m, p, pinv), vec8d_load(e + j), p, pinv));

        sd_ifft_trunc(Q, d, P->depth, W->itrunc);

        if (W->clear_to > W->itrunc)
            memset(d + W->itrunc, 0, (W->clear_to - W->itrunc) * sizeof(double));
    }
}

/* z = x y (xn + yn limbs) on the wide plan P from the transform of x
   (np lanes of xstride points), with the transform of y in ydata (the
   data of an operand of P); z may overlap y */
static void
_wide_mul(nn_ptr z, const fft_small_plan_struct * P, const double * xdata,
          ulong xstride, ulong xn, nn_srcptr y, ulong yn, double * ydata)
{
    ulong zn = yn + xn, nslots, need, en, itrunc, clear_to, i, nthreads;
    thread_pool_handle * handles = NULL;
    slong nworkers = 0;
    _wide_worker_struct args[MPN_CTX_NCRTS];
    fft_small_op_t Y;

    nslots = fft_small_plan_mpn_nslots(P, yn, xn);
    need = fft_small_export_mpn_limbs(P, 1, nslots);
    en = FLINT_MIN(zn, need);
    itrunc = n_min(n_round_up(nslots, BLK_SZ), n_pow2(P->depth));
    clear_to = n_min(n_cdiv(en * FLINT_BITS, P->bits), n_pow2(P->depth));

    /* the transform of y only as long as this product: a partner
       shorter than the plan's longest product costs only its own
       length (forming it at the full truncation made a pair with one
       product half as long slower than two products) */
    fft_small_op_init_borrowed(Y, P, ydata);
    _fft_small_fft_mpn_trunc(Y, y, yn, n_min(itrunc, P->ztrunc), P);

    if (itrunc >= 3072)
        nworkers = flint_request_threads(&handles, P->np);
    nthreads = FLINT_MIN((ulong) (nworkers + 1), P->np);

    for (i = 0; i < nthreads; i++)
    {
        args[i].P = P;
        args[i].y = Y->data;
        args[i].x = xdata;
        args[i].xstride = xstride;
        args[i].itrunc = itrunc;
        args[i].clear_to = clear_to;
        args[i].start_pi = (i + 0) * P->np / nthreads;
        args[i].stop_pi = (i + 1) * P->np / nthreads;
    }

    for (i = nthreads - 1; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0,
                         _wide_worker, args + i);
    _wide_worker(args + 0);
    for (i = nthreads - 1; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    flint_give_back_threads(handles, nworkers);

    Y->domain = FFT_SMALL_OP_PRODUCT;
    fft_small_export_mpn(z, en, Y, P);
    if (zn > en)
        flint_mpn_zero(z + en, zn - en);
}

void
_fft_small_mpn_mul_precomp(nn_ptr z, const fft_small_mpn_mul_precomp_t P,
                           nn_srcptr y, ulong yn)
{
    mpn_ctx_struct * R = get_default_mpn_ctx();
    ulong dsz;
    double * d;

    FLINT_ASSERT(P->d != NULL);
    FLINT_ASSERT(yn > 0 && yn <= P->yn_max);

    if (P->wide)
    {
        _wide_mul(z, P->P, P->d, P->xstride, P->xn, y, yn,
                  mpn_ctx_fit_buffer(R, fft_small_op_sizeof_data(P->P)));
        return;
    }

    dsz = sd_fft_ctx_data_size(P->depth);
    d = mpn_ctx_fit_buffer(R, 2 * dsz * sizeof(double));
    _np2_mul_lanes(R, z, y, yn, P->d, P->d + P->xstride, P->xn, P->bits,
                   P->depth, d, d + dsz);
}

/* the shared products in the wide geometry: the transform of x and
   that of the current partner in one region of the context buffer */
static void
_mul_shared_wide(nn_ptr z1, nn_srcptr a, ulong an, nn_ptr z2, nn_srcptr b,
                 ulong bn, nn_srcptr x, ulong xn)
{
    fft_small_plan_t P;
    fft_small_op_t X;
    ulong size;
    double * buf;

    _wide_plan(P, xn, FLINT_MAX(an, bn));

    size = fft_small_op_sizeof_data(P);
    buf = mpn_ctx_fit_buffer(get_default_mpn_ctx(), 2 * size);
    fft_small_op_init_borrowed(X, P, buf);
    fft_small_fft_mpn(X, x, xn, P);

    /* z1 is written before b is read */
    _wide_mul(z1, P, X->data, P->stride, xn, a, an,
              buf + size / sizeof(double));
    _wide_mul(z2, P, X->data, P->stride, xn, b, bn,
              buf + size / sizeof(double));

    fft_small_plan_clear(P);
}

void
_fft_small_mpn_mul_shared(nn_ptr z1, nn_srcptr a, ulong an,
                          nn_ptr z2, nn_srcptr b, ulong bn,
                          nn_srcptr x, ulong xn)
{
    mpn_ctx_struct * R = get_default_mpn_ctx();
    ulong m = FLINT_MAX(an, bn), bits, depth, dsz;
    double * x0;

    FLINT_ASSERT(an > 0 && bn > 0 && xn > 0);

    if (!_np2_plan(R, &bits, &depth, xn, m))
    {
        _mul_shared_wide(z1, a, an, z2, b, bn, x, xn);
        return;
    }

    dsz = sd_fft_ctx_data_size(depth);
    x0 = mpn_ctx_fit_buffer(R, 4 * dsz * sizeof(double));

    _np2_transform(R, x0, x0 + dsz, x, xn, m, bits, depth);
    /* z1 is written before b is read */
    _np2_mul_lanes(R, z1, a, an, x0, x0 + dsz, xn, bits, depth,
                   x0 + 2 * dsz, x0 + 3 * dsz);
    _np2_mul_lanes(R, z2, b, bn, x0, x0 + dsz, xn, bits, depth,
                   x0 + 2 * dsz, x0 + 3 * dsz);
}
