/*
    Copyright (C) 2022 Daniel Schultz
    Copyright (C) 2023, 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "mpn_mod.h"

#if FLINT_HAVE_FFT_SMALL

#include "thread_pool.h"
#include "thread_support.h"
#include "nmod.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "crt_helpers.h"
#include "fft_small.h"

static mp_limb_t
nmod_set_mpn_2(mp_srcptr ad, nmod_t mod)
{
    mp_limb_t r = 0;
    NMOD_RED2(r, r, ad[1], mod);
    NMOD_RED2(r, r, ad[0], mod);
    return r;
}

static mp_limb_t
nmod_set_mpn_3(mp_srcptr ad, nmod_t mod)
{
    mp_limb_t r = 0;
    NMOD_RED2(r, r, ad[2], mod);
    NMOD_RED2(r, r, ad[1], mod);
    NMOD_RED2(r, r, ad[0], mod);
    return r;
}

/* todo: precomputed inverse */
static mp_limb_t
nmod_set_mpn(mp_srcptr ad, mp_size_t an, nmod_t mod)
{
    return mpn_mod_1(ad, an, mod.n);
}

static void _mod(
    double* abuf, ulong atrunc,
    mp_srcptr a, ulong an,
    mp_size_t nlimbs,
    const sd_fft_ctx_struct* fft)
{
    double* aI;
    ulong i, j;
    nmod_t mod = fft->mod;

    if (atrunc < an)
    {
        flint_throw(FLINT_ERROR, "fft _mod: atrunc < an not handled\n");
    }

    for (i = 0; i + BLK_SZ <= an; i += BLK_SZ)
    {
        aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
        for (j = 0; j < BLK_SZ; j += 8)
        {
            ulong aa[8];

            FLINT_ASSERT(i+j < atrunc);

            /* todo: vectorize */
            if (nlimbs == 2)
            {
                aa[0] = nmod_set_mpn_2(a + (i + j + 0) * nlimbs, mod);
                aa[1] = nmod_set_mpn_2(a + (i + j + 1) * nlimbs, mod);
                aa[2] = nmod_set_mpn_2(a + (i + j + 2) * nlimbs, mod);
                aa[3] = nmod_set_mpn_2(a + (i + j + 3) * nlimbs, mod);
                aa[4] = nmod_set_mpn_2(a + (i + j + 4) * nlimbs, mod);
                aa[5] = nmod_set_mpn_2(a + (i + j + 5) * nlimbs, mod);
                aa[6] = nmod_set_mpn_2(a + (i + j + 6) * nlimbs, mod);
                aa[7] = nmod_set_mpn_2(a + (i + j + 7) * nlimbs, mod);
            }
            else if (nlimbs == 3)
            {
                aa[0] = nmod_set_mpn_3(a + (i + j + 0) * nlimbs, mod);
                aa[1] = nmod_set_mpn_3(a + (i + j + 1) * nlimbs, mod);
                aa[2] = nmod_set_mpn_3(a + (i + j + 2) * nlimbs, mod);
                aa[3] = nmod_set_mpn_3(a + (i + j + 3) * nlimbs, mod);
                aa[4] = nmod_set_mpn_3(a + (i + j + 4) * nlimbs, mod);
                aa[5] = nmod_set_mpn_3(a + (i + j + 5) * nlimbs, mod);
                aa[6] = nmod_set_mpn_3(a + (i + j + 6) * nlimbs, mod);
                aa[7] = nmod_set_mpn_3(a + (i + j + 7) * nlimbs, mod);
            }
            else
            {
                aa[0] = nmod_set_mpn(a + (i + j + 0) * nlimbs, nlimbs, mod);
                aa[1] = nmod_set_mpn(a + (i + j + 1) * nlimbs, nlimbs, mod);
                aa[2] = nmod_set_mpn(a + (i + j + 2) * nlimbs, nlimbs, mod);
                aa[3] = nmod_set_mpn(a + (i + j + 3) * nlimbs, nlimbs, mod);
                aa[4] = nmod_set_mpn(a + (i + j + 4) * nlimbs, nlimbs, mod);
                aa[5] = nmod_set_mpn(a + (i + j + 5) * nlimbs, nlimbs, mod);
                aa[6] = nmod_set_mpn(a + (i + j + 6) * nlimbs, nlimbs, mod);
                aa[7] = nmod_set_mpn(a + (i + j + 7) * nlimbs, nlimbs, mod);
            }

            vec8n t = vec8n_load_unaligned(aa);
            vec8d_store_aligned(aI + j, vec8n_convert_limited_vec8d(t));
        }
    }

    aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
    for (j = 0; j < an - i; j++)
        aI[j] = (double) nmod_set_mpn(a + (i + j) * nlimbs, nlimbs, mod);

    for (i = an; i < atrunc; i++)
        sd_fft_ctx_set_index(abuf, i, 0);
}

#define DEFINE_IT(NP, N, M) \
static void CAT(_crt, NP)( \
    mp_ptr z, ulong zl, ulong zi_start, ulong zi_stop, \
    sd_fft_ctx_struct* Rffts, double* d, ulong dstride, \
    crt_data_struct* Rcrts, mp_size_t nlimbs, gr_ctx_t ctx) \
{ \
    ulong np = NP; \
    ulong n = N; \
    ulong m = M; \
 \
    FLINT_ASSERT(n == Rcrts[np-1].coeff_len); \
    FLINT_ASSERT(1 <= N && N <= 7); \
 \
    if (n == m + 1) \
    { \
        for (ulong l = 0; l < np; l++) { \
            FLINT_ASSERT(crt_data_co_prime(Rcrts + np - 1, l)[m] == 0); \
        } \
    } \
    else \
    { \
        FLINT_ASSERT(n == m); \
    } \
 \
    ulong Mhalf[N]; \
    mpn_rshift(Mhalf, crt_data_prod_primes(Rcrts + np - 1), N, 1); \
 \
    ulong Xs[BLK_SZ*NP]; \
 \
    for (ulong i = n_round_down(zi_start, BLK_SZ); i < zi_stop; i += BLK_SZ) \
    { \
        _convert_block(Xs, Rffts, d, dstride, np, i/BLK_SZ); \
 \
        ulong jstart = (i < zi_start) ? zi_start - i : 0; \
        ulong jstop = FLINT_MIN(BLK_SZ, zi_stop - i); \
        for (ulong j = jstart; j < jstop; j += 1) \
        { \
            ulong r[N]; \
            ulong t[N]; \
            ulong l = 0; \
            ulong nn; \
 \
            CAT3(_big_mul, N, M)(r, t, _crt_data_co_prime(Rcrts + np - 1, l, n), Xs[l*BLK_SZ + j]); \
            for (l++; l < np; l++) \
                CAT3(_big_addmul, N, M)(r, t, _crt_data_co_prime(Rcrts + np - 1, l, n), Xs[l*BLK_SZ + j]); \
 \
            CAT(_reduce_big_sum, N)(r, t, crt_data_prod_primes(Rcrts + np - 1)); \
 \
            FLINT_ASSERT(mpn_cmp(r, Mhalf, N) <= 0); /* coefficients must be unsigned */ \
            nn = N; \
            MPN_NORM(r, nn); \
            mpn_mod_set_mpn(z + (i + j - zl) * nlimbs, r, nn, ctx); \
        } \
    } \
}

DEFINE_IT(2, 2, 1)  /* 100 bits (unused) */
DEFINE_IT(3, 3, 2)  /* 150 bits */
DEFINE_IT(4, 4, 3)  /* 200 bits */
DEFINE_IT(5, 4, 4)  /* 250 bits */
DEFINE_IT(6, 5, 4)  /* 300 bits */
DEFINE_IT(7, 6, 5)  /* 350 bits */
DEFINE_IT(8, 7, 6)  /* 400 bits */
#undef DEFINE_IT

/* 50 bits (unused) */
static void _crt_1(
    mp_ptr FLINT_UNUSED(z), ulong FLINT_UNUSED(zl), ulong FLINT_UNUSED(zi_start), ulong FLINT_UNUSED(zi_stop),
    sd_fft_ctx_struct* FLINT_UNUSED(Rffts), double* FLINT_UNUSED(d), ulong FLINT_UNUSED(dstride),
    crt_data_struct* FLINT_UNUSED(Rcrts), mp_size_t FLINT_UNUSED(nlimbs), gr_ctx_t FLINT_UNUSED(ctx))
{
    flint_abort();
}

typedef struct {
    ulong np;
    ulong start_pi;
    ulong stop_pi;
    ulong offset;
    double* abuf;
    double* bbuf;
    ulong depth;
    ulong stride;
    ulong atrunc;
    ulong btrunc;
    ulong ztrunc;
    mp_srcptr a;
    ulong an;
    mp_srcptr b;
    ulong bn;
    mp_size_t nlimbs;
    gr_ctx_struct * mpn_mod_ctx;
    sd_fft_ctx_struct* ffts;
    crt_data_struct* crts;
    ulong ioff;
    int squaring;
    int want_worker;
} s1worker_struct;


static void extra_func(void* varg)
{
    s1worker_struct* X = (s1worker_struct*) varg;
    sd_fft_lctx_t Q;

    sd_fft_lctx_init(Q, X->ffts + X->ioff, X->depth);
    _mod(X->bbuf, X->btrunc, X->b, X->bn, X->nlimbs, X->ffts + X->ioff);
    sd_fft_lctx_fft_trunc(Q, X->bbuf, X->depth, X->btrunc, X->ztrunc);
    sd_fft_lctx_clear(Q, X->ffts + X->ioff);
}

static void s1worker_func(void* varg)
{
    s1worker_struct* X = (s1worker_struct*) varg;
    sd_fft_lctx_t Q;
    ulong i, m;
    thread_pool_handle* handles = NULL;
    slong nworkers = 0;

    if (X->want_worker)
        nworkers = flint_request_threads(&handles, 2);

    for (i = X->start_pi; i < X->stop_pi; i++)
    {
        ulong ioff = i + X->offset;
        double* abuf = X->abuf + X->stride*i;
        double* bbuf = X->bbuf;

        sd_fft_lctx_init(Q, X->ffts + ioff, X->depth);

        if (!X->squaring)
        {
            if (nworkers > 0)
            {
                X->ioff = ioff;
                thread_pool_wake(global_thread_pool, handles[0], 0, extra_func, X);
            }
            else
            {
                _mod(bbuf, X->btrunc, X->b, X->bn, X->nlimbs, X->ffts + ioff);
                sd_fft_lctx_fft_trunc(Q, bbuf, X->depth, X->btrunc, X->ztrunc);
            }
        }

        _mod(abuf, X->atrunc, X->a, X->an, X->nlimbs, X->ffts + ioff);
        sd_fft_lctx_fft_trunc(Q, abuf, X->depth, X->atrunc, X->ztrunc);

        if (!X->squaring)
        {
            if (nworkers > 0)
                thread_pool_wait(global_thread_pool, handles[0]);
        }

        ulong cop = X->np == 1 ? 1 : *crt_data_co_prime_red(X->crts + X->np - 1, ioff);
        NMOD_RED2(m, cop >> (FLINT_BITS - X->depth), cop << X->depth, X->ffts[ioff].mod);
        m = nmod_inv(m, X->ffts[ioff].mod);

        if (X->squaring)
            sd_fft_lctx_point_sqr(Q, abuf, m, X->depth);
        else
            sd_fft_lctx_point_mul(Q, abuf, bbuf, m, X->depth);

        sd_fft_lctx_ifft_trunc(Q, abuf, X->depth, X->ztrunc);

        sd_fft_lctx_clear(Q, X->ffts + ioff);
    }

    flint_give_back_threads(handles, nworkers);
}

typedef struct {
    mp_ptr z;
    ulong zl;
    ulong start_zi;
    ulong stop_zi;
    double* buf;
    ulong offset;
    ulong stride;
    sd_fft_ctx_struct* ffts;
    crt_data_struct* crts;
    nmod_t mod;
    mp_size_t nlimbs;
    gr_ctx_struct * mpn_mod_ctx;
    void (*f)(
        mp_ptr z, ulong zl, ulong zi_start, ulong zi_stop,
        sd_fft_ctx_struct* Rffts, double* d, ulong dstride,
        crt_data_struct* Rcrts, mp_size_t nlimbs, gr_ctx_t ctx);
} s2worker_struct;

static void s2worker_func(void* varg)
{
    s2worker_struct* X = (s2worker_struct*) varg;

    X->f(X->z, X->zl, X->start_zi, X->stop_zi, X->ffts + X->offset, X->buf,
         X->stride, X->crts + X->offset, X->nlimbs, X->mpn_mod_ctx);
}

int _mpn_mod_poly_mulmid_fft_small_internal(mp_ptr z, ulong zl, ulong zh,
    mp_srcptr a, ulong an,
    mp_srcptr b, ulong bn,
    mpn_ctx_t R, gr_ctx_t ctx)
{
    ulong modbits;
    ulong offset = 0;
    ulong zn = an + bn - 1;
    ulong atrunc, btrunc, ztrunc;
    ulong i, np, depth, stride;
    double* buf;
    int squaring;
    slong bits1, bits2;
    int sign = 0;
    mp_size_t nlimbs = MPN_MOD_CTX_NLIMBS(ctx);
    mp_bitcnt_t nbits = MPN_MOD_CTX_MODULUS_BITS(ctx);

    FLINT_ASSERT(an > 0);
    FLINT_ASSERT(bn > 0);

    if (zl >= zh)
        return GR_SUCCESS;

    if (zh > zn)
    {
        if (zl >= zn)
        {
            _mpn_mod_vec_zero(z, zh - zl, ctx);
            return GR_SUCCESS;
        }

        _mpn_mod_vec_zero(z + (zn - zl) * nlimbs, zh - zn, ctx);
        zh = zn;
    }

    squaring = (a == b) && (an == bn);

    /* TODO: consider counting bits, in case the inputs are small */
    bits1 = bits2 = nbits;

    (void) sign;
    /* should be +sign instead of +1, but currently the CRT code doesn't
       distinguish between the signed and unsigned cases */
    modbits = FLINT_ABS(bits1) + FLINT_ABS(bits2) + 1;

    FLINT_ASSERT(zl < zh);
    FLINT_ASSERT(zh <= zn);

    /* need prod_of_primes >= blen * 2^modbits */
    for (np = 1; ; np++)
    {
        if (np > MPN_CTX_NCRTS)
            return GR_UNABLE;

        if (flint_mpn_cmp_ui_2exp(crt_data_prod_primes(R->crts + np - 1),
              R->crts[np - 1].coeff_len, bn, modbits) >= 0)
        {
            break;
        }
    }

    FLINT_ASSERT(0 <= flint_mpn_cmp_ui_2exp(
                                  crt_data_prod_primes(R->crts + np - 1),
                                  R->crts[np - 1].coeff_len, bn, modbits));

    atrunc = n_round_up(an, BLK_SZ);
    btrunc = n_round_up(bn, BLK_SZ);
    ztrunc = n_round_up(zn, BLK_SZ);
    /*
        if there is a power of two 2^d between zh and zn with good wrap around
            i.e. max(an, bn, zh) <= 2^d <= zn with zn - 2^d <= zl
        then use d as the depth, otherwise the usual with no wrap around
    */
    depth = n_flog2(zn);
    i = n_pow2(depth);
    if (atrunc <= i && btrunc <= i && zh <= i && i <= zn && zn <= zl + i)
    {
        ztrunc = i;
    }
    else
    {
        depth = n_max(LG_BLK_SZ, n_clog2(ztrunc));
    }

    stride = n_round_up(sd_fft_ctx_data_size(depth), 128);

    ulong want_threads;

    if ((np >= 2 && bn >= 1000) || (np >= 4 && bn >= 300))
        want_threads = np;
    else
        want_threads = 1;

    thread_pool_handle* handles;
    slong nworkers = flint_request_threads(&handles, want_threads);
    ulong nthreads = nworkers + 1;

    buf = (double*) mpn_ctx_fit_buffer(R, (np+nthreads)*stride*sizeof(double));

    s1worker_struct s1args[8];
    FLINT_ASSERT(nthreads <= 8);
    for (i = 0; i < nthreads; i++)
    {
        s1worker_struct* X = s1args + i;
        X->np = np;
        X->start_pi = (i+0)*np/nthreads;
        X->stop_pi  = (i+1)*np/nthreads;
        X->offset = offset;
        X->abuf = buf;
        X->bbuf = buf + (np+i)*stride;
        X->depth = depth;
        X->stride = stride;
        X->atrunc = atrunc;
        X->btrunc = btrunc;
        X->ztrunc = ztrunc;
        X->a = a;
        X->an = an;
        X->b = b;
        X->bn = bn;
        X->nlimbs = nlimbs;
        X->mpn_mod_ctx = ctx;
        X->ffts = R->ffts;
        X->crts = R->crts;
        X->squaring = squaring;
        X->want_worker = !squaring && ((np == 1 && bn > 5000) || (np >= 2 && bn >= 1000));
    }

    for (i = nworkers; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0, s1worker_func, s1args + i);
    s1worker_func(s1args + 0);
    for (i = nworkers; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    if (zn > 50000 || (np >= 2 && zn > 20000) || (np >= 4 && zn > 800))
    {
        flint_give_back_threads(handles, nworkers);
        nworkers = flint_request_threads(&handles, 8);
        nthreads = nworkers + 1;
    }

    s2worker_struct s2args[8];
    FLINT_ASSERT(nthreads <= 8);

    ulong o = zl;
    for (i = 0; i < nthreads; i++)
    {
        s2worker_struct* X = s2args + i;
        X->z = z;
        X->zl = zl;
        X->start_zi = o;
        ulong newo = n_round_down(zl + (i+1)*(zh-zl)/nthreads, BLK_SZ);
        o = i+1 < nthreads ? FLINT_MAX(o, newo) : zh;
        X->stop_zi = o;
        X->buf = buf;
        X->offset = offset;
        X->stride = stride;
        X->ffts = R->ffts;
        X->crts = R->crts;
        X->nlimbs = nlimbs;
        X->mpn_mod_ctx = ctx;
        X->f =  np == 1 ? _crt_1 :
                np == 2 ? _crt_2 :
                np == 3 ? _crt_3 :
                np == 4 ? _crt_4 :
                np == 5 ? _crt_5 :
                np == 6 ? _crt_6 :
                np == 7 ? _crt_7 :
                          _crt_8;
    }

    for (i = nworkers; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0, s2worker_func, s2args + i);
    s2worker_func(s2args + 0);
    for (i = nworkers; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    flint_give_back_threads(handles, nworkers);

    return GR_SUCCESS;
}

int
_mpn_mod_poly_mullow_fft_small(mp_ptr res, mp_srcptr poly1, slong len1, mp_srcptr poly2, slong len2, slong len, gr_ctx_t ctx)
{
    if (len1 >= len2)
        return _mpn_mod_poly_mulmid_fft_small_internal(res, 0, len, poly1, len1, poly2, len2, get_default_mpn_ctx(), ctx);
    else
        return _mpn_mod_poly_mulmid_fft_small_internal(res, 0, len, poly2, len2, poly1, len1, get_default_mpn_ctx(), ctx);
}

#else /* FLINT_HAVE_FFT_SMALL */

int
_mpn_mod_poly_mullow_fft_small(mp_ptr res,  mp_srcptr poly1, slong len1, mp_srcptr poly2, slong len2, slong len, gr_ctx_t ctx)
{
    return GR_UNABLE;
}

#endif
