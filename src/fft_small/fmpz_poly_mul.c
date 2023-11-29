/*
    Copyright (C) 2022 Daniel Schultz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "fft_small.h"
#include "crt_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

static void _mod(
    double* abuf, ulong atrunc,
    const fmpz * a, ulong an,
    slong abits,
    const sd_fft_ctx_struct* fft)
{
    double* aI;
    ulong i, j;
    nmod_t mod = fft->mod;
    ulong p = mod.n;

    if (atrunc < an)
    {
        flint_throw(FLINT_ERROR, "fft _mod: atrunc < an not handled\n");
    }

    if (FLINT_ABS(abits) < FLINT_BIT_COUNT(p))
    {
        if (abits >= 0)
        {
            for (i = 0; i + BLK_SZ <= an; i += BLK_SZ)
            {
                aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
                for (j = 0; j < BLK_SZ; j += 8)
                {
                    vec8n t = vec8n_load_unaligned((ulong *) (a + i + j));
                    FLINT_ASSERT(i+j < atrunc);
                    vec8d_store_aligned(aI + j, vec8n_convert_limited_vec8d(t));
                }
            }

            aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
            for (j = 0; j < an - i; j++)
                aI[j] = a[i + j];
        }
        else
        {
            for (i = 0; i + BLK_SZ <= an; i += BLK_SZ)
            {
                aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
                for (j = 0; j < BLK_SZ; j += 8)
                {
                    ulong aa[8];

                    /* todo: vectorize */
                    aa[0] = a[i + j + 0] >= 0 ? a[i + j + 0] : a[i + j + 0] + p;
                    aa[1] = a[i + j + 1] >= 0 ? a[i + j + 1] : a[i + j + 1] + p;
                    aa[2] = a[i + j + 2] >= 0 ? a[i + j + 2] : a[i + j + 2] + p;
                    aa[3] = a[i + j + 3] >= 0 ? a[i + j + 3] : a[i + j + 3] + p;
                    aa[4] = a[i + j + 4] >= 0 ? a[i + j + 4] : a[i + j + 4] + p;
                    aa[5] = a[i + j + 5] >= 0 ? a[i + j + 5] : a[i + j + 5] + p;
                    aa[6] = a[i + j + 6] >= 0 ? a[i + j + 6] : a[i + j + 6] + p;
                    aa[7] = a[i + j + 7] >= 0 ? a[i + j + 7] : a[i + j + 7] + p;

                    vec8n t = vec8n_load_unaligned(aa);
                    FLINT_ASSERT(i+j < atrunc);
                    vec8d_store_aligned(aI + j, vec8n_convert_limited_vec8d(t));
                }
            }

            aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
            for (j = 0; j < an - i; j++)
                aI[j] = a[i + j] >= 0 ? a[i + j] : a[i + j] + p;
        }
    }
    else if (FLINT_ABS(abits) <= SMALL_FMPZ_BITCOUNT_MAX)
    {
        if (abits >= 0)
        {
            for (i = 0; i + BLK_SZ <= an; i += BLK_SZ)
            {
                aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
                for (j = 0; j < BLK_SZ; j += 8)
                {
                    ulong aa[8];

                    FLINT_ASSERT(i+j < atrunc);

                    /* todo: vectorize */
                    aa[0] = nmod_set_ui(a[i + j + 0], mod);
                    aa[1] = nmod_set_ui(a[i + j + 1], mod);
                    aa[2] = nmod_set_ui(a[i + j + 2], mod);
                    aa[3] = nmod_set_ui(a[i + j + 3], mod);
                    aa[4] = nmod_set_ui(a[i + j + 4], mod);
                    aa[5] = nmod_set_ui(a[i + j + 5], mod);
                    aa[6] = nmod_set_ui(a[i + j + 6], mod);
                    aa[7] = nmod_set_ui(a[i + j + 7], mod);

                    vec8n t = vec8n_load_unaligned(aa);
                    vec8d_store_aligned(aI + j, vec8n_convert_limited_vec8d(t));
                }
            }

            aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
            for (j = 0; j < an - i; j++)
                aI[j] = nmod_set_ui(a[i + j], mod);
        }
        else
        {
            for (i = 0; i + BLK_SZ <= an; i += BLK_SZ)
            {
                aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
                for (j = 0; j < BLK_SZ; j += 8)
                {
                    ulong aa[8];

                    FLINT_ASSERT(i+j < atrunc);

                    /* todo: vectorize */
                    aa[0] = nmod_set_si(a[i + j + 0], mod);
                    aa[1] = nmod_set_si(a[i + j + 1], mod);
                    aa[2] = nmod_set_si(a[i + j + 2], mod);
                    aa[3] = nmod_set_si(a[i + j + 3], mod);
                    aa[4] = nmod_set_si(a[i + j + 4], mod);
                    aa[5] = nmod_set_si(a[i + j + 5], mod);
                    aa[6] = nmod_set_si(a[i + j + 6], mod);
                    aa[7] = nmod_set_si(a[i + j + 7], mod);

                    vec8n t = vec8n_load_unaligned(aa);
                    vec8d_store_aligned(aI + j, vec8n_convert_limited_vec8d(t));
                }
            }

            aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
            for (j = 0; j < an - i; j++)
                aI[j] = nmod_set_si(a[i + j], mod);
        }
    }
    else
    {
        for (i = 0; i + BLK_SZ <= an; i += BLK_SZ)
        {
            aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
            for (j = 0; j < BLK_SZ; j += 8)
            {
                ulong aa[8];

                FLINT_ASSERT(i+j < atrunc);

                /* todo: vectorize? */
                aa[0] = fmpz_get_nmod(&a[i + j + 0], mod);
                aa[1] = fmpz_get_nmod(&a[i + j + 1], mod);
                aa[2] = fmpz_get_nmod(&a[i + j + 2], mod);
                aa[3] = fmpz_get_nmod(&a[i + j + 3], mod);
                aa[4] = fmpz_get_nmod(&a[i + j + 4], mod);
                aa[5] = fmpz_get_nmod(&a[i + j + 5], mod);
                aa[6] = fmpz_get_nmod(&a[i + j + 6], mod);
                aa[7] = fmpz_get_nmod(&a[i + j + 7], mod);

                vec8n t = vec8n_load_unaligned(aa);
                vec8d_store_aligned(aI + j, vec8n_convert_limited_vec8d(t));
            }
        }

        aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
        for (j = 0; j < an - i; j++)
            aI[j] = (double) fmpz_get_nmod(&a[i + j], mod);
    }

    for (i = an; i < atrunc; i++)
        sd_fft_ctx_set_index(abuf, i, 0);
}


void fmpz_neg_ui_array(fmpz_t out, const ulong * in, slong in_len)
{
    slong size = in_len;
    FLINT_ASSERT(in_len > 0);

    /* find end of zero extension */
    while (size > WORD(1) && in[size - 1] == UWORD(0))
        size--;

    /* copy limbs */
    if (size == WORD(1))
    {
        fmpz_neg_ui(out, in[0]);
    }
    else
    {
        __mpz_struct * mpz = _fmpz_promote(out);
        if (mpz->_mp_alloc < size)
            mpz_realloc2(mpz, FLINT_BITS * size);
        mpz->_mp_size = -size;
        flint_mpn_copyi(mpz->_mp_d, in, size);
    }
}

#define DEFINE_IT(NP, N, M) \
static void CAT(_crt, NP)( \
    fmpz * z, ulong zl, ulong zi_start, ulong zi_stop, \
    sd_fft_ctx_struct* Rffts, double* d, ulong dstride, \
    crt_data_struct* Rcrts) \
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
 \
            CAT3(_big_mul, N, M)(r, t, _crt_data_co_prime(Rcrts + np - 1, l, n), Xs[l*BLK_SZ + j]); \
            for (l++; l < np; l++) \
                CAT3(_big_addmul, N, M)(r, t, _crt_data_co_prime(Rcrts + np - 1, l, n), Xs[l*BLK_SZ + j]); \
 \
            CAT(_reduce_big_sum, N)(r, t, crt_data_prod_primes(Rcrts + np - 1)); \
 \
            if (mpn_cmp(r, Mhalf, N) > 0) \
            { \
                multi_rsub_ ## N(r, crt_data_prod_primes(Rcrts + np - 1)); \
                fmpz_neg_ui_array(&z[i+j-zl], r, N); \
            } \
            else \
                fmpz_set_ui_array(&z[i+j-zl], r, N); \
        } \
    } \
}

DEFINE_IT(2, 2, 1)  /* 100 bits */
DEFINE_IT(3, 3, 2)  /* 150 bits */
DEFINE_IT(4, 4, 3)  /* 200 bits */
DEFINE_IT(5, 4, 4)  /* 250 bits */
DEFINE_IT(6, 5, 4)  /* 300 bits */
DEFINE_IT(7, 6, 5)  /* 350 bits */
DEFINE_IT(8, 7, 6)  /* 400 bits */
#undef DEFINE_IT

/* 50 bits */
static void _crt_1(
    fmpz * z, ulong zl, ulong zi_start, ulong zi_stop,
    sd_fft_ctx_struct* Rffts, double* d, ulong dstride,
    crt_data_struct* Rcrts)
{
    ulong i, j, jstart, jstop;
    ulong Xs[BLK_SZ*1];
    ulong p = Rffts[0].mod.n;

    for (i = n_round_down(zi_start, BLK_SZ); i < zi_stop; i += BLK_SZ)
    {
        _convert_block(Xs, Rffts, d, dstride, 1, i/BLK_SZ);

        jstart = (i < zi_start) ? zi_start - i : 0;
        jstop = FLINT_MIN(BLK_SZ, zi_stop - i);

        for (j = jstart; j < jstop; j += 1)
        {
            if (COEFF_IS_MPZ(z[i+j-zl]))
                _fmpz_clear_mpz(z[i+j-zl]);

            z[i+j-zl] = (Xs[j] <= p / 2) ? Xs[j] : Xs[j] - p;
        }
    }
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
    const fmpz * a;
    ulong an;
    slong abits;
    const fmpz * b;
    ulong bn;
    slong bbits;
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
    _mod(X->bbuf, X->btrunc, X->b, X->bn, X->bbits, X->ffts + X->ioff);
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
                _mod(bbuf, X->btrunc, X->b, X->bn, X->bbits, X->ffts + ioff);
                sd_fft_lctx_fft_trunc(Q, bbuf, X->depth, X->btrunc, X->ztrunc);
            }
        }

        _mod(abuf, X->atrunc, X->a, X->an, X->abits, X->ffts + ioff);
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
    fmpz * z;
    ulong zl;
    ulong start_zi;
    ulong stop_zi;
    double* buf;
    ulong offset;
    ulong stride;
    sd_fft_ctx_struct* ffts;
    crt_data_struct* crts;
    nmod_t mod;
    void (*f)(
        fmpz * z, ulong zl, ulong zi_start, ulong zi_stop,
        sd_fft_ctx_struct* Rffts, double* d, ulong dstride,
        crt_data_struct* Rcrts);
} s2worker_struct;

static void s2worker_func(void* varg)
{
    s2worker_struct* X = (s2worker_struct*) varg;

    X->f(X->z, X->zl, X->start_zi, X->stop_zi, X->ffts + X->offset, X->buf,
         X->stride, X->crts + X->offset);
}

int _fmpz_poly_mul_mid_mpn_ctx(
    fmpz * z, ulong zl, ulong zh,
    const fmpz * a, ulong an,
    const fmpz * b, ulong bn,
    mpn_ctx_t R)
{
    ulong modbits;
    ulong offset = 0;
    ulong zn = an + bn - 1;
    ulong atrunc, btrunc, ztrunc;
    ulong i, np, depth, stride;
    double* buf;
    int squaring;
    slong bits1, bits2;
    int sign;

    FLINT_ASSERT(an > 0);
    FLINT_ASSERT(bn > 0);

    if (zl >= zh)
        return 1;

    if (zh > zn)
    {
        if (zl >= zn)
        {
            _fmpz_vec_zero(z, zh - zl);
            return 1;
        }

        _fmpz_vec_zero(z + zn - zl, zh - zn);
        zh = zn;
    }

    squaring = (a == b) && (an == bn);

    bits1 = _fmpz_vec_max_bits(a, an);

    if (squaring)
        bits2 = bits1;
    else
        bits2 = _fmpz_vec_max_bits(b, bn);
    sign = (bits1 < 0) || (bits2 < 0);
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
            return 0;

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
        X->abits = bits1;
        X->b = b;
        X->bn = bn;
        X->bbits = bits2;
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

    return 1;
}
