/*
    Copyright (C) 2022 Daniel Schultz

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

static void _mod_red(
    double* abuf, ulong atrunc,
    const ulong* a, ulong an,
    const sd_fft_ctx_struct* fft,
    nmod_t mod)
{
    double* aI;
    ulong i, j;

    FLINT_ASSERT(atrunc < an);
    FLINT_ASSERT(atrunc%BLK_SZ == 0);

#if 1

    ulong tt = an%atrunc;

#define UNROLL 8

    for (i = 0; i < atrunc; i += BLK_SZ)
    {
        aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);

        vec8n nn = vec8n_set_n(mod.n);
        vec8d n = vec8d_set_d(fft->p);
        vec8d ninv = vec8d_set_d(fft->pinv);

        for (j = 0; j < BLK_SZ; j += UNROLL)
        {
            if (i+j+UNROLL <= tt || i+j >= tt)
            {
                ulong k = i+j;
                FLINT_ASSERT(k+UNROLL-1 < an);
                vec8n t = vec8n_load_unaligned(a + k);

                if (mod.norm == 0)
                    for (k += atrunc; k < an; k += atrunc)
                        t = vec8n_addmod(t, vec8n_load_unaligned(a + k), nn);
                else
                    for (k += atrunc; k < an; k += atrunc)
                        t = vec8n_addmod_limited(t, vec8n_load_unaligned(a + k), nn);


                vec8d tlo = vec8n_convert_limited_vec8d(vec8n_bit_and(t, vec8n_set_n(n_pow2(32)-1)));
                vec8d thi = vec8n_convert_limited_vec8d(vec8n_bit_shift_right_32(t));
                vec8d_store_aligned(aI + j, vec8d_add(tlo, vec8d_mulmod(thi, vec8d_set_d(n_pow2(32)), n, ninv)));
            }
            else
            {
                for (ulong l = 0; l < UNROLL; l++)
                {
                    ulong k = i+j+l;
                    ulong c = a[k];
                    for (k += atrunc; k < an; k += atrunc)
                        c = nmod_add(c, a[k], mod);

                    aI[j+l] = (slong)(nmod_set_ui(c, fft->mod));
                }
            }
        }
    }

#else
    // wrong way!!
    if (modn <= fft->mod.n)
    {
        /* first pass fill in */

        for (i = 0; i < atrunc; i += BLK_SZ)
        {
            aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
            for (j = 0; j < BLK_SZ; j += 8)
            {
                vec8n t = vec8n_load_unaligned(a + i + j);
                vec8d_store_aligned(aI + j, vec8n_convert_limited_vec8d(t));
            }
        }

        vec8d n = vec8d_set_d(fft->p);

        /* second pass add to existing */

        for (i = atrunc; i + BLK_SZ <= an; i += BLK_SZ)
        {
            aI = sd_fft_ctx_blk_index(abuf, (i%atrunc)/BLK_SZ);
            for (j = 0; j < BLK_SZ; j += 8)
            {
                vec8n t = vec8n_load_unaligned(a + i + j);
                vec8d s = vec8n_convert_limited_vec8d(t);
                s = vec8d_add(s, vec8d_load_aligned(aI + j));
                s = vec8d_reduce_2n_to_n(s, n);
                vec8d_store_aligned(aI + j, s);
            }
        }

        aI = sd_fft_ctx_blk_index(abuf, (i%atrunc)/BLK_SZ);
        for (j = 0; j < an - i; j++)
            aI[j] = vec1d_reduce_2n_to_n(aI[j] + (slong)a[i + j], fft->p);
    }
    else
    {
        vec8d n = vec8d_set_d(fft->p);
        vec8d ninv = vec8d_set_d(fft->pinv);
        for (i = 0; i < atrunc; i += BLK_SZ)
        {
            aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
            for (j = 0; j < BLK_SZ; j += 8)
            {
                vec8n t = vec8n_load_unaligned(a + i + j);
                vec8d tlo = vec8n_convert_limited_vec8d(vec8n_bit_and(t, vec8n_set_n(n_pow2(32)-1)));
                vec8d thi = vec8n_convert_limited_vec8d(vec8n_bit_shift_right_32(t));
                vec8d_store_aligned(aI + j, vec8d_add(tlo, vec8d_mulmod(thi, vec8d_set_d(n_pow2(32)), n, ninv)));
            }
        }

        for (i = atrunc; i + BLK_SZ <= an; i += BLK_SZ)
        {
            aI = sd_fft_ctx_blk_index(abuf, (i%atrunc)/BLK_SZ);
            for (j = 0; j < BLK_SZ; j += 8)
            {
                vec8n t = vec8n_load_unaligned(a + i + j);
                vec8d tlo = vec8n_convert_limited_vec8d(vec8n_bit_and(t, vec8n_set_n(n_pow2(32)-1)));
                vec8d thi = vec8n_convert_limited_vec8d(vec8n_bit_shift_right_32(t));
                vec8d s = vec8d_add(tlo, vec8d_mulmod(thi, vec8d_set_d(n_pow2(32)), n, ninv));
                s = vec8d_add(s, vec8d_load_aligned(aI + j));
                s = vec8d_reduce_2n_to_n(s, n);
                vec8d_store_aligned(aI + j, s);
            }
        }

        aI = sd_fft_ctx_blk_index(abuf, (i%atrunc)/BLK_SZ);
        for (j = 0; j < an - i; j++)
            aI[j] = vec1d_reduce_2n_to_n(aI[j] + (slong)nmod_set_ui(a[i+j], fft->mod), fft->p);
    }
#endif
}

static void _mod(
    double* abuf, ulong atrunc,
    const ulong* a, ulong an,
    const sd_fft_ctx_struct* fft,
    nmod_t mod)
{
    double* aI;
    ulong i, j;

    if (atrunc < an)
    {
        _mod_red(abuf, atrunc, a, an, fft, mod);
        return;
    }

    if (mod.n <= fft->mod.n)
    {
        for (i = 0; i + BLK_SZ <= an; i += BLK_SZ)
        {
            aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
            for (j = 0; j < BLK_SZ; j += 8)
            {
                vec8n t = vec8n_load_unaligned(a + i + j);
FLINT_ASSERT(i+j < atrunc);
                vec8d_store_aligned(aI + j, vec8n_convert_limited_vec8d(t));
            }
        }

        aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
        for (j = 0; j < an - i; j++)
            aI[j] = (slong)a[i + j];
    }
    else
    {
        vec8d n = vec8d_set_d(fft->p);
        vec8d ninv = vec8d_set_d(fft->pinv);
        for (i = 0; i + BLK_SZ <= an; i += BLK_SZ)
        {
            aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
#if 1
            for (j = 0; j < BLK_SZ; j += 8)
            {
                vec8n t = vec8n_load_unaligned(a + i + j);
                vec8d tlo = vec8n_convert_limited_vec8d(vec8n_bit_and(t, vec8n_set_n(n_pow2(32)-1)));
                vec8d thi = vec8n_convert_limited_vec8d(vec8n_bit_shift_right_32(t));
                vec8d_store_aligned(aI + j, vec8d_add(tlo, vec8d_mulmod(thi, vec8d_set_d(n_pow2(32)), n, ninv)));
            }
#else
            for (j = 0; j < BLK_SZ; j += 1)
            {
                aI[j] = (slong)nmod_set_ui(a[i+j], fft->mod);
            }
#endif
        }

        aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
        for (j = 0; j < an - i; j++)
            aI[j] = (slong)nmod_set_ui(a[i+j], fft->mod);
    }

    for (i = an; i < atrunc; i++)
        sd_fft_ctx_set_index(abuf, i, 0);
}


#define DEFINE_IT(NP, N, M) \
static void CAT(_crt, NP)( \
    ulong* z, ulong zl, ulong zi_start, ulong zi_stop, \
    sd_fft_ctx_struct* Rffts, double* d, ulong dstride, \
    crt_data_struct* Rcrts, \
    nmod_t mod) \
{ \
    ulong np = NP; \
    ulong n = N; \
    ulong m = M; \
 \
    FLINT_ASSERT(n == Rcrts[np-1].coeff_len); \
    FLINT_ASSERT(1 <= N && N <= 3); \
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
            if (N == 1) \
            { \
                NMOD_RED(z[i+j-zl], r[0], mod); \
            } \
            else if (N == 2) \
            { \
                NMOD2_RED2(z[i+j-zl], r[1], r[0], mod); \
            } \
            else \
            { \
                FLINT_ASSERT(N < 4 || r[3] == 0); \
                NMOD_RED3(z[i+j-zl], r[2], r[1], r[0], mod); \
            } \
        } \
    } \
}

DEFINE_IT(2, 2, 1)
DEFINE_IT(3, 3, 2)
DEFINE_IT(4, 4, 3)
#undef DEFINE_IT

static void _crt_1(
    ulong* z, ulong zl, ulong zi_start, ulong zi_stop,
    sd_fft_ctx_struct* Rffts, double* d, ulong dstride,
    crt_data_struct* Rcrts,
    nmod_t mod)
{
    ulong i, j, jstart, jstop;
    ulong Xs[BLK_SZ*1];

    if (mod.n == Rffts[0].mod.n)
    {
        for (i = n_round_down(zi_start, BLK_SZ); i < zi_stop; i += BLK_SZ)
        {
            _convert_block(Xs, Rffts, d, dstride, 1, i/BLK_SZ);

            jstart = (i < zi_start) ? zi_start - i : 0; \
            jstop = FLINT_MIN(BLK_SZ, zi_stop - i);
            for (j = jstart; j < jstop; j += 1)
            {
                z[i+j-zl] = Xs[j];
            }
        }
    }
    else
    {
        for (i = n_round_down(zi_start, BLK_SZ); i < zi_stop; i += BLK_SZ)
        {
            _convert_block(Xs, Rffts, d, dstride, 1, i/BLK_SZ);

            jstart = (i < zi_start) ? zi_start - i : 0; \
            jstop = FLINT_MIN(BLK_SZ, zi_stop - i);

            for (j = jstart; j < jstop; j += 1)
            {
                NMOD_RED(z[i+j-zl], Xs[j], mod);
            }
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
    const ulong* a;
    ulong an;
    const ulong* b;
    ulong bn;
    sd_fft_ctx_struct* ffts;
    crt_data_struct* crts;
    nmod_t mod;
    ulong ioff;
    int squaring;
} s1worker_struct;


static void extra_func(void* varg)
{
    s1worker_struct* X = (s1worker_struct*) varg;
    sd_fft_lctx_t Q;

    sd_fft_lctx_init(Q, X->ffts + X->ioff, X->depth);
    _mod(X->bbuf, X->btrunc, X->b, X->bn, X->ffts + X->ioff, X->mod);
    sd_fft_lctx_fft_trunc(Q, X->bbuf, X->depth, X->btrunc, X->ztrunc);
    sd_fft_lctx_clear(Q, X->ffts + X->ioff);
}

void s1worker_func(void* varg)
{
    s1worker_struct* X = (s1worker_struct*) varg;
    sd_fft_lctx_t Q;
    ulong i, m;
    thread_pool_handle* handles = NULL;
    slong nworkers = 0;

    if (X->bn > 5000 && !X->squaring)
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
                _mod(bbuf, X->btrunc, X->b, X->bn, X->ffts + ioff, X->mod);
                sd_fft_lctx_fft_trunc(Q, bbuf, X->depth, X->btrunc, X->ztrunc);
            }
        }

        _mod(abuf, X->atrunc, X->a, X->an, X->ffts + ioff, X->mod);
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
    ulong* z;
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
        ulong* z, ulong zl, ulong zi_start, ulong zi_stop,
        sd_fft_ctx_struct* Rffts, double* d, ulong dstride,
        crt_data_struct* Rcrts,
        nmod_t mod);
} s2worker_struct;

void s2worker_func(void* varg)
{
    s2worker_struct* X = (s2worker_struct*) varg;

    X->f(X->z, X->zl, X->start_zi, X->stop_zi, X->ffts + X->offset, X->buf,
         X->stride, X->crts + X->offset, X->mod);
}

void _nmod_poly_mul_mid_mpn_ctx(
    ulong* z, ulong zl, ulong zh,
    const ulong* a, ulong an,
    const ulong* b, ulong bn,
    nmod_t mod,
    mpn_ctx_t R)
{
    ulong modbits = FLINT_BITS - mod.norm;
    ulong offset = 0;
    ulong zn = an + bn - 1;
    ulong atrunc, btrunc, ztrunc;
    ulong i, np, depth, stride;
    double* buf;
    int squaring;

    FLINT_ASSERT(an > 0);
    FLINT_ASSERT(bn > 0);

    if (zl >= zh)
        return;

    if (zh > zn)
    {
        if (zl >= zn)
        {
            flint_mpn_zero(z, zh - zl);
            return;
        }

        flint_mpn_zero(z + zn - zl, zh - zn);
        zh = zn;
    }

    squaring = (a == b) && (an == bn);

    FLINT_ASSERT(zl < zh);
    FLINT_ASSERT(zh <= zn);

    /* first see if mod.n is on of R->ffts[i].mod.n */
    if (modbits == 50)
    {
        for (i = 0; i < MPN_CTX_NCRTS; i++)
        {
            if (mod.n == R->ffts[i].mod.n)
            {
                offset = i;
                np = 1;
                goto got_np_and_offset;
            }
        }
    }

    /* need prod_of_primes >= blen * 4^modbits */
    for (np = 1; np < 4; np++)
    {
        if (flint_mpn_cmp_ui_2exp(crt_data_prod_primes(R->crts + np - 1),
              R->crts[np - 1].coeff_len, bn, 2*modbits) >= 0)
        {
            break;
        }
    }

    FLINT_ASSERT(0 <= flint_mpn_cmp_ui_2exp(
                                  crt_data_prod_primes(R->crts + np - 1),
                                  R->crts[np - 1].coeff_len, bn, 2*modbits));


got_np_and_offset:

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

    if (bn < 1500)
        want_threads = 1;
    else
        want_threads = np;

    thread_pool_handle* handles;
    slong nworkers = flint_request_threads(&handles, want_threads);
    ulong nthreads = nworkers + 1;

    buf = (double*) mpn_ctx_fit_buffer(R, (np+nthreads)*stride*sizeof(double));

    s1worker_struct s1args[4];
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
        X->ffts = R->ffts;
        X->crts = R->crts;
        X->mod = mod;
        X->squaring = squaring;
    }

    for (i = nworkers; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0, s1worker_func, s1args + i);
    s1worker_func(s1args + 0);
    for (i = nworkers; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    if (np*zn > 10000)
    {
        flint_give_back_threads(handles, nworkers);
        nworkers = flint_request_threads(&handles, 8);
        nthreads = nworkers + 1;
    }

    s2worker_struct s2args[8];
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
        X->mod = mod;
        X->f = np == 1 ? _crt_1 : np == 2 ? _crt_2 : np == 3 ? _crt_3 : _crt_4;
    }

    for (i = nworkers; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0, s2worker_func, s2args + i);
    s2worker_func(s2args + 0);
    for (i = nworkers; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    flint_give_back_threads(handles, nworkers);
}

void _nmod_poly_mul_mod_xpnm1_naive(
    ulong* z, ulong zn,
    const ulong* a, ulong an,
    const ulong* b, ulong bn,
    ulong lgN,
    nmod_t mod,
    mpn_ctx_t R)
{
    ulong N = n_pow2(lgN);
    FLINT_ASSERT(zn <= N);

    ulong* t = FLINT_ARRAY_ALLOC(an + bn - 1, ulong);

    if (an >= bn)
        _nmod_poly_mul(t, a, an, b, bn, mod);
    else
        _nmod_poly_mul(t, b, bn, a, an, mod);

    for (ulong i = 0; i < zn; i++)
    {
        ulong c = 0;
        for (ulong j = i; j < an + bn - 1; j += N)
            c = nmod_add(c, t[j], mod);
        z[i] = c;
    }

    flint_free(t);
}


void _nmod_poly_mul_mod_xpnm1(
    ulong* z, ulong ztrunc,
    const ulong* a, ulong an,
    const ulong* b, ulong bn,
    ulong depth,
    nmod_t mod,
    mpn_ctx_t R)
{
    ulong N = n_pow2(depth);
    ulong modbits = FLINT_BITS - mod.norm;
    ulong offset = 0;
    ulong i, np, stride;
    double* buf;

    FLINT_ASSERT(an > 0);
    FLINT_ASSERT(bn > 0);
    FLINT_ASSERT(ztrunc <= N);

    /* first see if mod.n is on of R->ffts[i].mod.n */

    if (modbits == 50)
    {
        for (i = 0; i < MPN_CTX_NCRTS; i++)
        {
            if (mod.n == R->ffts[i].mod.n)
            {
                offset = i;
                np = 1;
                goto got_np_and_offset;
            }
        }
    }

    /* need prod_of_primes >= N * 4^modbits */
    for (np = 1; np < 4; np++)
    {
        if (flint_mpn_cmp_ui_2exp(crt_data_prod_primes(R->crts + np - 1),
              R->crts[np - 1].coeff_len, N, 2*modbits) >= 0)
        {
            break;
        }
    }

    FLINT_ASSERT(0 <= flint_mpn_cmp_ui_2exp(
                                  crt_data_prod_primes(R->crts + np - 1),
                                  R->crts[np - 1].coeff_len, N, 2*modbits));


got_np_and_offset:

    stride = n_round_up(sd_fft_ctx_data_size(depth), 128);

    ulong want_threads;

    if (bn < 1500)
        want_threads = 1;
    else
        want_threads = np;

    thread_pool_handle* handles;
    slong nworkers = flint_request_threads(&handles, want_threads);
    ulong nthreads = nworkers + 1;

    buf = (double*) mpn_ctx_fit_buffer(R, (np+nthreads)*stride*sizeof(double));

    s1worker_struct s1args[4];
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
        X->atrunc = N;
        X->btrunc = N;
        X->ztrunc = N;
        X->a = a;
        X->an = an;
        X->b = b;
        X->bn = bn;
        X->ffts = R->ffts;
        X->crts = R->crts;
        X->mod = mod;
        X->squaring = 0;
    }

    for (i = nworkers; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0, s1worker_func, s1args + i);
    s1worker_func(s1args + 0);
    for (i = nworkers; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    s2worker_struct s2args[8];
    ulong zl = 0;
    ulong zh = ztrunc;
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
        X->mod = mod;
        X->f = np == 1 ? _crt_1 : np == 2 ? _crt_2 : np == 3 ? _crt_3 : _crt_4;
    }

    for (i = nworkers; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0, s2worker_func, s2args + i);
    s2worker_func(s2args + 0);
    for (i = nworkers; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    flint_give_back_threads(handles, nworkers);
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
    ulong ztrunc;
    const ulong* a;
    ulong an;
    sd_fft_ctx_struct* ffts;
    crt_data_struct* crts;
    nmod_t mod;
} s1pworker_struct;

void s1pworker_func(void* varg)
{
    s1pworker_struct* X = (s1pworker_struct*) varg;
    sd_fft_lctx_t Q;
    ulong i, m;

    for (i = X->start_pi; i < X->stop_pi; i++)
    {
        ulong ioff = i + X->offset;
        double* abuf = X->abuf + X->stride*i;

        sd_fft_lctx_init(Q, X->ffts + ioff, X->depth);

        _mod(abuf, X->atrunc, X->a, X->an, X->ffts + ioff, X->mod);
        sd_fft_lctx_fft_trunc(Q, abuf, X->depth, X->atrunc, X->ztrunc);

        ulong cop = X->np == 1 ? 1 : *crt_data_co_prime_red(X->crts + X->np - 1, ioff);
        NMOD_RED2(m, cop >> (FLINT_BITS - X->depth), cop << X->depth, X->ffts[ioff].mod);
        m = nmod_inv(m, X->ffts[ioff].mod);
        sd_fft_lctx_point_mul(Q, abuf, X->bbuf + X->stride*i, m, X->depth);

        sd_fft_lctx_ifft_trunc(Q, abuf, X->depth, X->ztrunc);

        sd_fft_lctx_clear(Q, X->ffts + ioff);
    }
}

void _mul_precomp_init(
    mul_precomp_struct* M,
    const ulong * b, ulong bn, ulong btrunc,
    ulong depth,
    nmod_t mod,
    mpn_ctx_t R)
{
    ulong N = n_pow2(depth);
    ulong modbits = FLINT_BITS - mod.norm;
    ulong offset = 0;
    ulong i, np, stride;
    sd_fft_lctx_t Q;

    btrunc = n_round_up(btrunc, BLK_SZ);

    FLINT_ASSERT(bn > 0);

    /* first see if mod.n is on of R->ffts[i].mod.n */

    if (modbits == 50)
    {
        for (i = 0; i < MPN_CTX_NCRTS; i++)
        {
            if (mod.n == R->ffts[i].mod.n)
            {
                offset = i;
                np = 1;
                goto got_np_and_offset;
            }
        }
    }

    /* need prod_of_primes >= N * 4^modbits */
    for (np = 1; np < 4; np++)
    {
        if (flint_mpn_cmp_ui_2exp(crt_data_prod_primes(R->crts + np - 1),
              R->crts[np - 1].coeff_len, N, 2*modbits) >= 0)
        {
            break;
        }
    }

    FLINT_ASSERT(0 <= flint_mpn_cmp_ui_2exp(
                                  crt_data_prod_primes(R->crts + np - 1),
                                  R->crts[np - 1].coeff_len, N, 2*modbits));

got_np_and_offset:

    stride = n_round_up(sd_fft_ctx_data_size(depth), 128);

    M->depth = depth;
    M->N = N;
    M->offset = offset;
    M->np = np;
    M->stride = stride;
    M->bn = bn;
    M->btrunc = btrunc;
    M->bbuf = flint_aligned_alloc(4096, n_round_up(np*stride*sizeof(double), 4096));

    for (i = 0; i < np; i++)
    {
        ulong ioff = i + offset;
        double* bbuf = M->bbuf + stride*i;

        sd_fft_lctx_init(Q, R->ffts + ioff, depth);

        _mod(bbuf, N, b, bn, R->ffts + ioff, mod);
        sd_fft_lctx_fft_trunc(Q, bbuf, depth, N, N);

        sd_fft_lctx_clear(Q, R->ffts + ioff);
    }
}


int _nmod_poly_mul_mid_precomp(
    ulong* z, ulong zl, ulong zh,
    const ulong* a, ulong an,
    mul_precomp_struct* M,
    nmod_t mod,
    mpn_ctx_t R)
{
    ulong bn = M->bn;
    ulong zn = an + bn - 1;
    ulong atrunc, ztrunc;
    ulong depth = M->depth;
    ulong N = n_pow2(depth);
    ulong i, np = M->np;
    double* buf;

    FLINT_ASSERT(an > 0);
    FLINT_ASSERT(bn > 0);
    FLINT_ASSERT(M->btrunc <= N);

    if (zl >= zh)
        return 1;

    if (zh > zn)
    {
        if (zl >= zn)
        {
            flint_mpn_zero(z, zh - zl);
            return 1;
        }

        flint_mpn_zero(z + zn - zl, zh - zn);
        zh = zn;
    }

    FLINT_ASSERT(zl < zh);
    FLINT_ASSERT(zh <= zn);

    atrunc = n_round_up(an, BLK_SZ);
    ztrunc = n_round_up(zn, BLK_SZ);

    if (atrunc <= N && zh <= N && N <= zn && zn <= zl + N)
    {
        ztrunc = N;
    }
    else if (ztrunc <= N)
    {
    }
    else
    {
        return 0;
    }

    ulong want_threads;

    /* todo: even with np == 1, we may want threads for the conversion to fft */
    if (bn < 1500)
        want_threads = 1;
    else
        want_threads = np;

    thread_pool_handle* handles;
    slong nworkers = flint_request_threads(&handles, want_threads);
    ulong nthreads = nworkers + 1;

    buf = (double*) mpn_ctx_fit_buffer(R, np*M->stride*sizeof(double));

    s1pworker_struct s1pargs[4];
    for (i = 0; i < nthreads; i++)
    {
        s1pworker_struct* X = s1pargs + i;
        X->np = np;
        X->start_pi = (i+0)*np/nthreads;
        X->stop_pi  = (i+1)*np/nthreads;
        X->offset = M->offset;
        X->abuf = buf;
        X->bbuf = M->bbuf;
        X->depth = depth;
        X->stride = M->stride;
        X->atrunc = atrunc;
        X->ztrunc = ztrunc;
        X->a = a;
        X->an = an;
        X->ffts = R->ffts;
        X->crts = R->crts;
        X->mod = mod;
    }

    for (i = nworkers; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0, s1pworker_func, s1pargs + i);
    s1pworker_func(s1pargs + 0);
    for (i = nworkers; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    s2worker_struct s2args[8];
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
        X->offset = M->offset;
        X->stride = M->stride;
        X->ffts = R->ffts;
        X->crts = R->crts;
        X->mod = mod;
        X->f = np == 1 ? _crt_1 : np == 2 ? _crt_2 : np == 3 ? _crt_3 : _crt_4;
    }

    for (i = nworkers; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0, s2worker_func, s2args + i);
    s2worker_func(s2args + 0);
    for (i = nworkers; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    flint_give_back_threads(handles, nworkers);
    return 1;
}

void _nmod_poly_mul_mod_xpnm1_precomp(
    ulong* z, ulong ztrunc,
    const ulong* a, ulong an,
    mul_precomp_struct* M,
    nmod_t mod,
    mpn_ctx_t R)
{
    ulong depth = M->depth;
    ulong i, np = M->np;
    double* buf;

    FLINT_ASSERT(an > 0);

    ulong want_threads;

    if (an < 1500)
        want_threads = 1;
    else
        want_threads = np;

    thread_pool_handle* handles;
    slong nworkers = flint_request_threads(&handles, want_threads);
    ulong nthreads = nworkers + 1;

    buf = (double*) mpn_ctx_fit_buffer(R, np*M->stride*sizeof(double));

    s1pworker_struct s1pargs[4];
    for (i = 0; i < nthreads; i++)
    {
        s1pworker_struct* X = s1pargs + i;
        X->np = np;
        X->start_pi = (i+0)*np/nthreads;
        X->stop_pi  = (i+1)*np/nthreads;
        X->offset = M->offset;
        X->abuf = buf;
        X->bbuf = M->bbuf;
        X->depth = depth;
        X->stride = M->stride;
        X->atrunc = M->btrunc;
        X->ztrunc = M->btrunc;
        X->a = a;
        X->an = an;
        X->ffts = R->ffts;
        X->crts = R->crts;
        X->mod = mod;
    }

    for (i = nworkers; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0, s1pworker_func, s1pargs + i);
    s1pworker_func(s1pargs + 0);
    for (i = nworkers; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    s2worker_struct s2args[8];
    ulong zl = 0;
    ulong zh = ztrunc;
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
        X->offset = M->offset;
        X->stride = M->stride;
        X->ffts = R->ffts;
        X->crts = R->crts;
        X->mod = mod;
        X->f = np == 1 ? _crt_1 : np == 2 ? _crt_2 : np == 3 ? _crt_3 : _crt_4;
    }

    for (i = nworkers; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0, s2worker_func, s2args + i);
    s2worker_func(s2args + 0);
    for (i = nworkers; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    flint_give_back_threads(handles, nworkers);
}


/* z -= a mod x^N-1, write coeffs [0,ztrunc) */
void _nmod_poly_sub_mod_xpNm1(
    ulong* z, ulong ztrunc,
    const ulong* a, ulong an,
    ulong N, nmod_t mod)
{
    FLINT_ASSERT(ztrunc <= an);
    FLINT_ASSERT(ztrunc <= N);

    for (ulong i = 0; i < ztrunc; i++)
    {
        ulong k = i;
        ulong c = nmod_sub(a[k], z[i], mod);
        for (k += N; k < an; k += N)
            c = nmod_add(c, a[k], mod);
        z[i] = c;
    }
}

/*
    for divrem(a, b)

    an = length(a)
    bn = length(b)
    qn = length(q) = an - bn + 1

    choose a precision Bn of B(x) = B[0] + ... + B[Bn-1]*x^(Bn-1) with

        rev(B) = rev(b)^-1 mod x^Bn = B[Bn-1] + B[Bn-2]*x + ... + B[0]*x^(Bn-1)

    then
        (a[an-1] + a[an-2]*x + ... + a[an-qn]*x^(qn-1))
       *
        (B[Bn-1] + B[Bn-2]*x + ... + B[0]*x^(Bn-1))
       =
        q[qn-1] + q[qn-2]*x + ... + q[0]*x^(qn-1)  mod x^qn

    therefore need Bn >= qn, or, the same thing, Bn >= an - bn + 1

    in terms of non-reversed polys,

        _mul_mid(q, an+Bn-1-qn, an+Bn-1, a, an, B, Bn)

    or, the same thing,

        _mul_mid(q, Bn-1, Bn-1+qn, a+an-qn, qn, B, Bn)

    will calculate q. Then, find r via

        r = a - b*q mod x^N-1 where N >= bn - 1
*/
void _nmod_poly_divrem_mpn_ctx(
    ulong* q,
    ulong* r,
    const ulong* a, ulong an,
    const ulong* b, ulong bn,
    nmod_t mod,
    mpn_ctx_t R)
{
    ulong qn = an - bn + 1;

    FLINT_ASSERT(an >= bn);
    FLINT_ASSERT(bn > 1);
    FLINT_ASSERT(qn > 0);

    /* choose precision */
    ulong Bn = qn;

    ulong lgN = n_max(LG_BLK_SZ, n_clog2(bn-1));
    ulong N = n_pow2(lgN);

    ulong* B = FLINT_ARRAY_ALLOC(Bn, ulong);
    ulong* t = FLINT_ARRAY_ALLOC(FLINT_MAX(N, bn), ulong);

    _nmod_poly_reverse(t, b, bn, bn);
    _nmod_poly_inv_series(B, t, bn, Bn, mod);
    _nmod_poly_reverse(B, B, Bn, Bn);

    _nmod_poly_mul_mid_mpn_ctx(q, Bn-1, Bn-1+qn, a+an-qn, qn, B, Bn, mod, R);
    _nmod_poly_mul_mod_xpnm1(r, bn-1, q, qn, b, bn, lgN, mod, R);
    _nmod_poly_sub_mod_xpNm1(r, bn-1, a, an, N, mod);

    flint_free(B);
    flint_free(t);
}

void _nmod_poly_divrem_precomp_init(
    nmod_poly_divrem_precomp_struct* M,
    const ulong* b, ulong bn,
    ulong Bn,
    nmod_t mod,
    mpn_ctx_t R)
{
    ulong* B = FLINT_ARRAY_ALLOC(Bn, ulong);
    ulong* t = FLINT_ARRAY_ALLOC(bn, ulong);

    _nmod_poly_reverse(t, b, bn, bn);
    _nmod_poly_inv_series(B, t, bn, Bn, mod);
    _nmod_poly_reverse(B, B, Bn, Bn);


    _mul_precomp_init(M->quo_maker, B, Bn, Bn, n_max(LG_BLK_SZ, n_clog2(2*Bn-1)), mod, R);

    ulong lgN = n_max(LG_BLK_SZ, n_clog2(bn-1));
    ulong N = n_pow2(lgN);
    _mul_precomp_init(M->rem_maker, b, bn, N, lgN, mod, R);

    flint_free(B);
    flint_free(t);
}


int _nmod_poly_divrem_precomp(
    ulong* q,
    ulong* r,
    const ulong* a, ulong an,
    nmod_poly_divrem_precomp_struct* M,
    nmod_t mod,
    mpn_ctx_t R)
{
    ulong Bn = M->quo_maker->bn;
    ulong bn = M->rem_maker->bn;
    ulong qn = an-bn+1;

    FLINT_ASSERT(an >= bn);
    FLINT_ASSERT(bn > 1);
    FLINT_ASSERT(qn > 0);

    if (!_nmod_poly_mul_mid_precomp(q, Bn-1, Bn-1+qn, a+an-qn, qn, M->quo_maker, mod, R))
        return 0;
    _nmod_poly_mul_mod_xpnm1_precomp(r, bn-1, q, qn, M->rem_maker, mod, R);
    _nmod_poly_sub_mod_xpNm1(r, bn-1, a, an, M->rem_maker->N, mod);
    return 1;
}



/*
definition of _mul_mid(z, zl, zh, a, an, b, bn)

              h
[sum z[i]*x^i]  :=  sum  z[i]*x^(i-l)
  i           l    l<=i<h

i.e. the coeffs in [zl, zh)
*/
void _nmod_poly_mul_mid_classical(
    ulong* z, slong zl, slong zh,
    const ulong* a, slong an,
    const ulong* b, slong bn,
    nmod_t mod)
{
    for (slong i = zl; i < zh; i++)
    {
        slong jstart = z_max(0, i - (bn - 1));
        slong jstop = z_min(i + 1, an);
        ulong zi = 0;
        for (slong j = jstart; j < jstop; j++)
            zi = nmod_addmul(zi, a[j], b[i - j], mod);
        z[i - zl] = zi;
    }
}

/*
**** Karasuba ****

with deg(Ai), deg(Bi) < k, consider

P := (A0 + A1*x^k)*(B0 + B1*x^k)

define

P2 = A1*B1
P1 = (A0 + A1)*(B0 + B1)
P0 = A0*B0

then

# P = P0 + (P1 - P2 - P0)*x^k + P2*x^2k
*/




/*
**** Karasuba for middle product ****

with deg(Ai), deg(Bi) < k, consider
P := (A0 + A1*x^k + A2*x^2k + A3*x^3k)*(B0 + B1*x^k)

                            h
we would like to compute [P]
                            l
define

P0 := ((A0 + A1) + (A1 + A2)*x^k)*(B1)
P1 := (A1 + A2*x^k)*(B0 - B1)
P2 := ((A1 + A2) + (A2 + A3)*x^k)*(B0)

so that

P0 = (A0*B1 + A1*B1) + (A1*B1 + A2*B1)*x^k
P1 = (A1*B0 - A1*B1) + (A2*B0 - A2*B1)*x^k
P2 = (A1*B0 + A2*B0) + (A2*B0 + A3*B0)*x^k

and

   h            2k             h-2k
[P]  = [P0 + P1]    + [P2 - P1]    * x^(3k-l)
   l            l-k            k

In order to calculate the rhs, we need

    2k         h-2k        max(2k,h-2k)
[P0]     , [P2]     ,  [P1]
    l-k        k           min(l-k,k)


requies k <= l < 3k < h <= 4k

*/


void _nmod_poly_mul_mid_unbalanced(
    ulong* z, slong zl, slong zh,
    const ulong* a, slong an,
    const ulong* b, slong bn,
    nmod_t mod)
{
    FLINT_ASSERT(zl < zh);
    FLINT_ASSERT(bn < an);
    flint_mpn_zero(z, zh - zl);

    ulong* t = FLINT_ARRAY_ALLOC(2*bn, ulong);

    slong i;
    for (i = 0; i*bn < an; i++)
    {
        slong zl_new, zh_new, an_new;

        // produces coefficient for powers x^[zl_new + i*k, zh_new + i*k)

        zl_new = z_max(zl - i*bn, 0);
        an_new = z_min(bn, an - i*bn);
        zh_new = z_min(zh - i*bn, an_new + bn - 1);

        _nmod_poly_mul_mid(t, zl_new, zh_new, a + i*bn, an_new, b, bn, mod);
        ulong* Z = z + zl_new + i*bn - zl;
        _nmod_vec_add(Z, Z, t, zh_new - zl_new, mod);
    }

    flint_free(t);
}


void _nmod_poly_mul_mid(
    ulong* z, slong zl, slong zh,
    const ulong* a, slong an,
    const ulong* b, slong bn,
    nmod_t mod)
{
    if (zl >= zh)
        return;

    if (an < bn)
    {
        FLINT_SWAP(const ulong *, a, b);
        FLINT_SWAP(ulong, an, bn);
    }

    if (zl > bn - 1)
    {
        if (an > zl - (bn - 1))
        {
            an -= zl - (bn - 1);
            a  += zl - (bn - 1);
            zh -= zl - (bn - 1);
            zl -= zl - (bn - 1);
            _nmod_poly_mul_mid(z, zl, zh, a, an, b, bn, mod);
        }
        else
        {
            flint_mpn_zero(z, zh - zl);
        }
        return;
    }

    if (zh < an)
    {
        an = zh;
        _nmod_poly_mul_mid(z, zl, zh, a, an, b, bn, mod);
        return;
    }

    FLINT_ASSERT(zl < bn && bn <= an && an <= zh);

    if (an >= 2*bn)
    {
        _nmod_poly_mul_mid_unbalanced(z, zl, zh, a, an, b, bn, mod);
        return;
    }

    if (zl < an + bn + 1)
    {
        if (zh > 0)
        {
            /*
                middle product or three pieces
                +----------------+
                |             |\ |
                |      1      | \|
                |             |3 |
                |-------------+--|
                |\     2         |
                +----------------+
            */
        }
        else if (0)
        {
            /*
                two pieces
                +----------------+
                |             |\ |
                |             | \|
                |      1      |2 |
                |             |  |
                |             |  |
                +-------------+--+
            */
        }
        else
        {
            /*
                two pieces
                +----------------+
                |           \    |
                |            \   |
                |------------ \  |
                |              \ |
                |               \|
                +----------------+
            */
        }
    }
    else
    {
        if (zh > 0)
        {
            /*
                two pieces
                +----------------+
                | |              |
                | |              |
                |1|      2       |
                | |              |
                |\|              |
                +----------------+
            */

        }
    }

    _nmod_poly_mul_mid_classical(z, zl, zh, a, an, b, bn, mod);
    return;
}


