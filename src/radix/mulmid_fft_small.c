/*
    Copyright (C) 2022 Daniel Schultz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"

/* allow changing this from the test code */
FLINT_DLL int radix_mulmid_force_threading = 0;

#if FLINT_HAVE_FFT_SMALL

#include "thread_pool.h"
#include "thread_support.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "crt_helpers.h"
#include "fft_small.h"

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
        abuf[i] = 0;
}

#define DEFINE_IT(NP, N, M) \
static void CAT(_crt, NP)( \
    ulong* z, ulong zl, ulong zi_start, ulong zi_stop, \
    sd_fft_ctx_struct* Rffts, double* d, ulong dstride, \
    crt_data_struct* Rcrts, ulong * carry_out, \
    nmod_t mod) \
{ \
    ulong np = NP; \
    ulong n = N; \
    ulong m = M; \
    ulong cy[3] = { 0, 0, 0 }; \
 \
    FLINT_ASSERT(n == Rcrts[np-1].coeff_len); \
    FLINT_ASSERT(1 <= N && N <= 4); \
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
            if (N == 2) \
            { \
                add_ssaaaa(r[1], r[0], r[1], r[0], cy[1], cy[0]); \
                if (mod.norm == 0) \
                    z[i+j-zl] = flint_mpn_divrem_2_1_preinv_norm(cy, r, mod.n, mod.ninv); \
                else \
                    z[i+j-zl] = flint_mpn_divrem_2_1_preinv_unnorm(cy, r, mod.n, mod.ninv, mod.norm); \
            } \
            else \
            { \
                FLINT_ASSERT(N < 4 || r[3] == 0); \
                add_sssaaaaaa(r[2], r[1], r[0], r[2], r[1], r[0], cy[2], cy[1], cy[0]); \
                if (mod.norm == 0) \
                    z[i+j-zl] = flint_mpn_divrem_3_1_preinv_norm(cy, r, mod.n, mod.ninv); \
                else \
                    z[i+j-zl] = flint_mpn_divrem_3_1_preinv_unnorm(cy, r, mod.n, mod.ninv, mod.norm); \
            } \
        } \
    } \
    carry_out[0] = cy[0]; \
    carry_out[1] = cy[1]; \
    carry_out[2] = cy[2]; \
}

DEFINE_IT(2, 2, 1)
DEFINE_IT(3, 3, 2)
DEFINE_IT(4, 4, 3)
#undef DEFINE_IT

static void _crt_1(
    ulong* z, ulong zl, ulong zi_start, ulong zi_stop,
    sd_fft_ctx_struct* Rffts, double* d, ulong dstride,
    crt_data_struct* FLINT_UNUSED(Rcrts),
    ulong * carry_out,
    nmod_t mod)
{
    ulong i, j, jstart, jstop;
    ulong Xs[BLK_SZ*1];
    ulong cy = 0;

    for (i = n_round_down(zi_start, BLK_SZ); i < zi_stop; i += BLK_SZ)
    {
        _convert_block(Xs, Rffts, d, dstride, 1, i/BLK_SZ);

        jstart = (i < zi_start) ? zi_start - i : 0; \
        jstop = FLINT_MIN(BLK_SZ, zi_stop - i);

        FLINT_ASSERT(mod.norm != 0);

        for (j = jstart; j < jstop; j += 1)
        {
            z[i+j-zl] = n_divrem_preinv_unnorm(&cy, Xs[j] + cy, mod.n, mod.ninv, mod.norm);
        }
    }

    carry_out[0] = cy;
    carry_out[1] = 0;
    carry_out[2] = 0;
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
    sd_fft_ctx_struct* Q = X->ffts + X->ioff;

    _mod(X->bbuf, X->btrunc, X->b, X->bn, Q, X->mod);
    sd_fft_trunc(Q, X->bbuf, X->depth, X->btrunc, X->ztrunc);
}

static void s1worker_func(void* varg)
{
    s1worker_struct* X = (s1worker_struct*) varg;
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
        sd_fft_ctx_struct* Q = X->ffts + ioff;

        if (!X->squaring)
        {
            if (nworkers > 0)
            {
                X->ioff = ioff;
                thread_pool_wake(global_thread_pool, handles[0], 0, extra_func, X);
            }
            else
            {
                _mod(bbuf, X->btrunc, X->b, X->bn, Q, X->mod);
                sd_fft_trunc(Q, bbuf, X->depth, X->btrunc, X->ztrunc);
            }
        }

        _mod(abuf, X->atrunc, X->a, X->an, Q, X->mod);
        sd_fft_trunc(Q, abuf, X->depth, X->atrunc, X->ztrunc);

        if (!X->squaring)
        {
            if (nworkers > 0)
                thread_pool_wait(global_thread_pool, handles[0]);
        }

        ulong cop = X->np == 1 ? 1 : *crt_data_co_prime_red(X->crts + X->np - 1, ioff);
        NMOD_RED2(m, cop >> (FLINT_BITS - X->depth), cop << X->depth, Q->mod);
        m = nmod_inv(m, Q->mod);

        if (X->squaring)
            sd_fft_ctx_point_sqr(Q, abuf, m, X->depth);
        else
            sd_fft_ctx_point_mul(Q, abuf, bbuf, m, X->depth);

        sd_ifft_trunc(Q, abuf, X->depth, X->ztrunc);
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
        crt_data_struct* Rcrts, ulong * carry_out,
        nmod_t mod);
    ulong carry_out[3];
} s2worker_struct;

static void s2worker_func(void* varg)
{
    s2worker_struct* X = (s2worker_struct*) varg;

    X->f(X->z, X->zl, X->start_zi, X->stop_zi, X->ffts + X->offset, X->buf,
         X->stride, X->crts + X->offset, X->carry_out, X->mod);
}

static void _radix_mul_mpn_ctx(
    ulong * z, ulong zl, ulong zh,
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
    int write_carry_out = 0;

    an = FLINT_MIN(an, zh);
    bn = FLINT_MIN(bn, zh);

//    flint_printf("have a = %{ulong*} b = %{ulong*} mod %wu\n", a, an, b, bn, mod.n);

    FLINT_ASSERT(an > 0);
    FLINT_ASSERT(bn > 0);

    if (zl >= zh)
        return;

    if (zh == an + bn)
        write_carry_out = 1;

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

    if (bn < 1500 || radix_mulmid_force_threading)
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

    if (np*zn > 10000 || radix_mulmid_force_threading)
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

    for (i = 0; i < nthreads; i++)
    {
        slong stop = s2args[i].stop_zi;

        /* Final segment */
        if (i == nthreads - 1)
        {
            /* Write carry to the top limb of the full product. */
            if (write_carry_out)
            {
                z[stop - zl] = s2args[i].carry_out[0];

                FLINT_ASSERT(s2args[i].carry_out[1] == 0 && s2args[i].carry_out[2] == 0);
            }
            /* Otherwise, the product was truncated. */
        }
        else
        {
            /* Propagate carries from the previous segment */
            slong stop2 = s2args[i + 1].stop_zi;

            // flint_printf("overhang %wd / %wd: stop = %wd  stop2 = %wd,  zh = %wd\n", i, nthreads, stop, stop2, zh);

            ulong cy[3];
            slong j;

            cy[0] = s2args[i].carry_out[0];
            cy[1] = s2args[i].carry_out[1];
            cy[2] = s2args[i].carry_out[2];

            for (j = stop; j < stop2; j++)
            {
                add_sssaaaaaa(cy[2], cy[1], cy[0], cy[2], cy[1], cy[0], 0, 0, z[j - zl]);

                if (mod.norm == 0)
                    z[j - zl] = flint_mpn_divrem_3_1_preinv_norm(cy, cy, mod.n, mod.ninv);
                else
                    z[j - zl] = flint_mpn_divrem_3_1_preinv_unnorm(cy, cy, mod.n, mod.ninv, mod.norm);

                if (cy[0] == 0 && cy[1] == 0 && cy[2] == 0)
                    break;
            }

            /* Propagate carry to next segment */
            add_sssaaaaaa(s2args[i + 1].carry_out[2], s2args[i + 1].carry_out[1], s2args[i + 1].carry_out[0],
                          s2args[i + 1].carry_out[2], s2args[i + 1].carry_out[1], s2args[i + 1].carry_out[0],
                          cy[2], cy[1], cy[0]);

            // flint_printf("next: %wu %wu %wu\n", s2args[i + 1].carry_out[2], s2args[i + 1].carry_out[1], s2args[i + 1].carry_out[0]);
        }
    }
}

void
radix_mulmid_fft_small(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong zlo, slong zhi, const radix_t radix)
{
    if (an < bn)
    {
        FLINT_SWAP(nn_srcptr, a, b);
        FLINT_SWAP(slong, an, bn);
    }

    _radix_mul_mpn_ctx(res, zlo, zhi, a, an, b, bn, radix->B, get_default_mpn_ctx());
}

#else

void
radix_mulmid_fft_small(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong zlo, slong zhi, const radix_t radix)
{
    flint_throw(FLINT_ERROR, "radix_mulmid_fft_small: fft_small not available");
}

#endif

