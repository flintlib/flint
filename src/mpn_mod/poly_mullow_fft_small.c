/*
    Copyright (C) 2022 Daniel Schultz
    Copyright (C) 2023, 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_mod.h"

#if FLINT_HAVE_FFT_SMALL

#include "thread_pool.h"
#include "thread_support.h"
#include "nmod.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "crt_helpers.h"
#include "fft_small.h"

static ulong
nmod_set_mpn_2(nn_srcptr ad, nmod_t mod)
{
    ulong r = 0;
    NMOD_RED2(r, r, ad[1], mod);
    NMOD_RED2(r, r, ad[0], mod);
    return r;
}

static ulong
nmod_set_mpn_3(nn_srcptr ad, nmod_t mod)
{
    ulong r = 0;
    NMOD_RED2(r, r, ad[2], mod);
    NMOD_RED2(r, r, ad[1], mod);
    NMOD_RED2(r, r, ad[0], mod);
    return r;
}

/* todo: precomputed inverse */
static ulong
nmod_set_mpn(nn_srcptr ad, slong an, nmod_t mod)
{
    return mpn_mod_1(ad, an, mod.n);
}

static void _mod(
    double* abuf, ulong atrunc,
    nn_srcptr a, ulong an,
    slong nlimbs,
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
        abuf[i] = 0;
}

typedef struct {
    slong nlimbs;
    gr_ctx_struct * ctx;
} _conv_params_struct;

#define CRT_FN(NP) CAT3(_crt, NP, fn)
#define CRT_Z_TYPE nn_ptr
#define CRT_HEAD \
    const _conv_params_struct* par = (const _conv_params_struct*) params; \
    slong nlimbs = par->nlimbs;
#define CRT_EMIT(zi, r, NP, N, M) \
    { \
        ulong nn; \
        FLINT_ASSERT(mpn_cmp(r, Mhalf, N) <= 0); /* coefficients must be unsigned */ \
        nn = N; \
        MPN_NORM(r, nn); \
        mpn_mod_set_mpn(z + (zi) * nlimbs, r, nn, par->ctx); \
    }
#define CRT_TAIL
CRT_DEFINE(2, 2, 1)  /* 100 bits (unused) */
CRT_DEFINE(3, 3, 2)  /* 150 bits */
CRT_DEFINE(4, 4, 3)  /* 200 bits */
CRT_DEFINE(5, 4, 4)  /* 250 bits */
CRT_DEFINE(6, 5, 4)  /* 300 bits */
CRT_DEFINE(7, 6, 5)  /* 350 bits */
CRT_DEFINE(8, 7, 6)  /* 400 bits */
#undef CRT_FN
#undef CRT_Z_TYPE
#undef CRT_HEAD
#undef CRT_EMIT
#undef CRT_TAIL
#undef CRT_DEFINE

/* 50 bits (unused) */
static void _crt_1(
    nn_ptr FLINT_UNUSED(z), ulong FLINT_UNUSED(zl), ulong FLINT_UNUSED(zi_start), ulong FLINT_UNUSED(zi_stop),
    sd_fft_ctx_struct* FLINT_UNUSED(Rffts), double* FLINT_UNUSED(d), ulong FLINT_UNUSED(dstride),
    crt_data_struct* FLINT_UNUSED(Rcrts), slong FLINT_UNUSED(nlimbs), gr_ctx_t FLINT_UNUSED(ctx))
{
    flint_abort();
}


/* standalone forward transform of a coefficient vector into an operand,
   for the transformed polynomial representation */
void fft_small_fft_mpn_mod(fft_small_op_t X, nn_srcptr a, ulong an,
                    ulong itrunc, gr_ctx_t ctx, const fft_small_plan_t P)
{
    slong nlimbs = MPN_MOD_CTX_NLIMBS(ctx);
    ulong i;

    FLINT_ASSERT(itrunc <= P->ztrunc);

    for (i = 0; i < P->np; i++)
    {
        sd_fft_ctx_struct * Q = P->ffts + P->offset + i;
        double * d = X->data + P->stride * i;

        _mod(d, itrunc, a, an, nlimbs, Q);
        sd_fft_trunc(Q, d, P->depth, itrunc, P->ztrunc);
    }

    X->itrunc = itrunc;
    X->domain = FFT_SMALL_OP_PRIMAL;
}

/* chinese remaindering of an inverse-transformed operand restricted to the
   window [zl, zh) inside [P->zl, P->zh); z receives zh - zl coefficients */

static fft_small_crt_func _mpn_mod_crt_fn(ulong np);

typedef struct {
    nn_ptr z;
    const fft_small_op_struct * X;
    _conv_params_struct * par;
    const fft_small_plan_struct * P;
    fft_small_crt_func crt_fn;
    ulong start_zi;
    ulong stop_zi;
} _op_export_mpn_mod_worker_struct;

static void _op_export_mpn_mod_worker_func(void * varg)
{
    _op_export_mpn_mod_worker_struct * W =
        (_op_export_mpn_mod_worker_struct *) varg;
    const fft_small_plan_struct * P = W->P;

    W->crt_fn((void *) W->z, P->zl, W->start_zi, W->stop_zi,
              P->ffts + P->offset, W->X->data, P->stride,
              P->crts + P->offset, P->bound_c, NULL, W->par);
}

void fft_small_export_mpn_mod_range(nn_ptr z, const fft_small_op_t X,
                    ulong zl, ulong zh, gr_ctx_t ctx, const fft_small_plan_t P)
{
    slong nlimbs = MPN_MOD_CTX_NLIMBS(ctx);
    ulong i, o;
    thread_pool_handle * handles = NULL;
    slong nworkers = 0;
    ulong nthreads;
    _conv_params_struct par;
    _op_export_mpn_mod_worker_struct args[8];

    FLINT_ASSERT(P->zl <= zl && zh <= P->zh);

    if (zl >= zh)
        return;

    par.nlimbs = nlimbs;
    par.ctx = ctx;

    if (P->np * (zh - zl) > 10000)
        nworkers = flint_request_threads(&handles, 8);
    nthreads = nworkers + 1;

    o = zl;
    for (i = 0; i < nthreads; i++)
    {
        _op_export_mpn_mod_worker_struct * W = args + i;
        W->z = z - (zl - P->zl) * nlimbs;   /* worker offsets by P->zl */
        W->X = X;
        W->par = &par;
        W->P = P;
        W->crt_fn = _mpn_mod_crt_fn(P->np);
        W->start_zi = o;
        ulong newo = n_round_down(zl + (i + 1) * (zh - zl) / nthreads, BLK_SZ);
        o = i + 1 < nthreads ? FLINT_MAX(o, newo) : zh;
        W->stop_zi = o;
    }

    for (i = nworkers; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0,
                         _op_export_mpn_mod_worker_func, args + i);
    _op_export_mpn_mod_worker_func(args + 0);
    for (i = nworkers; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);

    flint_give_back_threads(handles, nworkers);
}

/* ------------------------------------------------------------------------ */
/* glue for the generic convolution engine                                  */
/* ------------------------------------------------------------------------ */

static void _mod_fn(double* xbuf, ulong xtrunc,
    const void* x, ulong xn, slong FLINT_UNUSED(xaux),
    const sd_fft_ctx_struct* fft, const void* params)
{
    const _conv_params_struct* par = (const _conv_params_struct*) params;

    _mod(xbuf, xtrunc, (nn_srcptr) x, xn, par->nlimbs, fft);
}

#define DEFINE_IT(NP) \
static void CAT3(_crt, NP, fn)(void* z, ulong zl, \
    ulong zi_start, ulong zi_stop, \
    sd_fft_ctx_struct* Rffts, double* d, ulong dstride, \
    crt_data_struct* Rcrts, ulong FLINT_UNUSED(min_an_bn), \
    ulong* FLINT_UNUSED(local), const void* params) \
{ \
    const _conv_params_struct* par = (const _conv_params_struct*) params; \
 \
    CAT(_crt, NP)((nn_ptr) z, zl, zi_start, zi_stop, Rffts, d, dstride, \
                  Rcrts, par->nlimbs, par->ctx); \
}
DEFINE_IT(1)
#undef DEFINE_IT

static fft_small_crt_func _mpn_mod_crt_fn(ulong np)
{
    return np == 1 ? _crt_1_fn :
           np == 2 ? _crt_2_fn :
           np == 3 ? _crt_3_fn :
           np == 4 ? _crt_4_fn :
           np == 5 ? _crt_5_fn :
           np == 6 ? _crt_6_fn :
           np == 7 ? _crt_7_fn :
                     _crt_8_fn;
}

/* currently the same thresholds as the fmpz driver, but kept separate so
   that they can be tuned independently (the per-coefficient conversions
   here are more expensive) */
static ulong _mpn_mod_conv_s1_threads(ulong np, ulong tune_n)
{
    return ((np >= 2 && tune_n >= 1000) || (np >= 4 && tune_n >= 300)) ? np : 1;
}

static int _mpn_mod_conv_s1_b_worker(ulong np, ulong tune_n, int squaring)
{
    return !squaring && ((np == 1 && tune_n > 5000) ||
                         (np >= 2 && tune_n >= 1000));
}

static int _mpn_mod_conv_s2_rethread(ulong np, ulong zn)
{
    return zn > 50000 || (np >= 2 && zn > 20000) || (np >= 4 && zn > 800);
}

static const fft_small_conv_tuning _mpn_mod_conv_tuning =
    { _mpn_mod_conv_s1_threads, _mpn_mod_conv_s1_b_worker,
      _mpn_mod_conv_s2_rethread };

static int _mpn_mod_poly_mulmid_fft_small_internal(nn_ptr z, ulong zl, ulong zh,
    nn_srcptr a, ulong an,
    nn_srcptr b, ulong bn,
    mpn_ctx_t R, gr_ctx_t ctx)
{
    ulong modbits;
    ulong zn = an + bn - 1;
    ulong atrunc, btrunc;
    int squaring;
    slong bits1, bits2;
    slong nlimbs = MPN_MOD_CTX_NLIMBS(ctx);
    flint_bitcnt_t nbits = MPN_MOD_CTX_MODULUS_BITS(ctx);
    fft_small_plan_t P;
    fft_small_conv_arg_struct A;
    _conv_params_struct par;

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

    /* should be +sign instead of +1, but currently the CRT code doesn't
       distinguish between the signed and unsigned cases */
    modbits = FLINT_ABS(bits1) + FLINT_ABS(bits2) + 1;

    FLINT_ASSERT(zl < zh);
    FLINT_ASSERT(zh <= zn);

    atrunc = n_round_up(an, BLK_SZ);
    btrunc = n_round_up(bn, BLK_SZ);

    P->R = R;
    P->sign = 0;

    _fft_small_plan_set_window(P, zl, zh, zn, n_max(atrunc, btrunc));

    /* need prod_of_primes >= bn * 2^modbits */
    if (!_fft_small_plan_set_bound(P, bn, modbits, MPN_CTX_NCRTS))
        return GR_UNABLE;

    _fft_small_plan_set_normalizers(P);

    par.nlimbs = nlimbs;
    par.ctx = ctx;

    A.a = a; A.an = an; A.aaux = 0; A.atrunc = atrunc;
    A.b = b; A.bn = bn; A.baux = 0; A.btrunc = btrunc;
    A.bfft = NULL;
    A.bfft_stride = 0;
    A.squaring = squaring;
    A.tune_n = bn;
    A.min_an_bn = FLINT_MIN(an, bn);
    A.mod_fn = _mod_fn;
    A.s2_finish = NULL;
    A.crt_fn = _mpn_mod_crt_fn(P->np);
    A.params = &par;
    A.tuning = &_mpn_mod_conv_tuning;

    _fft_small_conv(z, P, &A);

    fft_small_plan_clear(P);

    return GR_SUCCESS;
}


int
_mpn_mod_poly_mullow_fft_small(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong len, gr_ctx_t ctx)
{
    if (len1 >= len2)
        return _mpn_mod_poly_mulmid_fft_small_internal(res, 0, len, poly1, len1, poly2, len2, get_default_mpn_ctx(), ctx);
    else
        return _mpn_mod_poly_mulmid_fft_small_internal(res, 0, len, poly2, len2, poly1, len1, get_default_mpn_ctx(), ctx);
}

int
_mpn_mod_poly_mulmid_fft_small(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, gr_ctx_t ctx)
{
    if (len1 >= len2)
        return _mpn_mod_poly_mulmid_fft_small_internal(res, nlo, nhi, poly1, len1, poly2, len2, get_default_mpn_ctx(), ctx);
    else
        return _mpn_mod_poly_mulmid_fft_small_internal(res, nlo, nhi, poly2, len2, poly1, len1, get_default_mpn_ctx(), ctx);
}

#else /* FLINT_HAVE_FFT_SMALL */

int
_mpn_mod_poly_mullow_fft_small(nn_ptr res,  nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong len, gr_ctx_t ctx)
{
    return GR_UNABLE;
}

int
_mpn_mod_poly_mulmid_fft_small(nn_ptr res, nn_srcptr poly1, slong len1, nn_srcptr poly2, slong len2, slong nlo, slong nhi, gr_ctx_t ctx)
{
    return GR_UNABLE;
}

#endif
