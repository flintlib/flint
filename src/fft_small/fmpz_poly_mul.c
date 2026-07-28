/*
    Copyright (C) 2022 Daniel Schultz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "thread_pool.h"
#include "thread_support.h"
#include "nmod.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "crt_helpers.h"
#include "fft_small.h"

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
                    aa[0] = a[i + j + 0] >= 0 ? (ulong) a[i + j + 0] : a[i + j + 0] + p;
                    aa[1] = a[i + j + 1] >= 0 ? (ulong) a[i + j + 1] : a[i + j + 1] + p;
                    aa[2] = a[i + j + 2] >= 0 ? (ulong) a[i + j + 2] : a[i + j + 2] + p;
                    aa[3] = a[i + j + 3] >= 0 ? (ulong) a[i + j + 3] : a[i + j + 3] + p;
                    aa[4] = a[i + j + 4] >= 0 ? (ulong) a[i + j + 4] : a[i + j + 4] + p;
                    aa[5] = a[i + j + 5] >= 0 ? (ulong) a[i + j + 5] : a[i + j + 5] + p;
                    aa[6] = a[i + j + 6] >= 0 ? (ulong) a[i + j + 6] : a[i + j + 6] + p;
                    aa[7] = a[i + j + 7] >= 0 ? (ulong) a[i + j + 7] : a[i + j + 7] + p;

                    vec8n t = vec8n_load_unaligned(aa);
                    FLINT_ASSERT(i+j < atrunc);
                    vec8d_store_aligned(aI + j, vec8n_convert_limited_vec8d(t));
                }
            }

            aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
            for (j = 0; j < an - i; j++)
                aI[j] = a[i + j] >= 0 ? (ulong) a[i + j] : a[i + j] + p;
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
        abuf[i] = 0;
}


static void fmpz_neg_ui_array(fmpz_t out, const ulong * in, slong in_len)
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
        mpz_ptr mpz = _fmpz_promote(out);
        mp_ptr mp = FLINT_MPZ_REALLOC(mpz, size);
        mpz->_mp_size = -size;
        flint_mpn_copyi(mp, in, size);
    }
}

#define CRT_FN(NP) CAT3(_crt, NP, fn)
#define CRT_Z_TYPE fmpz*
#define CRT_HEAD
#define CRT_EMIT(zi, r, NP, N, M) \
    if (mpn_cmp(r, Mhalf, N) > 0) \
    { \
        CAT(multi_rsub, N)(r, crt_data_prod_primes(Rcrts + np - 1)); \
        fmpz_neg_ui_array(&z[zi], r, N); \
    } \
    else \
        fmpz_set_ui_array(&z[zi], r, N);
#define CRT_TAIL
CRT_DEFINE(2, 2, 1)  /* 100 bits */
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

/* 50 bits */
static void _crt_1(
    fmpz * z, ulong zl, ulong zi_start, ulong zi_stop,
    sd_fft_ctx_struct* Rffts, double* d, ulong dstride,
    crt_data_struct* FLINT_UNUSED(Rcrts))
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


/* ------------------------------------------------------------------------ */
/* glue for the generic convolution engine                                  */
/* ------------------------------------------------------------------------ */

static void _mod_fn(double* xbuf, ulong xtrunc,
    const void* x, ulong xn, slong xaux,
    const sd_fft_ctx_struct* fft, const void* FLINT_UNUSED(params))
{
    _mod(xbuf, xtrunc, (const fmpz *) x, xn, xaux, fft);
}

#define DEFINE_IT(NP) \
static void CAT3(_crt, NP, fn)(void* z, ulong zl, \
    ulong zi_start, ulong zi_stop, \
    sd_fft_ctx_struct* Rffts, double* d, ulong dstride, \
    crt_data_struct* Rcrts, ulong FLINT_UNUSED(min_an_bn), \
    ulong* FLINT_UNUSED(local), const void* FLINT_UNUSED(params)) \
{ \
    CAT(_crt, NP)((fmpz *) z, zl, zi_start, zi_stop, Rffts, d, dstride, \
                  Rcrts); \
}
DEFINE_IT(1)
#undef DEFINE_IT

static fft_small_crt_func _fmpz_crt_fn(ulong np)
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

static ulong _fmpz_conv_s1_threads(ulong np, ulong tune_n)
{
    return ((np >= 2 && tune_n >= 1000) || (np >= 4 && tune_n >= 300)) ? np : 1;
}

static int _fmpz_conv_s1_b_worker(ulong np, ulong tune_n, int squaring)
{
    return !squaring && ((np == 1 && tune_n > 5000) ||
                         (np >= 2 && tune_n >= 1000));
}

static int _fmpz_conv_s2_rethread(ulong np, ulong zn)
{
    return zn > 50000 || (np >= 2 && zn > 20000) || (np >= 4 && zn > 800);
}

static const fft_small_conv_tuning _fmpz_conv_tuning =
    { _fmpz_conv_s1_threads, _fmpz_conv_s1_b_worker, _fmpz_conv_s2_rethread };


int _fmpz_poly_mul_mid_mpn_ctx(
    fmpz * z, ulong zl, ulong zh,
    const fmpz * a, ulong an,
    const fmpz * b, ulong bn,
    mpn_ctx_t R)
{
    ulong modbits;
    ulong zn = an + bn - 1;
    ulong atrunc, btrunc;
    int squaring;
    slong bits1, bits2;
    int sign;
    fft_small_plan_t P;
    fft_small_conv_arg_struct A;

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
    /* should be +sign instead of +1, but currently the CRT code doesn't
       distinguish between the signed and unsigned cases */
    modbits = FLINT_ABS(bits1) + FLINT_ABS(bits2) + 1;

    FLINT_ASSERT(zl < zh);
    FLINT_ASSERT(zh <= zn);

    atrunc = n_round_up(an, BLK_SZ);
    btrunc = n_round_up(bn, BLK_SZ);

    P->R = R;
    P->sign = sign;

    _fft_small_plan_set_window(P, zl, zh, zn, n_max(atrunc, btrunc));

    /* need prod_of_primes >= bn * 2^modbits */
    if (!_fft_small_plan_set_bound(P, bn, modbits, MPN_CTX_NCRTS))
        return 0;

    _fft_small_plan_set_normalizers(P);

    A.a = a; A.an = an; A.aaux = bits1; A.atrunc = atrunc;
    A.b = b; A.bn = bn; A.baux = bits2; A.btrunc = btrunc;
    A.bfft = NULL;
    A.bfft_stride = 0;
    A.squaring = squaring;
    A.tune_n = bn;
    A.min_an_bn = FLINT_MIN(an, bn);
    A.mod_fn = _mod_fn;
    A.s2_finish = NULL;
    A.crt_fn = _fmpz_crt_fn(P->np);
    A.params = NULL;
    A.tuning = &_fmpz_conv_tuning;

    _fft_small_conv(z, P, &A);

    fft_small_plan_clear(P);

    return 1;
}
