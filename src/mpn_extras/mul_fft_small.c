/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpn_extras.h"

#include "fft.h"

/* FLINT's FFT can beat GMP below this threshold but apparently
   not consistently. Something needs retuning? */
#define FLINT_FFT_MUL_THRESHOLD 32000

#define FLINT_FFT_SMALL_MUL_THRESHOLD 400
#define FLINT_FFT_SMALL_SQR_THRESHOLD 800

#ifdef FLINT_HAVE_FFT_SMALL

#include "fft_small.h"

FLINT_TLS_PREFIX mpn_ctx_t default_mpn_ctx;
FLINT_TLS_PREFIX int default_mpn_ctx_initialized = 0;

void
mpn_ctx_cleanup(void)
{
    if (default_mpn_ctx_initialized)
    {
        default_mpn_ctx_initialized = 0;
        mpn_ctx_clear(default_mpn_ctx);
    }
}

#include "nmod.h"
#include "nmod_poly.h"

void
flint_nmod_poly_mul_mid_fft_small(mp_ptr res, slong zl, slong zh, mp_srcptr a, slong an, mp_srcptr b, slong bn, nmod_t mod)
{
    if (!default_mpn_ctx_initialized)
    {
        mpn_ctx_init(default_mpn_ctx, UWORD(0x0003f00000000001));
        flint_register_cleanup_function(mpn_ctx_cleanup);
        default_mpn_ctx_initialized = 1;
    }

    _nmod_poly_mul_mid_mpn_ctx(res, zl, zh, a, an, b, bn, mod, default_mpn_ctx);
}


mp_limb_t flint_mpn_mul_large(mp_ptr r1, mp_srcptr i1, mp_size_t n1,
                        mp_srcptr i2, mp_size_t n2)
{
    if (n2 < FLINT_FFT_SMALL_MUL_THRESHOLD || (i1 == i2 && n1 == n2 && n2 < FLINT_FFT_SMALL_SQR_THRESHOLD))
    {
        if (n1 == n2)
            if (i1 == i2)
                mpn_sqr(r1, i1, n1);
            else
                mpn_mul_n(r1, i1, i2, n2);
        else
            mpn_mul(r1, i1, n1, i2, n2);
    }
    else
    {
        if (!default_mpn_ctx_initialized)
        {
            mpn_ctx_init(default_mpn_ctx, UWORD(0x0003f00000000001));
            flint_register_cleanup_function(mpn_ctx_cleanup);
            default_mpn_ctx_initialized = 1;
        }

        mpn_ctx_mpn_mul(default_mpn_ctx, r1, i1, n1, i2, n2);
    }

    return r1[n1 + n2 - 1];
}

#else

mp_limb_t flint_mpn_mul_large(mp_ptr r1, mp_srcptr i1, mp_size_t n1,
                        mp_srcptr i2, mp_size_t n2)
{
    if (n2 < FLINT_FFT_MUL_THRESHOLD)
    {
        if (n1 == n2)
            if (i1 == i2)
                mpn_sqr(r1, i1, n1);
            else
                mpn_mul_n(r1, i1, i2, n2);
        else
            mpn_mul(r1, i1, n1, i2, n2);
    }
    else
    {
        flint_mpn_mul_fft_main(r1, i1, n1, i2, n2);
    }

    return r1[n1 + n2 - 1];
}

#endif
