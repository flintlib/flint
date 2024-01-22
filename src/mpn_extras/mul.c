/*
    Copyright (C) 2023 Albin Ahlbäck
    Copyright (C) 2023, 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

void __gmpn_mul_basecase(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);

#ifdef FLINT_HAVE_FFT_SMALL
#include "fft_small.h"
#define FFT_MUL mpn_mul_default_mpn_ctx
#else
#include "fft.h"
#define FFT_MUL flint_mpn_mul_fft_main
#endif

mp_limb_t _flint_mpn_mul(mp_ptr r, mp_srcptr x, mp_size_t xn, mp_srcptr y, mp_size_t yn)
{
    /* Experimental: strip trailing zeros. Normally this should
       be handled by the caller where appropriate, but there can be
       situations where it helps to do so here. */
    /*
    while (xn > 1 && x[0] == 0)
    {
        xn--;
        x++;
        r[0] = 0;
        r++;
    }

    while (yn > 1 && y[0] == 0)
    {
        yn--;
        y++;
        r[0] = 0;
        r++;
    }

    if (xn < yn)
    {
        FLINT_SWAP(mp_srcptr, x, y);
        FLINT_SWAP(mp_size_t, xn, yn);
    }
    */

    /* GMP's MUL_TOOM22_THRESHOLD is >= 16 on most machines */
    if (xn <= 16)
        __gmpn_mul_basecase(r, x, xn, y, yn);
    else if (yn == 1)
        r[xn + yn - 1] = mpn_mul_1(r, x, xn, y[0]);
    else if (yn < FLINT_FFT_MUL_THRESHOLD)
        return mpn_mul(r, x, xn, y, yn);
    else
        FFT_MUL(r, x, xn, y, yn);

    return r[xn + yn - 1];
}

void _flint_mpn_mul_n(mp_ptr r, mp_srcptr x, mp_srcptr y, mp_size_t n)
{
    /* GMP's MUL_TOOM22_THRESHOLD is >= 16 on most machines */
    if (n <= 16)
        __gmpn_mul_basecase(r, x, n, y, n);
    else if (n < FLINT_FFT_MUL_THRESHOLD)
        mpn_mul_n(r, x, y, n);
    else
        FFT_MUL(r, x, n, y, n);
}

void _flint_mpn_sqr(mp_ptr r, mp_srcptr x, mp_size_t n)
{
    /* We cannot call __gmpn_sqr_basecase directly because it
       may not support n above a certain size. */

    if (n < FLINT_FFT_SQR_THRESHOLD)
        mpn_sqr(r, x, n);
    else
        FFT_MUL(r, x, n, x, n);
}
