/*
    Copyright (C) 2023 Albin Ahlbäck
    Copyright (C) 2023, 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

void __gmpn_mul_basecase(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);

#if FLINT_HAVE_ADX
/* We don't call mpn_mul_basecase directly because our own basecases
   already cover small operands. */
# define WANT_GMP_MUL_BASECASE 0
# define BASECASE_LIMIT 0
/* toom22 on top our own basecase beats GMP up to some point */
# define WANT_TOOM22 1
# define TOOM22_LIMIT 230
#elif FLINT_HAVE_ARMV8
# define WANT_GMP_MUL_BASECASE 0
# define BASECASE_LIMIT 0
# define WANT_TOOM22 1
# define TOOM22_LIMIT 230
#else
/* As fallback where we don't have our own mul_basecase: */
/* GMP's MUL_TOOM22_THRESHOLD is >= 16 on most machines */
# define WANT_GMP_MUL_BASECASE 1
# define BASECASE_LIMIT 16
# define WANT_TOOM22 0
# define TOOM22_LIMIT 0
#endif


#ifdef FLINT_HAVE_FFT_SMALL
# include "fft_small.h"
# define FFT_MUL mpn_mul_default_mpn_ctx
#else
# include "fft.h"
# define FFT_MUL flint_mpn_mul_fft_main
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

#if FLINT_MPN_MUL_FUNC_N_TAB_WIDTH && FLINT_HAVE_NATIVE_MUL_2
    if (yn == 1)
        r[xn + yn - 1] = mpn_mul_1(r, x, xn, y[0]);
    else if (yn == 2)
        return flint_mpn_mul_2(r, x, xn, y);
    else if (yn <= FLINT_MPN_MUL_FUNC_N_TAB_WIDTH)
        return FLINT_MPN_MUL_HARD(r, y, yn, x, xn);
#else
    if (WANT_GMP_MUL_BASECASE && xn <= BASECASE_LIMIT)
        __gmpn_mul_basecase(r, x, xn, y, yn);
    else if (yn == 1)
        r[xn + yn - 1] = mpn_mul_1(r, x, xn, y[0]);
#endif
    else if (WANT_TOOM22 && yn <= TOOM22_LIMIT && 5 * yn >= 4 * xn)
        flint_mpn_mul_toom22(r, x, xn, y, yn, NULL);
    else if (yn < FLINT_FFT_MUL_THRESHOLD)
        mpn_mul(r, x, xn, y, yn);
    else
        FFT_MUL(r, x, xn, y, yn);

    return r[xn + yn - 1];
}

void _flint_mpn_mul_n(mp_ptr r, mp_srcptr x, mp_srcptr y, mp_size_t n)
{
    if (WANT_GMP_MUL_BASECASE && n <= BASECASE_LIMIT)
        __gmpn_mul_basecase(r, x, n, y, n);
    else if (WANT_TOOM22 && n <= TOOM22_LIMIT)
        flint_mpn_mul_toom22(r, x, n, y, n, NULL);
    else if (n < FLINT_FFT_MUL_THRESHOLD)
        mpn_mul_n(r, x, y, n);
    else
        FFT_MUL(r, x, n, y, n);
}

mp_limb_t _flint_mpn_sqr(mp_ptr r, mp_srcptr x, mp_size_t n)
{
    /* We cannot call __gmpn_sqr_basecase directly because it
       may not support n above a certain size. */

    if (n < FLINT_FFT_SQR_THRESHOLD)
        mpn_sqr(r, x, n);
    else
    {
        /* TODO: Don't do a full multiplication here. */
        FFT_MUL(r, x, n, x, n);
    }

    return r[2 * n - 1];
}
