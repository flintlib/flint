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
        mpn_mul_default_mpn_ctx(r1, i1, n1, i2, n2);
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

