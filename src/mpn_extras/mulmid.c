/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "ulong_extras.h"
#include "flint-mparam.h"

#define APPROX_EQUAL(mm, nn) (FLINT_ABS((mm) - (nn)) <= FLINT_MAX((mm), (nn)) / 8)

void
flint_mpn_mulmid(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn,
                 mp_size_t zlo, mp_size_t zhi)
{
    FLINT_ASSERT(an >= 1);
    FLINT_ASSERT(bn >= 1);
    FLINT_ASSERT(zlo >= 0);
    FLINT_ASSERT(zhi > zlo);
    FLINT_ASSERT(zhi <= an + bn);

    if (zlo == 0 && zhi == an + bn)
    {
        if (a == b && an == bn)
            flint_mpn_sqr(z, a, an);
        else if (an >= bn)
            flint_mpn_mul(z, a, an, b, bn);
        else
            flint_mpn_mul(z, b, bn, a, an);
        return;
    }

    an = FLINT_MIN(an, zhi);
    bn = FLINT_MIN(bn, zhi);

    if (FLINT_MIN(FLINT_MIN(an, bn), zhi - zlo) < 48)
    {
        flint_mpn_mulmid_classical(z, a, an, b, bn, zlo, zhi);
        return;
    }

#if FLINT_HAVE_FFT_SMALL
    if (FLINT_MIN(an, bn) >= ((a == b) ? FLINT_FFT_SMALL_SQR_THRESHOLD : FLINT_FFT_SMALL_MUL_THRESHOLD))
    {
        flint_mpn_mulmid_fft_small(z, a, an, b, bn, zlo, zhi);
        return;
    }
#endif

    if (zlo <= (an + bn) / 8 && FLINT_MIN(an, bn) >= 7 * zhi / 8)
    {
        flint_mpn_mulmid_via_mullow_n(z, a, an, b, bn, zlo, zhi);
        return;
    }

    if (zhi >= 7 * (an + bn) / 8 && FLINT_MIN(an, bn) >= 7 * (zhi - zlo) / 8)
    {
        flint_mpn_mulmid_via_mulhigh_n(z, a, an, b, bn, zlo, zhi);
        return;
    }

#if FLINT_HAVE_NATIVE_mpn_mulmid_n
    if (APPROX_EQUAL(an, 2 * bn) && APPROX_EQUAL(zlo, bn) && APPROX_EQUAL(zhi, an))
    {
        flint_mpn_mulmid_via_n_padded(z, a, an, b, bn, zlo, zhi);
        return;
    }

    if (APPROX_EQUAL(bn, 2 * an) && APPROX_EQUAL(zlo, an) && APPROX_EQUAL(zhi, bn))
    {
        flint_mpn_mulmid_via_n_padded(z, b, bn, a, an, zlo, zhi);
        return;
    }
#endif

    flint_mpn_mulmid_via_mul(z, a, an, b, bn, zlo, zhi);
    return;
}

