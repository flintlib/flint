/*
    Copyright (C) 2026 Fredrik Johansson

    Developed using Claude Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#if FLINT_HAVE_FFT_SMALL
# include "fft_small.h"
#endif

/*
    Two products with a common operand: transformed with the transform
    of x shared when both products are long enough, two products
    otherwise.

    The rule, with T = FLINT_MPN_MUL_SHARED_THRESHOLD: the longer
    product (by max(an, bn)) must be in the FFT range of flint_mpn_mul
    with threshold T, the shorter with threshold T / 2. The shared
    transform pays as soon as the longer product is transformed, and the
    shorter product, which then costs about its share of the
    transforms, need not reach the FFT range on its own. T defaults to
    FLINT_FFT_SMALL_MUL_THRESHOLD and is set per architecture in
    flint-mparam.h where measured otherwise.

    Measured with p-mpn_mul_shared (Zen 3 with T = 240, Emerald Rapids
    with T = 240): above the rule the shared products are 7-14% faster
    than two products, in the two-prime and the wide geometry alike,
    for x as long as, longer than or a quarter of the longer partner
    and a shorter partner of any length down to an eighth of it; below
    it they were at best even, and lost up to 10% for a short partner
    of 32-48 limbs.
*/

#ifndef FLINT_MPN_MUL_SHARED_THRESHOLD
#define FLINT_MPN_MUL_SHARED_THRESHOLD FLINT_FFT_SMALL_MUL_THRESHOLD
#endif

static int
_use_fft(slong l, slong s, slong t)
{
    if (l < s)
        FLINT_SWAP(slong, l, s);

    return s >= t || (l + s > 2 * t && 3 * s > t);
}

int
_flint_mpn_mul_shared_use_fft(mp_size_t an, mp_size_t bn, mp_size_t xn)
{
#if FLINT_HAVE_FFT_SMALL
    slong t = FLINT_MPN_MUL_SHARED_THRESHOLD;

    return _use_fft(FLINT_MAX(an, bn), xn, t)
        && _use_fft(FLINT_MIN(an, bn), xn, t / 2);
#else
    return 0;
#endif
}

static void
_mul(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
{
    if (an >= bn)
        flint_mpn_mul(z, a, an, b, bn);
    else
        flint_mpn_mul(z, b, bn, a, an);
}

void
flint_mpn_mul_shared(mp_ptr z1, mp_srcptr a, mp_size_t an,
    mp_ptr z2, mp_srcptr b, mp_size_t bn, mp_srcptr x, mp_size_t xn)
{
    FLINT_ASSERT(an >= 1 && bn >= 1 && xn >= 1);

#if FLINT_HAVE_FFT_SMALL
    if (_flint_mpn_mul_shared_use_fft(an, bn, xn))
    {
        _fft_small_mpn_mul_shared(z1, a, an, z2, b, bn, x, xn);
        return;
    }
#endif

    _mul(z1, a, an, x, xn);
    _mul(z2, b, bn, x, xn);
}
