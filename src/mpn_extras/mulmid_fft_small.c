/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

#if FLINT_HAVE_FFT_SMALL

#include "fft_small.h"

/*
    Windowed middle product through the small-prime FFT.  This is a thin wrapper
    over _mpn_ctx_mpn_mul_range, which returns the limb window [zlo, zhi) of a*b
    as the same lower approximation as the rest of the family (partial products
    landing entirely below zlo are dropped).  The default (thread-local) context
    is used, exactly as mpn_mul_default_mpn_ctx does for the full product.
*/
void
flint_mpn_mulmid_fft_small(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn,
                           mp_size_t zlo, mp_size_t zhi)
{
    _mpn_ctx_mpn_mul_range(get_default_mpn_ctx(), z,
                           (ulong) zlo, (ulong) zhi,
                           a, (ulong) an, b, (ulong) bn);
}

#endif
