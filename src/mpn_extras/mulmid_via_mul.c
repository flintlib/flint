/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

void
flint_mpn_mulmid_via_mul(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn,
                         mp_size_t zlo, mp_size_t zhi)
{
    mp_size_t zn = zhi - zlo;
    mp_ptr t;
    TMP_INIT;

    FLINT_ASSERT(an >= 1 && bn >= 1);
    FLINT_ASSERT(0 <= zlo && zlo < zhi && zhi <= an + bn);

    TMP_START;
    t = TMP_ARRAY_ALLOC(an + bn, mp_limb_t);

    if (a == b && an == bn)
        flint_mpn_sqr(t, a, an);
    else if (an >= bn)
        flint_mpn_mul(t, a, an, b, bn);
    else
        flint_mpn_mul(t, b, bn, a, an);

    flint_mpn_copyi(z, t + zlo, zn);
    TMP_END;
}
