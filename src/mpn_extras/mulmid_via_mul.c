/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/*
    Windowed middle product via the full product.  Computes the whole product
    a * b with the standard (unbalanced) flint_mpn_mul and copies out the limb
    window [zlo, zhi).  Correct for any an >= 1, bn >= 1, 0 <= zlo < zhi <= an+bn;
    the result is the *exact* window (every partial product and carry included),
    which is a valid instance of the lower-approximation contract with zero
    deficit.  Cheapest when the window is a large fraction of the product; for a
    small window it wastes work on the discarded limbs.
*/
void
flint_mpn_mulmid_via_mul(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn,
                         mp_size_t zlo, mp_size_t zhi)
{
    mp_size_t zn = zhi - zlo;
    mp_srcptr xp, yp;
    mp_size_t xn, yn;
    mp_ptr t;
    TMP_INIT;

    FLINT_ASSERT(an >= 1 && bn >= 1);
    FLINT_ASSERT(0 <= zlo && zlo < zhi && zhi <= an + bn);

    if (an >= bn) { xp = a; xn = an; yp = b; yn = bn; }
    else          { xp = b; xn = bn; yp = a; yn = an; }

    if (zlo == 0 && zhi == an + bn)     /* window is the entire product */
    {
        flint_mpn_mul(z, xp, xn, yp, yn);
        return;
    }

    TMP_START;
    t = TMP_ARRAY_ALLOC(an + bn, mp_limb_t);
    flint_mpn_mul(t, xp, xn, yp, yn);
    flint_mpn_copyi(z, t + zlo, zn);
    TMP_END;
}
