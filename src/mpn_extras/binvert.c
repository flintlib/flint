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

/* v = m^(-1) mod B^n for odd m, by Hensel lifting v <- v*(2 - m*v).
   No aliasing. About 2/3 the cost of two full-size mullows. */
void
_flint_mpn_binvert(nn_ptr v, nn_srcptr m, mp_size_t n)
{
    mp_size_t k, k2;
    nn_ptr u, w;
    TMP_INIT;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(m[0] & 1);

    v[0] = n_binvert(m[0]);

    if (n == 1)
        return;

    TMP_START;
    u = TMP_ALLOC(2 * n * sizeof(ulong));
    w = u + n;

    for (k = 1; k < n; k = k2)
    {
        k2 = FLINT_MIN(2 * k, n);

        /* v has k valid limbs; the lift reads it zero-padded to k2 */
        flint_mpn_zero(v + k, k2 - k);

        /* u = 2 - m*v mod B^k2 */
        flint_mpn_mullow_n(u, m, v, k2);
        mpn_neg(u, u, k2);
        mpn_add_1(u, u, k2, 2);

        /* v = v*u mod B^k2; low k limbs are unchanged by construction,
           but recomputing them costs nothing extra with mullow */
        flint_mpn_mullow_n(w, v, u, k2);
        flint_mpn_copyi(v, w, k2);
    }

    TMP_END;
}
