/*
    Copyright (C) 2012 William Hart
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

void flint_mpn_preinvn(mp_ptr dinv, mp_srcptr d, mp_size_t n)
{
    mp_ptr q, r, d1, tmp;
    TMP_INIT;
    TMP_START;

    tmp = TMP_ALLOC((n + (2 * n + 1) + (n + 2)) * sizeof(mp_limb_t));
    d1 = tmp;
    r = d1 + n;
    q = r + (2 * n + 1);

    if (mpn_add_1(d1, d, n, 1)) /* check for d + 1 == 0 */
    {
        flint_mpn_zero(dinv, n);
    }
    else
    {
        flint_mpn_zero(r, 2 * n);
        r[2*n] = 1;

        mpn_tdiv_qr(q, r, 0, r, 2 * n + 1, d1, n);
        flint_mpn_copyi(dinv, q, n);
    }

    TMP_END;
}
