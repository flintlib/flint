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
    mp_ptr q, d1;
    TMP_INIT;
    TMP_START;

    d1 = TMP_ALLOC((n + (n + 2)) * sizeof(mp_limb_t));
    q = d1 + n;

    if (mpn_add_1(d1, d, n, 1)) /* check for d + 1 == 0 */
    {
        flint_mpn_zero(dinv, n);
    }
    else
    {
        /* dinv = floor(B^(2n) / (d + 1)) mod B^n */
        flint_mpn_inv(q, d1, n, 2 * n);
        flint_mpn_copyi(dinv, q, n);
    }

    TMP_END;
}
