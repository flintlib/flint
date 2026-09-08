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
    Exact division test: if b divides a, sets (q, an - bn + 1) to the
    quotient and returns 1; otherwise returns 0 (q is then undefined).
    Requires an >= bn >= 1 and b[bn-1] != 0.
*/
int
flint_mpn_div(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
{
    mp_ptr r;
    int exact;
    TMP_INIT;

    FLINT_ASSERT(an >= bn);
    FLINT_ASSERT(bn >= 1);
    FLINT_ASSERT(b[bn - 1] != 0);

    if (bn == 1)
        return mpn_divrem_1(q, 0, a, an, b[0]) == 0;

    if (flint_mpn_zero_p(a + bn, an - bn) && mpn_cmp(a, b, bn) < 0)
    {
        if (flint_mpn_zero_p(a, bn))
        {
            flint_mpn_zero(q, an - bn + 1);
            return 1;
        }
        return 0;
    }

    TMP_START;
    r = TMP_ALLOC(bn * sizeof(mp_limb_t));
    _flint_mpn_tdiv_qr(q, r, a, an, b, bn);
    exact = flint_mpn_zero_p(r, bn);
    TMP_END;

    return exact;
}
