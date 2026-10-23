/*
    Copyright (C) 2024, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/* res += {t, m} mod 2^(FLINT_BITS * s), where s is one of m - 1, m, m + 1 */
FLINT_FORCE_INLINE void
_add_prod(nn_ptr res, nn_srcptr t, slong m, slong s)
{
    if (s == m)
        mpn_add_n(res, res, t, m);
    else if (s == m + 1)
        res[m] += mpn_add_n(res, res, t, m);
    else
        mpn_add_n(res, res, t, m - 1);
}

void
_flint_mpn_dot_rev_generic_signed(nn_ptr res, nn_srcptr a, slong n1, nn_srcptr b, slong n2, slong len, slong s, nn_ptr t)
{
    slong k, m = n1 + n2;
    int sa, sb;

    FLINT_ASSERT(n1 >= n2);
    FLINT_ASSERT(n2 >= 1);
    FLINT_ASSERT(len >= 1);
    FLINT_ASSERT(s >= m - 1 && s <= m + 1);

    b += (len - 1) * n2;

    flint_mpn_zero(res, s);

    for (k = 0; k < len; k++)
    {
        /* a * b = A * B - sa * B * 2^(64 n1) - sb * A * 2^(64 n2) + sa * sb * 2^(64 m) */
        flint_mpn_mul(t, a, n1, b, n2);
        _add_prod(res, t, m, s);

        sa = ((slong) a[n1 - 1]) < 0;
        sb = ((slong) b[n2 - 1]) < 0;

        if (sa)
            mpn_sub(res + n1, res + n1, s - n1, b, FLINT_MIN(n2, s - n1));
        if (sb)
            mpn_sub(res + n2, res + n2, s - n2, a, FLINT_MIN(n1, s - n2));
        if (sa && sb && s == m + 1)
            res[m] += 1;

        a += n1;
        b -= n2;
    }
}
