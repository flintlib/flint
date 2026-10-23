/*
    Copyright (C) 2012 William Hart
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

mp_limb_t
flint_mpn_divrem_preinv1(mp_ptr q, mp_ptr a, mp_size_t m,
                         mp_srcptr b, mp_size_t n, mp_limb_t dinv)
{
    mp_limb_t qh;
    mp_ptr tp;
    TMP_INIT;

    FLINT_ASSERT(n >= 2);
    FLINT_ASSERT(m >= n);

    /* the scratch space is only used by the divide and conquer code */
    if (n < FLINT_MPN_DIV_DC_CUTOFF || m - n < FLINT_MPN_DIV_DC_CUTOFF)
        return _flint_mpn_divrem_preinv1(q, a, m, b, n, dinv, NULL);

    TMP_START;
    tp = TMP_ALLOC(n * sizeof(mp_limb_t));
    qh = _flint_mpn_divrem_preinv1(q, a, m, b, n, dinv, tp);
    TMP_END;
    return qh;
}
