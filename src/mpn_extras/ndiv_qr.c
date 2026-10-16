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
    Division with rounding to nearest, ties to even: q = round(a/b) written
    to (q, an - bn + 2) (the top limb being 0 or 1), and |a - q b| <= b/2
    written to (r, bn). Returns the sign (-1, 0 or 1) of a - q b.

    With the truncated quotient t and remainder u = a - t b, we have q = t + 1
    iff 2u > b, or 2u = b and t is odd.
*/
static int
_round_up(mp_srcptr q, mp_srcptr u, mp_srcptr b, mp_size_t bn)
{
    int cmp;

    /* compare 2u with b */
    if (u[bn - 1] >> (FLINT_BITS - 1))
    {
        cmp = 1;
    }
    else
    {
        mp_ptr u2;
        TMP_INIT;
        TMP_START;
        u2 = TMP_ALLOC(bn * sizeof(mp_limb_t));
        mpn_lshift(u2, u, bn, 1);
        cmp = mpn_cmp(u2, b, bn);
        TMP_END;
    }

    return cmp > 0 || (cmp == 0 && (q[0] & 1));
}

int
flint_mpn_ndiv_qr(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
{
    mp_size_t n = an - bn + 1;

    _flint_mpn_tdiv_qr(q, r, a, an, b, bn);

    if (flint_mpn_zero_p(r, bn))
    {
        q[n] = 0;
        return 0;
    }

    if (_round_up(q, r, b, bn))
    {
        q[n] = mpn_add_1(q, q, n, 1);
        mpn_sub_n(r, b, r, bn);
        return -1;
    }

    q[n] = 0;
    return 1;
}

void
flint_mpn_ndiv_q(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
{
    mp_ptr r;
    TMP_INIT;
    TMP_START;
    r = TMP_ALLOC(bn * sizeof(mp_limb_t));
    flint_mpn_ndiv_qr(q, r, a, an, b, bn);
    TMP_END;
}

int
flint_mpn_ndiv_r(mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
{
    mp_ptr q;
    int sgn;
    TMP_INIT;
    TMP_START;
    q = TMP_ALLOC((an - bn + 2) * sizeof(mp_limb_t));
    sgn = flint_mpn_ndiv_qr(q, r, a, an, b, bn);
    TMP_END;
    return sgn;
}
