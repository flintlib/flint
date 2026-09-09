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
    Ceiling division: q = ceil(a/b) written to (q, an - bn + 2) (the top
    limb being 0 or 1), and r = q b - a in [0, b) written to (r, bn).
*/
void
flint_mpn_cdiv_qr(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
{
    mp_size_t n = an - bn + 1;

    _flint_mpn_tdiv_qr(q, r, a, an, b, bn);

    if (flint_mpn_zero_p(r, bn))
    {
        q[n] = 0;
    }
    else
    {
        q[n] = mpn_add_1(q, q, n, 1);
        mpn_sub_n(r, b, r, bn);
    }
}

void
flint_mpn_cdiv_q(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
{
    mp_size_t n = an - bn + 1;
    mp_ptr r;
    TMP_INIT;
    TMP_START;
    r = TMP_ALLOC(bn * sizeof(mp_limb_t));

    _flint_mpn_tdiv_qr(q, r, a, an, b, bn);

    if (flint_mpn_zero_p(r, bn))
        q[n] = 0;
    else
        q[n] = mpn_add_1(q, q, n, 1);

    TMP_END;
}

void
flint_mpn_cdiv_r(mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
{
    flint_mpn_tdiv_r(r, a, an, b, bn);

    if (!flint_mpn_zero_p(r, bn))
        mpn_sub_n(r, b, r, bn);
}
