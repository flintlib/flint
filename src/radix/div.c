/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "radix.h"

/* todo: optimise */
ulong
radix_divrem_1(nn_ptr res, nn_srcptr a, slong an, ulong d, const radix_t radix)
{
    slong i;
    ulong q, r, hi, lo, B = LIMB_RADIX(radix);

    q = a[an - 1] / d;
    r = a[an - 1] % d;
    res[an - 1] = q;

    for (i = an - 2; i >= 0; i--)
    {
        umul_ppmm(hi, lo, r, B);
        add_ssaaaa(hi, lo, hi, lo, 0, a[i]);
        udiv_qrnnd(q, r, hi, lo, d);
        res[i] = q;
    }

    return r;
}

/* todo: optimise */
void
radix_divexact_1(nn_ptr res, nn_srcptr a, slong an, ulong d, const radix_t radix)
{
    ulong r;
    r = radix_divrem_1(res, a, an, d, radix);
    FLINT_ASSERT(r == 0);
}

