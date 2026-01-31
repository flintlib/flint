/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"

ulong
radix_neg(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)
{
    slong i;

    for (i = 0; i < an; i++)
    {
        if (a[i] == 0)
            res[i] = 0;
        else
        {
            ulong B = LIMB_RADIX(radix);

            res[i] = B - a[i];

            for (i = i + 1; i < an; i++)
                res[i] = B - 1 - a[i];

            return 1;
        }
    }

    return 0;
}

ulong
radix_add(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix)
{
    ulong cy, hi, lo, B = LIMB_RADIX(radix);
    slong i;

    FLINT_ASSERT(bn >= 1);
    FLINT_ASSERT(an >= bn);

    cy = 0;
    /* for reasons that are not completely obvious, this runs faster with
       the true carry inverted */
    cy -= 1;

    for (i = 0; i < bn; i++)
    {
        sub_ddmmss(hi, lo, 0, a[i] + (cy + 1), 0, B - b[i]);
        res[i] = lo + (hi & B);
        cy = hi;
    }

    cy += 1;

    for (i = bn; i < an; i++)  /* improve */
    {
        res[i] = a[i] + cy;

        if (res[i] == B)
        {
            cy = 1;
            res[i] = 0;
        }
        else
        {
            cy = 0;
            if (res != a)
                flint_mpn_copyi(res + i + 1, a + i + 1, an - i - 1);
            break;
        }
    }

    return cy;
}

ulong
radix_sub(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix)
{
    ulong cy, hi, lo, B = LIMB_RADIX(radix);
    slong i;

    FLINT_ASSERT(bn >= 1);
    FLINT_ASSERT(an >= bn);

    cy = 0;

    for (i = 0; i < bn; i++)
    {
        sub_ddmmss(hi, lo, 0, a[i], 0, b[i] - cy);
        res[i] = lo + (hi & B);
        cy = hi;
    }

    cy = -cy;

    for (i = bn; i < an; i++)  /* improve */
    {
        res[i] = a[i] - cy;

        if (res[i] == UWORD_MAX)
        {
            cy = 1;
            res[i] = B - 1;
        }
        else
        {
            cy = 0;
            if (res != a)
                flint_mpn_copyi(res + i + 1, a + i + 1, an - i - 1);
            break;
        }
    }

    return cy;
}

