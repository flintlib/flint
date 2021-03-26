/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "mpn_extras.h"
#include "fmpz.h"


static int flint_mpn_cmp2abs(mp_srcptr x, slong xn, mp_srcptr a, slong an)
{
    mp_limb_t xhi, ahi;

    FLINT_ASSERT(an >= 0);
    FLINT_ASSERT(xn >= 0);
    FLINT_ASSERT(xn == 0 || x[xn - 1] != 0);
    FLINT_ASSERT(an == 0 || a[an - 1] != 0);

    if (an > xn)
        return -1;

    if (an + 1 < xn)
        return 1;

    xhi = an == xn ? 0 : x[an];
    ahi = 0;

    while (an > 0)
    {
        ahi = MPN_LEFT_SHIFT_HI(ahi, a[an - 1], 1);
        if (xhi != ahi)
            return xhi > ahi ? 1 : -1;
        ahi = a[an - 1];
        xhi = x[an - 1];
        an--;
    }

    ahi = MPN_LEFT_SHIFT_HI(ahi, UWORD(0), 1);
    if (xhi != ahi)
        return xhi > ahi ? 1 : -1;

    return 0;
}

/* compare |a| and 2|b| */
int fmpz_cmp2abs(const fmpz_t a, const fmpz_t b)
{
    if (!COEFF_IS_MPZ(*b))
    {
        mp_limb_t ub = FLINT_ABS(*b);

        if (!COEFF_IS_MPZ(*a))
        {
            mp_limb_t ua = FLINT_ABS(*a);
            return ua < 2*ub ? -1 : ua > 2*ub ? 1 : 0;
        }
        else
        {
            return flint_mpn_cmp2abs(COEFF_TO_PTR(*a)->_mp_d,
                                     FLINT_ABS(COEFF_TO_PTR(*a)->_mp_size),
                                     &ub, ub != 0);
        }
    }
    else
    {
        if (!COEFF_IS_MPZ(*a))
        {
            return -1;
        }
        else
        {
            return flint_mpn_cmp2abs(COEFF_TO_PTR(*a)->_mp_d,
                                     FLINT_ABS(COEFF_TO_PTR(*a)->_mp_size),
                                     COEFF_TO_PTR(*b)->_mp_d,
                                     FLINT_ABS(COEFF_TO_PTR(*b)->_mp_size));
        }
    }
}

