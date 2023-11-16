/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

void
arf_set_mpn(arf_t y, mp_srcptr x, mp_size_t xn, int sgnbit)
{
    unsigned int leading;
    mp_size_t yn, xn1;
    mp_ptr yptr;
    mp_limb_t bot;

    xn1 = xn;

    while (x[0] == 0)
    {
        x++;
        xn--;
    }

    leading = flint_clz(x[xn - 1]);

    bot = x[0];

    /* This works when leading == 0, since x[0] != 0. */
    yn = xn - ((bot << leading) == 0);

    ARF_GET_MPN_WRITE(yptr, yn, y);
    ARF_XSIZE(y) |= sgnbit;

    if (leading == 0)
    {
        flint_mpn_copyi(yptr, x, xn);
    }
    else if (xn == yn)
    {
        mpn_lshift(yptr, x, yn, leading);
    }
    else
    {
        mpn_lshift(yptr, x + 1, yn, leading);
        yptr[0] |= (bot >> (FLINT_BITS - leading));
    }

    fmpz_set_ui(ARF_EXPREF(y), xn1 * FLINT_BITS - leading);
}

