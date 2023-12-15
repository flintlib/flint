/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fmpz.h"
#include "mpn_extras.h"


mp_size_t flint_mpn_remove_power_ascending(mp_ptr x, mp_size_t xsize,
                                      mp_ptr p, mp_size_t psize, ulong *exp)
{
    int i, maxi;
    mp_ptr div;
    mp_ptr rem;
    mp_ptr square[FLINT_BITS];
    mp_size_t square_size[FLINT_BITS];
    mp_size_t sqsize;

    *exp = 0;

    if (psize > xsize)
        return xsize;

    maxi = 0;
    square[0] = p;
    square_size[0] = psize;

    /* Most likely less memory will be needed, but this way we
       avoid reallocations */
    div = flint_malloc(sizeof(mp_limb_t) * xsize);
    rem = flint_malloc(sizeof(mp_limb_t) * xsize);

    /* Remove ascending powers */
    for (i = 0; i < FLINT_BITS && xsize >= square_size[i]; i++)
    {
        mpn_tdiv_qr(div, rem, 0, x, xsize, square[i], square_size[i]);
        if (!flint_mpn_zero_p(rem, square_size[i]))
        {
            i -= 1;
            break;
        }

        *exp += (1 << i);
        xsize = xsize - square_size[i] + 1;
        if (div[xsize-1] == 0)
            xsize--;
        flint_mpn_copyi(x, div, xsize);

        /* Form next square if needed */
        sqsize = square_size[i] * 2;
        if (sqsize - 1 > xsize)
            break;
        maxi = i + 1;
        square[i + 1] = flint_malloc(sizeof(mp_limb_t) * sqsize);
        flint_mpn_sqr(square[i + 1], square[i], square_size[i]);
        if (square[i + 1][sqsize - 1] == 0)
            sqsize -= 1;
        square_size[i + 1] = sqsize;
   }

    /* Remove descending powers */
    for ( ; i >= 0; i--)
    {
        if (xsize >= square_size[i])
        {
            mpn_tdiv_qr(div, rem, 0, x, xsize, square[i], square_size[i]);
            if (flint_mpn_zero_p(rem, square_size[i]))
            {
                *exp += (1 << i);
                xsize = xsize - square_size[i] + 1;
                if (div[xsize-1] == 0)
                    xsize--;
                flint_mpn_copyi(x, div, xsize);
            }
        }
    }

    for (i = 1; i <= maxi; i++)
        flint_free(square[i]);
    flint_free(div);
    flint_free(rem);
    return xsize;
}
