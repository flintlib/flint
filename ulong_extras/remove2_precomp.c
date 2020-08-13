/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int
n_remove2_precomp(mp_limb_t * n, mp_limb_t p, double ppre)
{
    int exp = 0;
    mp_limb_t quot, rem = UWORD(0);

    if (p == 2)
    {
        count_trailing_zeros(exp, *n);
        if (exp)
            (*n) >>= exp;

        return exp;
    }

    do
    {
        if ((*n) < p)
            break;
        rem = n_divrem2_precomp(&quot, *n, p, ppre);
        if (rem)
            break;
        exp++;
        (*n) = quot;
    } while (1);

    return exp;
}
