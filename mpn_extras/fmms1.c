/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "mpn_extras.h"

/*
    Try to compute y = a1*x1 - a2*x2 where x1 and x2 both have n > 0 limbs.
    If y is negative or does not fit into n limbs, the return is -1. Otherwise,
    the return is the size of y.
*/
mp_size_t flint_mpn_fmms1(mp_ptr y, mp_limb_t a1, mp_srcptr x1,
                                    mp_limb_t a2, mp_srcptr x2, mp_size_t n)
{
    mp_limb_t h0, h1;

    FLINT_ASSERT(n > 0);
    h0 =    mpn_mul_1(y, x1, n, a1);
    h1 = mpn_submul_1(y, x2, n, a2);

    if (h1 != h0)
        return -1;

    while (n > 0 && y[n - 1] == 0)
        n--;

    return n;
}
