/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"

void
mag_randtest_special(mag_t x, flint_rand_t state, slong expbits)
{
    switch (n_randint(state, 32))
    {
        case 0:
            mag_zero(x);
            break;
        case 1:
            mag_inf(x);
            break;
        case 2:
            MAG_MAN(x) = (LIMB_ONE << MAG_BITS) - 1;
            fmpz_randtest(MAG_EXPREF(x), state, expbits);
            break;
        case 3:
            MAG_MAN(x) = LIMB_ONE << (MAG_BITS - 1);
            fmpz_randtest(MAG_EXPREF(x), state, expbits);
            break;
        default:
            MAG_MAN(x) = (n_randtest(state) >> (FLINT_BITS - MAG_BITS));
            MAG_MAN(x) |= (LIMB_ONE << (MAG_BITS - 1));
            fmpz_randtest(MAG_EXPREF(x), state, expbits);
    }
}

void
mag_randtest(mag_t x, flint_rand_t state, slong expbits)
{
    mag_randtest_special(x, state, expbits);
    if (mag_is_inf(x))
        mag_zero(x);
}

