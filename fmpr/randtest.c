/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

void
fmpr_randtest(fmpr_t x, flint_rand_t state, slong bits, slong mag_bits)
{
    fmpz_randtest(fmpr_manref(x), state, bits);
    fmpz_randtest(fmpr_expref(x), state, mag_bits);
    fmpz_sub_ui(fmpr_expref(x), fmpr_expref(x), fmpz_bits(fmpr_manref(x)));
    _fmpr_normalise(fmpr_manref(x), fmpr_expref(x), bits, FMPR_RND_DOWN);
}

void
fmpr_randtest_not_zero(fmpr_t x, flint_rand_t state, slong bits, slong mag_bits)
{
    fmpz_randtest_not_zero(fmpr_manref(x), state, bits);
    fmpz_randtest(fmpr_expref(x), state, mag_bits);
    fmpz_sub_ui(fmpr_expref(x), fmpr_expref(x), fmpz_bits(fmpr_manref(x)));
    _fmpr_normalise(fmpr_manref(x), fmpr_expref(x), bits, FMPR_RND_DOWN);
}

void
fmpr_randtest_special(fmpr_t x, flint_rand_t state, slong bits, slong mag_bits)
{
    switch (n_randint(state, 32))
    {
        case 0:
            fmpr_zero(x);
            break;
        case 1:
            fmpr_pos_inf(x);
            break;
        case 2:
            fmpr_neg_inf(x);
            break;
        case 3:
            fmpr_nan(x);
            break;
        default:
            fmpr_randtest_not_zero(x, state, bits, mag_bits);
    }
}
