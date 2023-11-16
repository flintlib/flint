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
arf_randtest(arf_t x, flint_rand_t state, slong bits, slong mag_bits)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_randtest(t, state, bits);
    arf_set_fmpz(x, t);
    if (!arf_is_zero(x))
        fmpz_randtest(ARF_EXPREF(x), state, mag_bits);
    fmpz_clear(t);
}

void
arf_randtest_not_zero(arf_t x, flint_rand_t state, slong bits, slong mag_bits)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_randtest_not_zero(t, state, bits);
    arf_set_fmpz(x, t);
    fmpz_randtest(ARF_EXPREF(x), state, mag_bits);
    fmpz_clear(t);
}

void
arf_randtest_special(arf_t x, flint_rand_t state, slong bits, slong mag_bits)
{
    switch (n_randint(state, 32))
    {
        case 0:
            arf_zero(x);
            break;
        case 1:
            arf_pos_inf(x);
            break;
        case 2:
            arf_neg_inf(x);
            break;
        case 3:
            arf_nan(x);
            break;
        default:
            arf_randtest_not_zero(x, state, bits, mag_bits);
    }
}

