/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_randtest_exact(arb_t x, flint_rand_t state, slong prec, slong mag_bits)
{
    arf_randtest(arb_midref(x), state, prec, mag_bits);
    mag_zero(arb_radref(x));
}

void
arb_randtest_wide(arb_t x, flint_rand_t state, slong prec, slong mag_bits)
{
    arf_randtest(arb_midref(x), state, prec, mag_bits);
    mag_randtest(arb_radref(x), state, mag_bits);
}

void
arb_randtest_precise(arb_t x, flint_rand_t state, slong prec, slong mag_bits)
{
    arf_randtest(arb_midref(x), state, prec, mag_bits);

    if (arf_is_zero(arb_midref(x)) || (n_randint(state, 8) == 0))
    {
        mag_zero(arb_radref(x));
    }
    else
    {
        mag_randtest(arb_radref(x), state, 0);

        if (!mag_is_zero(arb_radref(x)))
        {
            fmpz_add_si(MAG_EXPREF(arb_radref(x)),
                ARF_EXPREF(arb_midref(x)),
                -prec + 2 - n_randint(state, 8));
        }
    }
}

void
arb_randtest_positive(arb_t x, flint_rand_t state, slong prec, slong mag_bits)
{
    int stop = 0;

    while (!stop)
    {
        arb_randtest_precise(x, state, prec, mag_bits);
        stop = !arb_contains_zero(x);
    }
    if (!arb_is_positive(x))
    {
        arb_neg(x, x);
    }
}

void
arb_randtest(arb_t x, flint_rand_t state, slong prec, slong mag_bits)
{
    switch (n_randint(state, 8))
    {
        case 0:
            arb_randtest_exact(x, state, prec, mag_bits);
            break;
        case 1:
            arb_randtest_wide(x, state, prec, mag_bits);
            break;
        default:
            arb_randtest_precise(x, state, prec, mag_bits);
    }
}

void
arb_randtest_special(arb_t x, flint_rand_t state, slong prec, slong mag_bits)
{
    arb_randtest(x, state, prec, mag_bits);

    if (n_randint(state, 10) == 0)
        mag_inf(arb_radref(x));

    switch (n_randint(state, 10))
    {
        case 0:
            arf_pos_inf(arb_midref(x));
            break;
        case 1:
            arf_neg_inf(arb_midref(x));
            break;
        case 2:
            arf_nan(arb_midref(x));
            mag_inf(arb_radref(x));
            break;
        default:
            break;
    }
}

