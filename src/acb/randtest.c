/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_randtest(acb_t z, flint_rand_t state, slong prec, slong mag_bits)
{
    arb_randtest(acb_realref(z), state, prec, mag_bits);
    arb_randtest(acb_imagref(z), state, prec, mag_bits);
}

void
acb_randtest_special(acb_t z, flint_rand_t state, slong prec, slong mag_bits)
{
    arb_randtest_special(acb_realref(z), state, prec, mag_bits);
    arb_randtest_special(acb_imagref(z), state, prec, mag_bits);
}

void
acb_randtest_precise(acb_t z, flint_rand_t state, slong prec, slong mag_bits)
{
    arb_randtest_precise(acb_realref(z), state, prec, mag_bits);
    arb_randtest_precise(acb_imagref(z), state, prec, mag_bits);
}

void
acb_randtest_param(acb_t x, flint_rand_t state, slong prec, slong size)
{
    if (n_randint(state, 8) == 0)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_randtest(t, state, 1 + n_randint(state, prec));
        arb_set_fmpz(acb_realref(x), t);
        arb_zero(acb_imagref(x));
        acb_mul_2exp_si(x, x, -1);
        fmpz_clear(t);
    }
    else
    {
        acb_randtest(x, state, prec, size);
    }
}

