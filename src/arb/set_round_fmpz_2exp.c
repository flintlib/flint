/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_set_round_fmpz_2exp(arb_t y, const fmpz_t x, const fmpz_t exp, slong prec)
{
    int inexact;
    inexact = arf_set_round_fmpz_2exp(arb_midref(y), x, exp, prec, ARB_RND);

    if (inexact)
        arf_mag_set_ulp(arb_radref(y), arb_midref(y), prec);
    else
        mag_zero(arb_radref(y));
}

