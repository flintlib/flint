/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mag.h"
#include "arf.h"

void
mag_set_fmpr(mag_t x, const fmpr_t y)
{
    if (fmpr_is_special(y))
    {
        if (fmpr_is_zero(y))
            mag_zero(x);
        else
            mag_inf(x);
    }
    else
    {
        mag_set_fmpz_2exp_fmpz(x, fmpr_manref(y), fmpr_expref(y));
    }
}

