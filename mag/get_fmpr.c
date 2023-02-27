/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mag.h"

void
mag_get_fmpr(fmpr_t x, const mag_t r)
{
    if (mag_is_zero(r))
    {
        fmpr_zero(x);
    }
    else if (mag_is_inf(r))
    {
        fmpr_pos_inf(x);
    }
    else
    {
        fmpr_set_ui_2exp_si(x, MAG_MAN(r), -MAG_BITS);
        _fmpz_add2_fast(fmpr_expref(x), fmpr_expref(x), MAG_EXPREF(r), 0);
    }
}

