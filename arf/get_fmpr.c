/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arf.h"

void
arf_get_fmpr(fmpr_t y, const arf_t x)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            fmpr_zero(y);
        else if (arf_is_pos_inf(x))
            fmpr_pos_inf(y);
        else if (arf_is_neg_inf(x))
            fmpr_neg_inf(y);
        else
            fmpr_nan(y);
    }
    else
    {
        arf_get_fmpz_2exp(fmpr_manref(y), fmpr_expref(y), x);
    }
}

