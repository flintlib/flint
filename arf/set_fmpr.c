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
arf_set_fmpr(arf_t y, const fmpr_t x)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
            arf_zero(y);
        else if (fmpr_is_pos_inf(x))
            arf_pos_inf(y);
        else if (fmpr_is_neg_inf(x))
            arf_neg_inf(y);
        else
            arf_nan(y);
    }
    else
    {
        arf_set_fmpz(y, fmpr_manref(x));
        fmpz_add_inline(ARF_EXPREF(y), ARF_EXPREF(y), fmpr_expref(x));
    }
}

