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
arf_mul_special(arf_t z, const arf_t x, const arf_t y)
{
    if (arf_is_zero(x))
    {
        if (arf_is_finite(y))
            arf_zero(z);
        else
            arf_nan(z);
    }
    else if (arf_is_zero(y))
    {
        if (arf_is_finite(x))
            arf_zero(z);
        else
            arf_nan(z);
    }
    else if (arf_is_nan(x) || arf_is_nan(y))
        arf_nan(z);
    else if (arf_sgn(x) == arf_sgn(y))
        arf_pos_inf(z);
    else
        arf_neg_inf(z);
}

