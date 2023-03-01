/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fexpr.h"

void
fexpr_set_ui(fexpr_t res, ulong c)
{
    if (c <= FEXPR_COEFF_MAX)
    {
        res->data[0] = (c << FEXPR_TYPE_BITS);
    }
    else
    {
        fexpr_fit_size(res, 2);
        res->data[0] = FEXPR_TYPE_BIG_INT_POS | (2 << FEXPR_TYPE_BITS);
        res->data[1] = c;
    }
}
