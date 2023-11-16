/*
    Copyright (C) 2019 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

int
arb_sgn_nonzero(const arb_t x)
{
    if (arb_is_positive(x))
        return 1;
    else if (arb_is_negative(x))
        return -1;
    else
        return 0;
}
