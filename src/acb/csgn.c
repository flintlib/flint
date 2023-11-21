/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_csgn(arb_t res, const acb_t z)
{
    if (arb_is_zero(acb_realref(z)))
    {
        arb_sgn(res, acb_imagref(z));
    }
    else
    {
        arb_sgn(res, acb_realref(z));
    }
}

