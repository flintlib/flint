/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "acb_theta.h"

slong
acb_theta_jet_nb(slong ord, slong g)
{
    fmpz_t x;
    slong res;

    FLINT_ASSERT(g >= 0);

    if (ord < 0)
    {
        return 0;
    }

    fmpz_init(x);
    fmpz_bin_uiui(x, g + ord, g);
    res = fmpz_get_si(x);

    fmpz_clear(x);
    return res;
}
