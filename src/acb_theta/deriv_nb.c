/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

slong acb_theta_deriv_nb(slong ord, slong g)
{
    fmpz_t x;
    slong res;

    fmpz_init(x);
    fmpz_bin_uiui(x, g + ord - 1, g - 1);
    res = fmpz_get_si(x);

    fmpz_clear(x);
    return res;
}
