/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_get_rad_ubound_arf(arf_t u, const acb_t z, slong prec)
{
    /* fixme: this bound is very sloppy */

    if (mag_cmp(arb_radref(acb_realref(z)), arb_radref(acb_imagref(z))) >= 0)
        arf_set_mag(u, arb_radref(acb_realref(z)));
    else
        arf_set_mag(u, arb_radref(acb_imagref(z)));

    arf_mul_2exp_si(u, u, 1);
}
