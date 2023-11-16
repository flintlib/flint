/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"

void
mag_add_ui_2exp_si(mag_t z, const mag_t x, ulong y, slong e)
{
    mag_t t;
    mag_init(t);
    mag_set_ui_2exp_si(t, y, e);
    mag_add(z, x, t);
    mag_clear(t);
}

