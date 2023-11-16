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
mag_add_ui(mag_t y, const mag_t x, ulong k)
{
    mag_t t;
    mag_init(t); /* no need to free */
    mag_set_ui(t, k);
    mag_add(y, x, t);
}

void
mag_add_ui_lower(mag_t res, const mag_t x, ulong k)
{
    mag_t t;
    mag_init(t);
    mag_set_ui_lower(t, k);  /* no need to free */
    mag_add_lower(res, x, t);
}

