/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"
#include "arf.h"

void
mag_get_fmpq(fmpq_t y, const mag_t x)
{
    arf_t t;
    arf_init_set_mag_shallow(t, x);
    arf_get_fmpq(y, t);
}

