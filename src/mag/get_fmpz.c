/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mag.h"
#include "arf.h"

void
mag_get_fmpz(fmpz_t res, const mag_t x)
{
    arf_t t;
    arf_init_set_mag_shallow(t, x);
    arf_get_fmpz(res, t, ARF_RND_UP);
}

void
mag_get_fmpz_lower(fmpz_t res, const mag_t x)
{
    arf_t t;
    arf_init_set_mag_shallow(t, x);
    arf_get_fmpz(res, t, ARF_RND_DOWN);
}

