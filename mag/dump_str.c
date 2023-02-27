/*
    Copyright (C) 2019 Julian Rüth

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arf.h"
#include "mag.h"

char *
mag_dump_str(const mag_t x)
{
    char * res;
    arf_t y;

    arf_init_set_mag_shallow(y, x);

    res = arf_dump_str(y);
    return res;
}
