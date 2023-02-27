/*
    Copyright (C) 2019 Julian Rüth

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <string.h>

#include "arb.h"

int
arb_load_file(arb_t x, FILE* stream)
{
    int err;

    err = arf_load_file(arb_midref(x), stream);

    if (err) return err;

    err = mag_load_file(arb_radref(x), stream);

    return err;
}

