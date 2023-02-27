/*
    Copyright (C) 2019 Julian Rüth

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <string.h>

#include "arb.h"
#include "arf.h"
#include "mag.h"

char *
arb_dump_str(const arb_t x)
{
    char * mid;
    char * mag;
    size_t res_len;
    char * res;

    mid = arf_dump_str(arb_midref(x));
    mag = mag_dump_str(arb_radref(x));

    res_len = strlen(mid) + 1 + strlen(mag);
    res = (char*)flint_malloc(res_len + 1);
    strcpy(res, mid);
    strcat(res, " ");
    strcat(res, mag);

    flint_free(mid);
    flint_free(mag);

    return res;
}
