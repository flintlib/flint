/*
    Copyright (C) 2015 Fredrik Johansson
    Copyright (C) 2015 Arb authors

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_fprint(FILE * file, const arb_t x)
{
    arf_fprint(file, arb_midref(x));
    flint_fprintf(file, " +/- ");
    mag_fprint(file, arb_radref(x));
}

void
arb_fprintd(FILE * file, const arb_t x, slong digits)
{
    arf_fprintd(file, arb_midref(x), FLINT_MAX(digits, 1));
    flint_fprintf(file, " +/- ");
    mag_fprintd(file, arb_radref(x), 5);
}

void
arb_fprintn(FILE * file, const arb_t x, slong digits, ulong flags)
{
    char * s = arb_get_str(x, digits, flags);
    flint_fprintf(file, "%s", s);
    flint_free(s);
}

