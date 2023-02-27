/*
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2015 Arb authors

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_fprintd(FILE * file, const acb_t z, slong digits)
{
    flint_fprintf(file, "(");
    arf_fprintd(file, arb_midref(acb_realref(z)), digits);

    if (arf_sgn(arb_midref(acb_imagref(z))) < 0)
    {
        arf_t t;
        arf_init_neg_shallow(t, arb_midref(acb_imagref(z)));
        flint_fprintf(file, " - ");
        arf_fprintd(file, t, digits);
    }
    else
    {
        flint_fprintf(file, " + ");
        arf_fprintd(file, arb_midref(acb_imagref(z)), digits);
    }
    flint_fprintf(file, "j)");

    flint_fprintf(file, "  +/-  ");

    flint_fprintf(file, "(");
    mag_fprintd(file, arb_radref(acb_realref(z)), 3);
    flint_fprintf(file, ", ");
    mag_fprintd(file, arb_radref(acb_imagref(z)), 3);
    flint_fprintf(file, "j)");
}
