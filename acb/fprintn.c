/*
    Copyright (C) 2015 Fredrik Johansson
    Copyright (C) 2015 Arb authors

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_fprintn(FILE * file, const acb_t z, slong digits, ulong flags)
{
    if (arb_is_zero(acb_imagref(z)))
    {
        arb_fprintn(file, acb_realref(z), digits, flags);
    }
    else if (arb_is_zero(acb_realref(z)))
    {
        arb_fprintn(file, acb_imagref(z), digits, flags);
        flint_fprintf(file, "*I");
    }
    else
    {
        arb_fprintn(file, acb_realref(z), digits, flags);

        if ((arb_is_exact(acb_imagref(z)) || (flags & ARB_STR_NO_RADIUS))
                && arf_sgn(arb_midref(acb_imagref(z))) < 0)
        {
            arb_t t;
            arb_init(t);
            arb_neg(t, acb_imagref(z));
            flint_fprintf(file, " - ");
            arb_fprintn(file, t, digits, flags);
            arb_clear(t);
        }
        else
        {
            flint_fprintf(file, " + ");
            arb_fprintn(file, acb_imagref(z), digits, flags);
        }

        flint_fprintf(file, "*I");
    }
}

