/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "calcium.h"

void
calcium_write_acb(calcium_stream_t out, const acb_t z, slong digits, ulong flags)
{
    if (arb_is_zero(acb_imagref(z)))
    {
        CALCIUM_WRITE_FREE(out, arb_get_str(acb_realref(z), digits, flags));
    }
    else if (arb_is_zero(acb_realref(z)))
    {
        CALCIUM_WRITE_FREE(out, arb_get_str(acb_imagref(z), digits, flags));
        calcium_write(out, "*I");
    }
    else
    {
        CALCIUM_WRITE_FREE(out, arb_get_str(acb_realref(z), digits, flags));

        if ((arb_is_exact(acb_imagref(z)) || (flags & ARB_STR_NO_RADIUS))
                && arf_sgn(arb_midref(acb_imagref(z))) < 0)
        {
            arb_t t;
            arb_init(t);
            arb_neg(t, acb_imagref(z));
            calcium_write(out, " - ");
            CALCIUM_WRITE_FREE(out, arb_get_str(t, digits, flags));
            arb_clear(t);
        }
        else
        {
            calcium_write(out, " + ");
            CALCIUM_WRITE_FREE(out, arb_get_str(acb_imagref(z), digits, flags));
        }

        calcium_write(out, "*I");
    }
}
