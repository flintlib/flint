/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "calcium.h"

void
calcium_write_arb(calcium_stream_t out, const arb_t x, slong digits, ulong flags)
{
    calcium_write_free(out, arb_get_str(x, digits, flags));
}

void
calcium_write_acb(calcium_stream_t out, const acb_t z, slong digits, ulong flags)
{
    if (arb_is_zero(acb_imagref(z)))
    {
        calcium_write_arb(out, acb_realref(z), digits, flags);
    }
    else if (arb_is_zero(acb_realref(z)))
    {
        calcium_write_arb(out, acb_imagref(z), digits, flags);
        calcium_write(out, "*I");
    }
    else
    {
        calcium_write_arb(out, acb_realref(z), digits, flags);

        if ((arb_is_exact(acb_imagref(z)) || (flags & ARB_STR_NO_RADIUS))
                && arf_sgn(arb_midref(acb_imagref(z))) < 0)
        {
            arb_t t;
            arb_init(t);
            arb_neg(t, acb_imagref(z));
            calcium_write(out, " - ");
            calcium_write_arb(out, t, digits, flags);
            arb_clear(t);
        }
        else
        {
            calcium_write(out, " + ");
            calcium_write_arb(out, acb_imagref(z), digits, flags);
        }

        calcium_write(out, "*I");
    }
}
