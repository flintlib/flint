/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_get_abs_lbound_arf(arf_t u, const acb_t z, slong prec)
{
    if (arb_is_zero(acb_imagref(z)))
    {
        arb_get_abs_lbound_arf(u, acb_realref(z), prec);
    }
    else if (arb_is_zero(acb_realref(z)))
    {
        arb_get_abs_lbound_arf(u, acb_imagref(z), prec);
    }
    else
    {
        arf_t v;
        arf_init(v);

        arb_get_abs_lbound_arf(u, acb_realref(z), prec);
        arb_get_abs_lbound_arf(v, acb_imagref(z), prec);

        arf_mul(u, u, u, prec, ARF_RND_DOWN);
        arf_mul(v, v, v, prec, ARF_RND_DOWN);
        arf_add(u, u, v, prec, ARF_RND_DOWN);
        arf_sqrt(u, u, prec, ARF_RND_DOWN);

        arf_clear(v);
    }
}
