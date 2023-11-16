/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_set_d(ca_t res, double x, ca_ctx_t ctx)
{
    arf_t t;
    arf_init(t);
    arf_set_d(t, x);

    if (arf_is_finite(t))
    {
        _ca_make_fmpq(res, ctx);
        arf_get_fmpq(CA_FMPQ(res), t);
    }
    else
    {
        if (arf_is_pos_inf(t))
            ca_pos_inf(res, ctx);
        else if (arf_is_neg_inf(t))
            ca_neg_inf(res, ctx);
        else
            ca_unknown(res, ctx);   /* or undefined? */
    }

    /* no need to free t */
}

