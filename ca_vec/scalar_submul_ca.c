/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_vec.h"

void
_ca_vec_scalar_submul_ca(ca_ptr res, ca_srcptr vec, slong len, const ca_t c, ca_ctx_t ctx)
{
    slong i;
    ca_t t;

    if (len > 0)
    {
        ca_init(t, ctx);
        for (i = 0; i < len; i++)
        {
            ca_mul(t, vec + i, c, ctx);
            ca_sub(res + i, res + i, t, ctx);
        }
        ca_clear(t, ctx);
    }
}
