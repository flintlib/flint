/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_printn(const ca_t x, slong n, ulong flags, ca_ctx_t ctx)
{
    acb_t t;
    acb_init(t);
    ca_get_acb(t, x, n * 3.33 + 64, ctx);
    acb_printn(t, n, flags);
    acb_clear(t);
}
