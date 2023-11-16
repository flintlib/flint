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
ca_pi_i(ca_t res, ca_ctx_t ctx)
{
    ca_t t;
    ca_init(t, ctx);
    ca_pi(res, ctx);
    ca_i(t, ctx);
    ca_mul(res, res, t, ctx);
    ca_clear(t, ctx);
}

