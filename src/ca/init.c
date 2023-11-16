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
ca_init(ca_t x, ca_ctx_t ctx)
{
    x->field = (ulong) ctx->field_qq;
    *CA_FMPQ_NUMREF(x) = 0;
    *CA_FMPQ_DENREF(x) = 1;
}

