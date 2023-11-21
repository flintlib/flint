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
ca_set_fmpq(ca_t x, const fmpq_t v, ca_ctx_t ctx)
{
    _ca_make_fmpq(x, ctx);
    fmpq_set(CA_FMPQ(x), v);
}

