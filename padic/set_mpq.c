/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic.h"

void padic_set_mpq(padic_t rop, const mpq_t op, const padic_ctx_t ctx)
{
    fmpq_t t;

    fmpq_init(t);
    fmpq_set_mpq(t, op);
    padic_set_fmpq(rop, t, ctx);
    fmpq_clear(t);
}

