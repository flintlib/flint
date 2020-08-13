/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic.h"

void padic_get_mpz(mpz_t rop, const padic_t op, const padic_ctx_t ctx)
{
    fmpz_t t;

    fmpz_init(t);
    padic_get_fmpz(t, op, ctx);
    fmpz_get_mpz(rop, t);
    fmpz_clear(t);
}

