/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"

int fmpz_mod_is_invertible(const fmpz_t a, const fmpz_mod_ctx_t ctx)
{
    int res;
    fmpz_t d;
    fmpz_init(d);
    fmpz_gcd(d, a, ctx->n);
    res = fmpz_is_one(d);
    fmpz_clear(d);
    return res;
}
