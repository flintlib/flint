/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"

void fmpz_mod_addmul(
    fmpz_t a,
    const fmpz_t b,
    const fmpz_t c,
    const fmpz_t d,
    const fmpz_mod_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_mul(t, c, d);
    fmpz_add(t, t, b);
    fmpz_mod(a, t, ctx->n);
    fmpz_clear(t);
}
