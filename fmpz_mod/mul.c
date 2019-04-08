/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"

void fmpz_mod_mul(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
{
    fmpz_t t;
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));
    fmpz_init(t);
    fmpz_mul(t, b, c);
    fmpz_mod(a, t, ctx->n);
    fmpz_clear(t);
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}
