/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"

void fmpz_mod_inv(fmpz_t a, const fmpz_t b, const fmpz_mod_ctx_t ctx)
{
    fmpz_t d;
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    fmpz_init(d);
    fmpz_gcdinv(d, a, b, ctx->n);
    if (!fmpz_is_one(d))
    {
        fmpz_clear(d);
        /* printing b and n would entail leaking the string so no b nor n */
        flint_throw(FLINT_IMPINV, "Exception in fmpz_mod_inv: Cannot invert.\n");
    }
    fmpz_clear(d);
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}
