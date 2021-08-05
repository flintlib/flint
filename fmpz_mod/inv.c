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
        if (fmpz_is_zero(b))
        {
            fmpz_clear(d);
            flint_exception(FLINT_DIVZERO, "fmpz_mod_inv: Division by zero.",
                            FLINT_EXC_EXTRA_NONE, NULL);
        }
        else
        {
            fmpz * ex = FLINT_ARRAY_ALLOC(1, fmpz);
            *ex = *d;
            flint_exception(FLINT_IMPINV, "fmpz_mod_inv: Modulus has divisor ",
                            FLINT_EXC_EXTRA_FMPZ, ex);
        }
    }
    fmpz_clear(d);
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}
