/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"

int fmpz_mod_pow_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t pow,
                                                      const fmpz_mod_ctx_t ctx)
{
    int success = 1;
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));

    if (fmpz_sgn(pow) >= 0)
    {
        fmpz_powm(a, b, pow, ctx->n);
    }
    else
    {
        fmpz_t d;
        fmpz_init(d);
        fmpz_gcdinv(d, a, b, ctx->n);
        if (fmpz_is_one(d))
        {
            fmpz_neg(d, pow);
            fmpz_powm(a, a, d, ctx->n);
        }
        else
        {
            success = 0;
        }
        fmpz_clear(d);
    }

    FLINT_ASSERT(!success || fmpz_mod_is_canonical(a, ctx));
    return success;
}
