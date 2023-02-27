/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"

void fmpz_mod_rand(fmpz_t a, flint_rand_t state, const fmpz_mod_ctx_t ctx)
{
    fmpz_randm(a, state, fmpz_mod_ctx_modulus(ctx));
}

void fmpz_mod_rand_not_zero(fmpz_t a, flint_rand_t state, const fmpz_mod_ctx_t ctx)
{
    slong i = 3;
    for (i = 0; i < 3; i++)
    {
        fmpz_randm(a, state, fmpz_mod_ctx_modulus(ctx));
        if (!fmpz_is_zero(a))
            return;
    }
    fmpz_one(a);
}
