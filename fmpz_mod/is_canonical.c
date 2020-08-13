/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"

int fmpz_mod_is_canonical(const fmpz_t a, const fmpz_mod_ctx_t ctx)
{
    return fmpz_sgn(a) >= 0 && fmpz_cmp(a, ctx->n) < 0;
}

void fmpz_mod_assert_canonical(const fmpz_t a, const fmpz_mod_ctx_t ctx)
{
    if (!fmpz_mod_is_canonical(a, ctx))
        flint_throw(FLINT_ERROR, "Fmpz mod invalid");
}
