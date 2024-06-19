/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod.h"

int fmpz_mod_is_one(const fmpz_t a, const fmpz_mod_ctx_t ctx)
{
    return fmpz_is_one(a) || fmpz_is_one(ctx->n);
}
