/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"

void fmpz_mod_set_si(fmpz_t a, slong b, const fmpz_mod_ctx_t ctx)
{
    fmpz_set_si(a, b);
    fmpz_mod(a, a, ctx->n);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}
