/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"

void fmpz_mod_add_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                      const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));

    fmpz_add(a, b, c);
    fmpz_mod(a, a, ctx->n);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void fmpz_mod_add_ui(fmpz_t a, const fmpz_t b, ulong c,
                                                      const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));

    fmpz_add_ui(a, b, c);
    fmpz_mod(a, a, ctx->n);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void fmpz_mod_add_si(fmpz_t a, const fmpz_t b, slong c,
                                                      const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));

    fmpz_add_si(a, b, c);
    fmpz_mod(a, a, ctx->n);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}
