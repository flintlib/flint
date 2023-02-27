/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"

void fmpz_mod_sub_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                      const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));

    fmpz_sub(a, b, c);
    fmpz_mod(a, a, ctx->n);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void fmpz_mod_sub_ui(fmpz_t a, const fmpz_t b, ulong c,
                                                      const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));

    fmpz_sub_ui(a, b, c);
    fmpz_mod(a, a, ctx->n);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void fmpz_mod_sub_si(fmpz_t a, const fmpz_t b, slong c,
                                                      const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));

    fmpz_sub_si(a, b, c);
    fmpz_mod(a, a, ctx->n);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void fmpz_mod_fmpz_sub(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                      const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));

    fmpz_sub(a, b, c);
    fmpz_mod(a, a, ctx->n);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void fmpz_mod_ui_sub(fmpz_t a, ulong b, const fmpz_t c,
                                                      const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));

    fmpz_sub_ui(a, c, b);
    fmpz_neg(a, a);
    fmpz_mod(a, a, ctx->n);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

void fmpz_mod_si_sub(fmpz_t a, slong b, const fmpz_t c,
                                                      const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));

    fmpz_sub_si(a, c, b);
    fmpz_neg(a, a);
    fmpz_mod(a, a, ctx->n);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

