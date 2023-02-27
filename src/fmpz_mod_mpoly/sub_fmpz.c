/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_mpoly_sub_fmpz(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_t c,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t cc;
    fmpz_init(cc);
    fmpz_mod_set_fmpz(cc, c, ctx->ffinfo);
    fmpz_mod_neg(cc, cc, ctx->ffinfo);
    fmpz_mod_mpoly_add_fmpz_mod(A, B, cc, ctx);
    fmpz_clear(cc);
}

void fmpz_mod_mpoly_sub_ui(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    ulong c,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t cc;
    fmpz_init(cc);
    fmpz_mod_set_ui(cc, c, ctx->ffinfo);
    fmpz_mod_neg(cc, cc, ctx->ffinfo);
    fmpz_mod_mpoly_add_fmpz_mod(A, B, cc, ctx);
    fmpz_clear(cc);
}

void fmpz_mod_mpoly_sub_si(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    slong c,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t cc;
    fmpz_init(cc);
    fmpz_mod_set_si(cc, c, ctx->ffinfo);
    fmpz_mod_neg(cc, cc, ctx->ffinfo);
    fmpz_mod_mpoly_add_fmpz_mod(A, B, cc, ctx);
    fmpz_clear(cc);
}
