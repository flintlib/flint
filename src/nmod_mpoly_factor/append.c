/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"

void nmod_mpoly_factor_append_ui(
    nmod_mpoly_factor_t f,
    const nmod_mpoly_t A,
    ulong e,
    const nmod_mpoly_ctx_t ctx)
{
    slong i = f->num;
    nmod_mpoly_factor_fit_length(f, i + 1, ctx);
    nmod_mpoly_set(f->poly + i, A, ctx);
    fmpz_set_ui(f->exp + i, e);
    f->num = i + 1;
}

void nmod_mpoly_factor_append_fmpz(
    nmod_mpoly_factor_t f,
    const nmod_mpoly_t A,
    const fmpz_t e,
    const nmod_mpoly_ctx_t ctx)
{
    slong i = f->num;
    nmod_mpoly_factor_fit_length(f, i + 1, ctx);
    nmod_mpoly_set(f->poly + i, A, ctx);
    fmpz_set(f->exp + i, e);
    f->num = i + 1;
}

