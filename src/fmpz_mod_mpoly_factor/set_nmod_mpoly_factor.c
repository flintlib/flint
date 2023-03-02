/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_factor.h"

void _fmpz_mod_mpoly_factor_set_nmod_mpoly_factor(
    fmpz_mod_mpoly_factor_t f,
    const fmpz_mod_mpoly_ctx_t ctx,
    const nmod_mpoly_factor_t nf,
    const nmod_mpoly_ctx_t nctx)
{
    slong i;
    fmpz_set_ui(f->constant, nf->constant);
    fmpz_mod_mpoly_factor_fit_length(f, nf->num, ctx);
    f->num = nf->num;
    for (i = 0; i < nf->num; i++)
    {
        fmpz_set(f->exp + i, nf->exp + i);
        _fmpz_mod_mpoly_set_nmod_mpoly(f->poly + i, ctx, nf->poly + i, nctx);
    }
}
