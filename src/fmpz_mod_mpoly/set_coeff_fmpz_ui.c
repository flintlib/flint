/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_mpoly_set_coeff_fmpz_ui(fmpz_mod_mpoly_t poly,
             const fmpz_t c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, nvars = ctx->minfo->nvars;
    fmpz * newexp;
    TMP_INIT;

    TMP_START;
    newexp = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (i = 0; i < nvars; i++)
        fmpz_init_set_ui(newexp + i, exp[i]);

    _fmpz_mod_mpoly_set_coeff_fmpz_fmpz(poly, c, newexp, ctx);

    for (i = 0; i < nvars; i++)
        fmpz_clear(newexp + i);

    TMP_END;
}

void fmpz_mod_mpoly_set_coeff_ui_ui(fmpz_mod_mpoly_t poly,
                    ulong c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t C;
    fmpz_init_set_ui(C, c);
    fmpz_mod_mpoly_set_coeff_fmpz_ui(poly, C, exp, ctx);
    fmpz_clear(C);
}

void fmpz_mod_mpoly_set_coeff_si_ui(fmpz_mod_mpoly_t poly,
                    slong c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t C;
    fmpz_init_set_si(C, c);
    fmpz_mod_mpoly_set_coeff_fmpz_ui(poly, C, exp, ctx);
    fmpz_clear(C);
}

