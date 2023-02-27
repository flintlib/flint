/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_factor.h"


void fmpz_mod_mpoly_factor_realloc(
    fmpz_mod_mpoly_factor_t f,
    slong alloc,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;

    if (alloc < 1)
    {
        fmpz_mod_mpoly_factor_clear(f, ctx);
        fmpz_mod_mpoly_factor_init(f, ctx);
        return;
    }

    if (f->alloc > 0)
    {
        if (f->alloc > alloc)
        {
            for (i = alloc; i < f->alloc; i++)
            {
                fmpz_mod_mpoly_clear(f->poly + i, ctx);
                fmpz_clear(f->exp + i);
            }

            f->exp  = FLINT_ARRAY_REALLOC(f->exp, alloc, fmpz);
            f->poly = FLINT_ARRAY_REALLOC(f->poly, alloc, fmpz_mod_mpoly_struct);
        }
        else if (f->alloc < alloc)
        {
            f->exp  = FLINT_ARRAY_REALLOC(f->exp, alloc, fmpz);
            f->poly = FLINT_ARRAY_REALLOC(f->poly, alloc, fmpz_mod_mpoly_struct);

            for (i = f->alloc; i < alloc; i++)
            {
                fmpz_mod_mpoly_init(f->poly + i, ctx);
                fmpz_init(f->exp + i);
            }
        }
    }
    else
    {
        f->exp  = (fmpz *) flint_calloc(alloc, sizeof(fmpz));
        f->poly = FLINT_ARRAY_ALLOC(alloc, fmpz_mod_mpoly_struct);
        for (i = 0; i < alloc; i++)
            fmpz_mod_mpoly_init(f->poly + i, ctx);
    }

    f->alloc = alloc;
}
