/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"


void fmpz_mpoly_factor_init2(fmpz_mpoly_factor_t f, slong alloc,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    fmpz_init_set_ui(f->constant, 1);

    if (alloc > 0)
    {
        slong i;

        f->exp = (fmpz *) flint_calloc(alloc, sizeof(fmpz));
        f->poly = (fmpz_mpoly_struct *) flint_malloc(alloc *
                                                    sizeof(fmpz_mpoly_struct));
        for (i = 0; i < alloc; i++)
            fmpz_mpoly_init(f->poly + i, ctx);

        f->alloc = alloc;
    }
    else
    {
        f->exp = NULL;
        f->poly = NULL;
        f->alloc = 0;
    }

    f->num = 0;
}

