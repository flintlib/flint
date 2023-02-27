/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"

void nmod_mpoly_factor_init2(nmod_mpoly_factor_t f, slong alloc,
                                                    const nmod_mpoly_ctx_t ctx)
{
	f->constant = 1;

    if (alloc > 0)
    {
        slong i;

        f->exp  = (fmpz *) flint_malloc(alloc * sizeof(fmpz));
        f->poly = (nmod_mpoly_struct *) flint_malloc(alloc *
                                                    sizeof(nmod_mpoly_struct));
        for (i = 0; i < alloc; i++)
        {
            nmod_mpoly_init(f->poly + i, ctx);
			fmpz_init(f->exp + i);
        }

        f->alloc  = alloc;
    }
    else
    {
        f->poly = NULL;
        f->exp  = NULL;
        f->alloc = 0;
    }

    f->num = 0;
}

