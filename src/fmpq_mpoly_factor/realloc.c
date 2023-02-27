/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly_factor.h"


void fmpq_mpoly_factor_realloc(fmpq_mpoly_factor_t f, slong alloc,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    slong i;

    if (alloc <= 0)
    {
        fmpq_mpoly_factor_clear(f, ctx);
        fmpq_mpoly_factor_init(f, ctx);
        return;
    }

    if (f->alloc > 0)
    {
        if (f->alloc > alloc)
        {
            for (i = alloc; i < f->alloc; i++)
            {
                fmpq_mpoly_clear(f->poly + i, ctx);
                fmpz_clear(f->exp + i);
            }

            f->exp  = (fmpz *) flint_realloc(f->exp, alloc * sizeof(fmpz));
            f->poly = (fmpq_mpoly_struct *) flint_realloc(f->poly,
                                            alloc * sizeof(fmpq_mpoly_struct));
        }
        else if (f->alloc < alloc)
        {
            f->exp  = (fmpz *) flint_realloc(f->exp, alloc * sizeof(fmpz));
            f->poly = (fmpq_mpoly_struct *) flint_realloc(f->poly,
                                            alloc * sizeof(fmpq_mpoly_struct));

            for (i = f->alloc; i < alloc; i++)
            {
                fmpq_mpoly_init(f->poly + i, ctx);
                fmpz_init(f->exp + i);
            }
        }
    }
    else
    {
        f->exp  = (fmpz *) flint_calloc(alloc, sizeof(fmpz));
        f->poly = (fmpq_mpoly_struct *) flint_malloc(alloc *
                                                    sizeof(fmpq_mpoly_struct));
        for (i = 0; i < alloc; i++)
            fmpq_mpoly_init(f->poly + i, ctx);
    }

    f->alloc = alloc;
}

