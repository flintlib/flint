/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fq_nmod_mpoly_factor.h"

void fq_nmod_mpoly_factor_realloc(fq_nmod_mpoly_factor_t f, slong alloc,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    if (alloc <= 0)
    {
        fq_nmod_mpoly_factor_clear(f, ctx);
        fq_nmod_mpoly_factor_init(f, ctx);
        return;
    }

    if (f->alloc > 0)
    {
        if (f->alloc > alloc)
        {
            for (i = alloc; i < f->alloc; i++)
            {
                fq_nmod_mpoly_clear(f->poly + i, ctx);
                fmpz_clear(f->exp + i);
            }

            f->exp  = (fmpz *) flint_realloc(f->exp, alloc * sizeof(fmpz));
            f->poly = (fq_nmod_mpoly_struct *) flint_realloc(f->poly,
                                         alloc * sizeof(fq_nmod_mpoly_struct));
        }
        else if (f->alloc < alloc)
        {
            f->exp  = (fmpz *) flint_realloc(f->exp, alloc * sizeof(fmpz));
            f->poly = (fq_nmod_mpoly_struct *) flint_realloc(f->poly,
                                         alloc * sizeof(fq_nmod_mpoly_struct));

            for (i = f->alloc; i < alloc; i++)
            {
                fq_nmod_mpoly_init(f->poly + i, ctx);
                fmpz_init(f->exp + i);
            }
        }
    }
    else
    {
        f->exp  = (fmpz *) flint_calloc(alloc, sizeof(fmpz));
        f->poly = (fq_nmod_mpoly_struct *) flint_malloc(alloc *
                                                 sizeof(fq_nmod_mpoly_struct));
        for (i = 0; i < alloc; i++)
            fq_nmod_mpoly_init(f->poly + i, ctx);
    }

    f->alloc = alloc;
}
