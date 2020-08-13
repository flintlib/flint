/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"

void _fmpz_mpoly_realloc(fmpz ** poly, ulong ** exps,
                                             slong * alloc, slong len, slong N)
{
    (*poly) = (fmpz *) flint_realloc(*poly, len*sizeof(fmpz));
    (*exps) = (ulong *) flint_realloc(*exps, len*N*sizeof(ulong));

    if (len > *alloc)
        memset(*poly + *alloc, 0, (len - *alloc)*sizeof(fmpz));
    
    (*alloc) = len;
}

void fmpz_mpoly_realloc(fmpz_mpoly_t poly,
                                       slong alloc, const fmpz_mpoly_ctx_t ctx)
{
    slong N;

    if (alloc == 0)             /* Clear up, reinitialise */
    {
        fmpz_mpoly_clear(poly, ctx);
        fmpz_mpoly_init(poly, ctx);

        return;
    }

    N = mpoly_words_per_exp(poly->bits, ctx->minfo);

    if (poly->alloc != 0)            /* Realloc */
    {
        fmpz_mpoly_truncate(poly, alloc, ctx);

        poly->coeffs = (fmpz *) flint_realloc(poly->coeffs, alloc*sizeof(fmpz));
        poly->exps = (ulong *) flint_realloc(poly->exps, alloc*N*sizeof(ulong));

        if (alloc > poly->alloc)
            memset(poly->coeffs + poly->alloc, 0,
                                           (alloc - poly->alloc)*sizeof(fmpz));
    }
    else                        /* Nothing allocated already so do it now */
    {
        poly->coeffs = (fmpz *) flint_calloc(alloc, sizeof(fmpz));
        poly->exps   = (ulong *) flint_malloc(alloc*N*sizeof(ulong));
    }

    poly->alloc = alloc;
}

