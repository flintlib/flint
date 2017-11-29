/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void _nmod_mpoly_fit_length(ulong ** poly,
                              ulong ** exps, slong * alloc, slong len, slong N, const nmodf_ctx_t fctx)
{
    if (len > *alloc)
    {
        /* at least double size */
        len = FLINT_MAX(len, 2*(*alloc));
        _nmod_mpoly_realloc(poly, exps, alloc, len, N, fctx);
    }
}

void
nmod_mpoly_fit_length(nmod_mpoly_t poly, slong len, const nmod_mpoly_ctx_t ctx)
{
    if (len > poly->alloc)
    {
        /* At least double number of allocated coeffs */
        if (len < 2 * poly->alloc)
            len = 2 * poly->alloc;
        nmod_mpoly_realloc(poly, len, ctx);
    }
}
