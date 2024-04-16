/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mpoly.h"

void _fmpz_mpoly_fit_length(fmpz ** poly,
                              ulong ** exps, slong * alloc, slong len, slong N)
{
    if (len > *alloc)
    {
        /* at least double size */
        len = FLINT_MAX(len, 2*(*alloc));
        _fmpz_mpoly_realloc(poly, exps, alloc, len, N);
    }
}

void
fmpz_mpoly_fit_length(fmpz_mpoly_t poly, slong len, const fmpz_mpoly_ctx_t ctx)
{
    if (len > poly->alloc)
    {
        /* At least double number of allocated coeffs */
        if (len < 2 * poly->alloc)
            len = 2 * poly->alloc;
        fmpz_mpoly_realloc(poly, len, ctx);
    }
}

void _fmpz_mpoly_set_length(fmpz_mpoly_t A, slong newlen, const fmpz_mpoly_ctx_t FLINT_UNUSED(ctx))
{
    slong i;
    for (i = newlen; i < A->length; i++)
       _fmpz_demote(A->coeffs + i);

    A->length = newlen;
}

void fmpz_mpoly_truncate(fmpz_mpoly_t A, slong newlen, const fmpz_mpoly_ctx_t FLINT_UNUSED(ctx))
{
    if (A->length > newlen)
    {
        slong i;

        for (i = newlen; i < A->length; i++)
            _fmpz_demote(A->coeffs + i);

        A->length = newlen;
    }
}
