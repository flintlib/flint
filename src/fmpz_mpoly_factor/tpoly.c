/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"
#include "nmod_mpoly_factor.h"


void fmpz_tpoly_fit_length(fmpz_tpoly_t A, slong len)
{
    slong i;

    if (len <= A->alloc)
        return;

    if (len < 2 * A->alloc)
        len = 2 * A->alloc;

    if (A->alloc > 0)
        A->coeffs = (fmpz_bpoly_struct *) flint_realloc(A->coeffs,
                                              len * sizeof(fmpz_bpoly_struct));
    else
        A->coeffs = (fmpz_bpoly_struct *) flint_malloc(
                                              len * sizeof(fmpz_bpoly_struct));

    for (i = A->alloc; i < len; i++)
        fmpz_bpoly_init(A->coeffs + i);

    A->alloc = len;
}

void fmpz_tpoly_clear(fmpz_tpoly_t A)
{
    if (A->alloc > 0)
    {
        slong i;
        for (i = 0; i < A->alloc; i++)
            fmpz_bpoly_clear(A->coeffs + i);
        flint_free(A->coeffs);
    }
}

