/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_poly.h"
#include "fq_zech_mpoly_factor.h"

void fq_zech_polyun_clear(fq_zech_polyun_t A, const fq_zech_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fq_zech_poly_clear(A->coeffs + i, ctx);
    flint_free(A->coeffs);
    flint_free(A->exps);
}

void fq_zech_polyun_realloc(fq_zech_polyun_t A, slong len, const fq_zech_ctx_t ctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(len, old_alloc + 1 + old_alloc/2);

    FLINT_ASSERT(A->alloc >= 0);
    if (len <= A->alloc)
        return;

    A->coeffs = FLINT_ARRAY_REALLOC(A->coeffs, new_alloc, fq_zech_poly_struct);
    A->exps = FLINT_ARRAY_REALLOC(A->exps, new_alloc, ulong);

    for (i = old_alloc; i < new_alloc; i++)
        fq_zech_poly_init(A->coeffs + i, ctx);

    A->alloc = new_alloc;
}

int fq_zech_polyun_is_canonical(
    const fq_zech_polyun_t A,
    const fq_zech_ctx_t ctx)
{
    slong i;
    if (A->length < 0)
        return 0;
    for (i = 0; i < A->length; i++)
    {
        if (fq_zech_poly_is_zero(A->coeffs + i, ctx))
            return 0;
        if (i > 0 && A->exps[i] >= A->exps[i - 1])
            return 0;
    }
    return 1;
}
