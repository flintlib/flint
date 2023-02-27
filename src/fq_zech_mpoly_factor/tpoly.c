/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly_factor.h"


void fq_zech_tpoly_fit_length(fq_zech_tpoly_t A, slong len, const fq_zech_ctx_t ctx)
{
    slong i;

    if (len <= A->alloc)
        return;

    if (len < 2 * A->alloc)
        len = 2 * A->alloc;

    if (A->alloc > 0)
        A->coeffs = (fq_zech_bpoly_struct *) flint_realloc(A->coeffs,
                                           len * sizeof(fq_zech_bpoly_struct));
    else
        A->coeffs = (fq_zech_bpoly_struct *) flint_malloc(
                                           len * sizeof(fq_zech_bpoly_struct));

    for (i = A->alloc; i < len; i++)
        fq_zech_bpoly_init(A->coeffs + i, ctx);

    A->alloc = len;
}

void fq_zech_tpoly_clear(fq_zech_tpoly_t A, const fq_zech_ctx_t ctx)
{
    if (A->alloc > 0)
    {
        slong i;
        for (i = 0; i < A->alloc; i++)
            fq_zech_bpoly_clear(A->coeffs + i, ctx);
        flint_free(A->coeffs);
    }
}
