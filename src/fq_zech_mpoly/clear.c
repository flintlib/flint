/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"

void fq_zech_mpoly_clear(fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
{
    slong i;

    if (A->alloc > 0)
    {
        for (i = 0; i < A->alloc; i++)
            fq_zech_clear(A->coeffs + i, ctx->fqctx);

        flint_free(A->coeffs);
        flint_free(A->exps);
    }
}
