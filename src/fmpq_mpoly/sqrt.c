/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


int fmpq_mpoly_sqrt(fmpq_mpoly_t Q, const fmpq_mpoly_t A,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpz_t r;

    fmpz_init(r);

    if (fmpz_sgn(fmpq_numref(A->content)) < 0)
        goto fail;

    fmpz_sqrtrem(fmpq_numref(Q->content), r, fmpq_numref(A->content));
    if (!fmpz_is_zero(r))
        goto fail;

    fmpz_sqrtrem(fmpq_denref(Q->content), r, fmpq_denref(A->content));
    if (!fmpz_is_zero(r))
        goto fail;

    if (!fmpz_mpoly_sqrt(Q->zpoly, A->zpoly, ctx->zctx))
        goto fail;

    fmpz_clear(r);
    return 1;

fail:

    fmpq_mpoly_zero(Q, ctx);
    fmpz_clear(r);
    return 0;
}

