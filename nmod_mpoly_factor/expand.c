/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


int nmod_mpoly_factor_expand(
    nmod_mpoly_t A,
    const nmod_mpoly_factor_t f,
    const nmod_mpoly_ctx_t ctx)
{
    int success = 1;
    slong i;
    nmod_mpoly_t t1, t2;

    nmod_mpoly_init(t1, ctx);
    nmod_mpoly_init(t2, ctx);

    nmod_mpoly_set_ui(A, f->constant, ctx);

    for (i = 0; i < f->num; i++)
    {
        if (fmpz_sgn(f->exp + i) < 0 ||
            !nmod_mpoly_pow_fmpz(t1, f->poly + i, f->exp + i, ctx))
        {
            success = 0;
            goto cleanup;
        }
        nmod_mpoly_mul(t2, A, t1, ctx);
        nmod_mpoly_swap(A, t2, ctx);
    }

cleanup:

    nmod_mpoly_clear(t1, ctx);
    nmod_mpoly_clear(t2, ctx);

    return success;
}

