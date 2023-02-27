/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


int fmpq_mpoly_gcd_cofactors(fmpq_mpoly_t G, fmpq_mpoly_t Abar,
              fmpq_mpoly_t Bbar, const fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpq_t cAbar, cBbar;
    int success;

    success = fmpz_mpoly_gcd_cofactors(G->zpoly, Abar->zpoly, Bbar->zpoly,
                                                A->zpoly, B->zpoly, ctx->zctx);
    if (!success)
        return 0;

    fmpq_init(cAbar);
    fmpq_init(cBbar);

    if (G->zpoly->length > 0)
    {
        fmpq_mul_fmpz(cAbar, A->content, G->zpoly->coeffs + 0);
        fmpq_mul_fmpz(cBbar, B->content, G->zpoly->coeffs + 0);
        fmpz_set(fmpq_denref(G->content), G->zpoly->coeffs + 0);
        fmpz_one(fmpq_numref(G->content));
    }
    else
    {
        fmpq_zero(G->content);
    }

    fmpq_swap(Abar->content, cAbar);
    fmpq_swap(Bbar->content, cBbar);

    fmpq_clear(cAbar);
    fmpq_clear(cBbar);

    return 1;
}

