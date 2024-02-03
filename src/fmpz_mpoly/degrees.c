/*
    Copyright (C) 2016-2017 William Hart
    Copyright (C) 2017-2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "fmpz_mpoly.h"

int fmpz_mpoly_degrees_fit_si(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
               : mpoly_degrees_fit_si(A->exps, A->length, A->bits, ctx->minfo);
}

void fmpz_mpoly_degrees_fmpz(fmpz ** degs, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
    mpoly_degrees_pfmpz(degs, A->exps, A->length, A->bits, ctx->minfo);
}

void fmpz_mpoly_degrees_si(slong * degs, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
    mpoly_degrees_si(degs, A->exps, A->length, A->bits, ctx->minfo);
}

void fmpz_mpoly_degree_fmpz(fmpz_t deg, const fmpz_mpoly_t A, slong var, const fmpz_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(0 <= var && var < ctx->minfo->nvars);
    mpoly_degree_fmpz(deg, A->exps, A->length, A->bits, var, ctx->minfo);
}

slong fmpz_mpoly_degree_si(const fmpz_mpoly_t A, slong var, const fmpz_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(0 <= var && var < ctx->minfo->nvars);
    return mpoly_degree_si(A->exps, A->length, A->bits, var, ctx->minfo);
}

int fmpz_mpoly_total_degree_fits_si(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_fits_si(A->exps, A->length, A->bits, ctx->minfo);
}

void fmpz_mpoly_total_degree_fmpz(fmpz_t td, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
    mpoly_total_degree_fmpz(td, A->exps, A->length, A->bits, ctx->minfo);
}

slong fmpz_mpoly_total_degree_si(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_si(A->exps, A->length, A->bits, ctx->minfo);
}

void fmpz_mpoly_used_vars(int * used, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    for (i = 0; i < ctx->minfo->nvars; i++)
        used[i] = 0;

    mpoly_used_vars_or(used, A->exps, A->length, A->bits, ctx->minfo);
}
