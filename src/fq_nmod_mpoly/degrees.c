/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "fq_nmod_mpoly.h"

int fq_nmod_mpoly_degrees_fit_si(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
               : mpoly_degrees_fit_si(A->exps, A->length, A->bits, ctx->minfo);
}

void fq_nmod_mpoly_degrees_fmpz(fmpz ** degs, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    mpoly_degrees_pfmpz(degs, A->exps, A->length, A->bits, ctx->minfo);
}

void fq_nmod_mpoly_degrees_si(slong * degs, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    mpoly_degrees_si(degs, A->exps, A->length, A->bits, ctx->minfo);
}

void fq_nmod_mpoly_degree_fmpz(fmpz_t deg, const fq_nmod_mpoly_t A, slong var, const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(0 <= var && var < ctx->minfo->nvars);
    mpoly_degree_fmpz(deg, A->exps, A->length, A->bits, var, ctx->minfo);
}

slong fq_nmod_mpoly_degree_si(const fq_nmod_mpoly_t A, slong var, const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(0 <= var && var < ctx->minfo->nvars);
    return mpoly_degree_si(A->exps, A->length, A->bits, var, ctx->minfo);
}

int fq_nmod_mpoly_total_degree_fits_si(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_fits_si(A->exps, A->length, A->bits, ctx->minfo);
}

void fq_nmod_mpoly_total_degree_fmpz(fmpz_t td, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    mpoly_total_degree_fmpz(td, A->exps, A->length, A->bits, ctx->minfo);
}

slong fq_nmod_mpoly_total_degree_si(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_si(A->exps, A->length, A->bits, ctx->minfo);
}

void fq_nmod_mpoly_used_vars(int * used, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    for (i = 0; i < ctx->minfo->nvars; i++)
        used[i] = 0;

    mpoly_used_vars_or(used, A->exps, A->length, A->bits, ctx->minfo);
}
