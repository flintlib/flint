/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "fmpq_mpoly.h"

int fmpq_mpoly_degrees_fit_si(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
{
    return A->zpoly->bits <= FLINT_BITS ? 1
                               : mpoly_degrees_fit_si(A->zpoly->exps,
                                       A->zpoly->length, A->zpoly->bits,
                                                             ctx->zctx->minfo);
}

void fmpq_mpoly_degrees_fmpz(fmpz ** degs, const fmpq_mpoly_t A,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    mpoly_degrees_pfmpz(degs, A->zpoly->exps, A->zpoly->length,
                                             A->zpoly->bits, ctx->zctx->minfo);
}

void fmpq_mpoly_degrees_si(slong * degs, const fmpq_mpoly_t A,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    mpoly_degrees_si(degs, A->zpoly->exps, A->zpoly->length,
                                             A->zpoly->bits, ctx->zctx->minfo);
}

void fmpq_mpoly_degree_fmpz(fmpz_t deg, const fmpq_mpoly_t A, slong var,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    mpoly_degree_fmpz(deg, A->zpoly->exps, A->zpoly->length,
                                        A->zpoly->bits, var, ctx->zctx->minfo);
}

slong fmpq_mpoly_degree_si(const fmpq_mpoly_t A, slong var,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    return mpoly_degree_si(A->zpoly->exps, A->zpoly->length,
                                        A->zpoly->bits, var, ctx->zctx->minfo);
}

int fmpq_mpoly_total_degree_fits_si(const fmpq_mpoly_t A,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_fits_si(A->zpoly->exps, A->zpoly->length,
                                             A->zpoly->bits, ctx->zctx->minfo);
}

void fmpq_mpoly_total_degree_fmpz(fmpz_t tdeg, const fmpq_mpoly_t A,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    mpoly_total_degree_fmpz(tdeg, A->zpoly->exps, A->zpoly->length,
                                             A->zpoly->bits, ctx->zctx->minfo);
}

slong fmpq_mpoly_total_degree_si(const fmpq_mpoly_t A,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_si(A->zpoly->exps, A->zpoly->length,
                                             A->zpoly->bits, ctx->zctx->minfo);
}
