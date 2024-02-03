/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "fq_zech_mpoly.h"

int fq_zech_mpoly_degrees_fit_si(const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
               : mpoly_degrees_fit_si(A->exps, A->length, A->bits, ctx->minfo);
}

void fq_zech_mpoly_degrees_fmpz(fmpz ** degs, const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
{
    mpoly_degrees_pfmpz(degs, A->exps, A->length, A->bits, ctx->minfo);
}

void fq_zech_mpoly_degrees_si(slong * degs, const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
{
    mpoly_degrees_si(degs, A->exps, A->length, A->bits, ctx->minfo);
}

void fq_zech_mpoly_degree_fmpz(fmpz_t deg, const fq_zech_mpoly_t A, slong var, const fq_zech_mpoly_ctx_t ctx)
{
    mpoly_degree_fmpz(deg, A->exps, A->length, A->bits, var, ctx->minfo);
}

slong fq_zech_mpoly_degree_si(const fq_zech_mpoly_t A, slong var, const fq_zech_mpoly_ctx_t ctx)
{
    return mpoly_degree_si(A->exps, A->length, A->bits, var, ctx->minfo);
}

int fq_zech_mpoly_total_degree_fits_si(const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_fits_si(A->exps, A->length, A->bits, ctx->minfo);
}

void fq_zech_mpoly_total_degree_fmpz(fmpz_t td, const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
{
    mpoly_total_degree_fmpz(td, A->exps, A->length, A->bits, ctx->minfo);
}

slong fq_zech_mpoly_total_degree_si(const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_si(A->exps, A->length, A->bits, ctx->minfo);
}
