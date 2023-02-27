/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly_factor.h"
#include "fq_nmod_mpoly_factor.h"


/* TODO */
int fq_zech_mpoly_factor_algo(
    fq_zech_mpoly_factor_t f,
    const fq_zech_mpoly_t A,
    const fq_zech_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong i;
    fq_nmod_mpoly_ctx_t ctx2;
    fq_nmod_mpoly_t A2;
    fq_nmod_mpoly_factor_t f2;

    *ctx2->minfo = *ctx->minfo;
    *ctx2->fqctx = *ctx->fqctx->fq_nmod_ctx;

    fq_nmod_mpoly_init(A2, ctx2);
    fq_nmod_mpoly_factor_init(f2, ctx2);

    _fq_zech_mpoly_get_fq_nmod_mpoly(A2, ctx2, A, ctx);
    success = fq_nmod_mpoly_factor_algo(f2, A2, ctx2, algo);
    if (success)
    {
        fq_zech_set_fq_nmod(f->constant, f2->constant, ctx->fqctx);
        fq_zech_mpoly_factor_fit_length(f, f2->num, ctx);
        for (i = 0; i < f2->num; i++)
        {
            _fq_zech_mpoly_set_fq_nmod_mpoly(f->poly + i, ctx, f2->poly + i, ctx2);
            fmpz_swap(f->exp + i, f2->exp + i);
        }
        f->num = f2->num;
    }

    fq_nmod_mpoly_clear(A2, ctx2);
    fq_nmod_mpoly_factor_clear(f2, ctx2);

    return success;
}


int fq_zech_mpoly_factor(
    fq_zech_mpoly_factor_t f,
    const fq_zech_mpoly_t A,
    const fq_zech_mpoly_ctx_t ctx)
{
    return fq_zech_mpoly_factor_algo(f, A, ctx, MPOLY_FACTOR_USE_ALL);
}


int fq_zech_mpoly_factor_zassenhaus(
    fq_zech_mpoly_factor_t f,
    const fq_zech_mpoly_t A,
    const fq_zech_mpoly_ctx_t ctx)
{
    return fq_zech_mpoly_factor_algo(f, A, ctx, MPOLY_FACTOR_USE_ZAS);
}


int fq_zech_mpoly_factor_wang(
    fq_zech_mpoly_factor_t f,
    const fq_zech_mpoly_t A,
    const fq_zech_mpoly_ctx_t ctx)
{
    return fq_zech_mpoly_factor_algo(f, A, ctx, MPOLY_FACTOR_USE_WANG);
}


int fq_zech_mpoly_factor_zippel(
    fq_zech_mpoly_factor_t f,
    const fq_zech_mpoly_t A,
    const fq_zech_mpoly_ctx_t ctx)
{
    return fq_zech_mpoly_factor_algo(f, A, ctx, MPOLY_FACTOR_USE_ZIP);
}

