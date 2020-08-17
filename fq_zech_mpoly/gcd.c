/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"


int fq_zech_mpoly_gcd(
    fq_zech_mpoly_t G,
    const fq_zech_mpoly_t A,
    const fq_zech_mpoly_t B,
    const fq_zech_mpoly_ctx_t ctx)
{
    int success;
    fq_nmod_mpoly_ctx_t ctx2;
    fq_nmod_mpoly_t A2, B2, G2;

    if (fq_zech_mpoly_is_zero(A, ctx))
    {
        if (fq_zech_mpoly_is_zero(B, ctx))
            fq_zech_mpoly_zero(G, ctx);
        else
            fq_zech_mpoly_make_monic(G, B, ctx);
        return 1;
    }

    if (fq_zech_mpoly_is_zero(B, ctx))
    {
        fq_zech_mpoly_make_monic(G, A, ctx);
        return 1;
    }

    *ctx2->minfo = *ctx->minfo;
    *ctx2->fqctx = *ctx->fqctx->fq_nmod_ctx;

    fq_nmod_mpoly_init(A2, ctx2);
    fq_nmod_mpoly_init(B2, ctx2);
    fq_nmod_mpoly_init(G2, ctx2);

    _fq_zech_mpoly_get_fq_nmod_mpoly(A2, ctx2, A, ctx);
    _fq_zech_mpoly_get_fq_nmod_mpoly(B2, ctx2, B, ctx);
    success = fq_nmod_mpoly_gcd(G2, A2, B2, ctx2);
    if (success)
        _fq_zech_mpoly_set_fq_nmod_mpoly(G, ctx, G2, ctx2);

    fq_nmod_mpoly_clear(A2, ctx2);
    fq_nmod_mpoly_clear(B2, ctx2);
    fq_nmod_mpoly_clear(G2, ctx2);

    return success;
}

