/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

/*
    Try to set G to the gcd of A and B using Brown's alogrithm M.
    This function switches to a big primes version if needed.
    It should only really fail if the dense size of the inputs is too large.
*/
int fq_nmod_mpoly_gcd_brown(fq_nmod_mpoly_t G,
                            const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    fq_nmod_mpolyd_t Ad, Bd, Gd, Abar, Bbar;
    fq_nmod_mpolyd_ctx_t dctx;
    slong nvars = ctx->minfo->nvars;

    success = 1;

    if (fq_nmod_mpoly_is_zero(A, ctx))
    {
        fq_nmod_mpoly_set(G, B, ctx);
        goto cleanup_stage0;
    }
    if (fq_nmod_mpoly_is_zero(B, ctx))
    {
        fq_nmod_mpoly_set(G, A, ctx);
        goto cleanup_stage0;
    }

    fq_nmod_mpolyd_ctx_init2(dctx, nvars, ctx->fqctx);
    success = fq_nmod_mpolyd_ctx_set_for_gcd(dctx, A, B, ctx);
    if (!success)
    {
        fq_nmod_mpoly_zero(G, ctx);
        goto cleanup_stage1;
    }

    fq_nmod_mpolyd_init(Ad, nvars, ctx->fqctx);
    fq_nmod_mpolyd_init(Bd, nvars, ctx->fqctx);
    fq_nmod_mpolyd_init(Gd, nvars, ctx->fqctx);
    fq_nmod_mpolyd_init(Abar, nvars, ctx->fqctx);
    fq_nmod_mpolyd_init(Bbar, nvars, ctx->fqctx);

    fq_nmod_mpoly_convert_to_fq_nmod_mpolyd(Ad, dctx, A, ctx);
    fq_nmod_mpoly_convert_to_fq_nmod_mpolyd(Bd, dctx, B, ctx);
    success = fq_nmod_mpolyd_gcd_brown_smprime(Gd, Abar, Bbar, Ad, Bd, dctx);
    if (!success)
    {
        fq_nmod_mpoly_convert_to_fq_nmod_mpolyd(Ad, dctx, A, ctx);
        fq_nmod_mpoly_convert_to_fq_nmod_mpolyd(Bd, dctx, B, ctx);
        success = fq_nmod_mpolyd_gcd_brown_lgprime(Gd, Abar, Bbar, Ad, Bd, dctx);
    }
    if (success)
    {
        fq_nmod_mpoly_convert_from_fq_nmod_mpolyd(G, ctx, Gd, dctx);
    }

    fq_nmod_mpolyd_clear(Bbar, ctx->fqctx);
    fq_nmod_mpolyd_clear(Abar, ctx->fqctx);
    fq_nmod_mpolyd_clear(Gd, ctx->fqctx);
    fq_nmod_mpolyd_clear(Bd, ctx->fqctx);
    fq_nmod_mpolyd_clear(Ad, ctx->fqctx);

cleanup_stage1:

    fq_nmod_mpolyd_ctx_clear(dctx);

cleanup_stage0:

    if (success && !fq_nmod_mpoly_is_zero(G, ctx))
    {
        fq_nmod_mpoly_make_monic(G, G, ctx);
    }

    return success;
}
