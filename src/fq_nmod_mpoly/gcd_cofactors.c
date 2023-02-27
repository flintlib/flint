/*
    Copyright (C) 2019-2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"


int fq_nmod_mpoly_gcd_cofactors(
    fq_nmod_mpoly_t G,
    fq_nmod_mpoly_t Abar,
    fq_nmod_mpoly_t Bbar,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    if (fq_nmod_mpoly_is_zero(A, ctx))
    {
        if (fq_nmod_mpoly_is_zero(B, ctx))
        {
            fq_nmod_mpoly_zero(G, ctx);
            fq_nmod_mpoly_zero(Abar, ctx);
            fq_nmod_mpoly_zero(Bbar, ctx);
            return 1;
        }
        fq_nmod_mpoly_set(G, B, ctx);
        fq_nmod_mpoly_zero(Abar, ctx);
        fq_nmod_mpoly_one(Bbar, ctx);
        if (!n_fq_is_one(G->coeffs + 0, ctx->fqctx))
        {
            fq_nmod_mpoly_scalar_mul_n_fq(Bbar, Bbar, G->coeffs + 0, ctx);
            fq_nmod_mpoly_make_monic(G, G, ctx);
        }
        return 1;
    }

    if (fq_nmod_mpoly_is_zero(B, ctx))
    {
        fq_nmod_mpoly_set(G, A, ctx);
        fq_nmod_mpoly_zero(Bbar, ctx);
        fq_nmod_mpoly_one(Abar, ctx);
        if (!n_fq_is_one(G->coeffs + 0, ctx->fqctx))
        {
            fq_nmod_mpoly_scalar_mul_n_fq(Abar, Abar, G->coeffs + 0, ctx);
            fq_nmod_mpoly_make_monic(G, G, ctx);
        }
        return 1;
    }

    return _fq_nmod_mpoly_gcd_algo(G, Abar, Bbar, A, B, ctx, MPOLY_GCD_USE_ALL);
}

