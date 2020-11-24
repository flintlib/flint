/*
    Copyright (C) 2019-2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"


int nmod_mpoly_gcd_cofactors(
    nmod_mpoly_t G,
    nmod_mpoly_t Abar,
    nmod_mpoly_t Bbar,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    if (nmod_mpoly_is_zero(A, ctx))
    {
        if (nmod_mpoly_is_zero(B, ctx))
        {
            nmod_mpoly_zero(G, ctx);
            nmod_mpoly_zero(Abar, ctx);
            nmod_mpoly_zero(Bbar, ctx);
            return 1;
        }
        nmod_mpoly_set(G, B, ctx);
        nmod_mpoly_zero(Abar, ctx);
        nmod_mpoly_one(Bbar, ctx);
        if (G->coeffs[0] != 1)
        {
            _nmod_vec_scalar_mul_nmod(Bbar->coeffs, Bbar->coeffs,
                                         Bbar->length, G->coeffs[0], ctx->mod);
            _nmod_vec_scalar_mul_nmod(G->coeffs, G->coeffs, G->length,
                                   nmod_inv(G->coeffs[0], ctx->mod), ctx->mod);
        }
        return 1;
    }

    if (nmod_mpoly_is_zero(B, ctx))
    {
        nmod_mpoly_set(G, A, ctx);
        nmod_mpoly_zero(Bbar, ctx);
        nmod_mpoly_one(Abar, ctx);
        if (G->coeffs[0] != 1)
        {
            _nmod_vec_scalar_mul_nmod(Abar->coeffs, Abar->coeffs,
                                         Abar->length, G->coeffs[0], ctx->mod);
            _nmod_vec_scalar_mul_nmod(G->coeffs, G->coeffs, G->length,
                                   nmod_inv(G->coeffs[0], ctx->mod), ctx->mod);
        }
        return 1;
    }

    return _nmod_mpoly_gcd_algo(G, Abar, Bbar, A, B, ctx, MPOLY_GCD_USE_ALL);
}

