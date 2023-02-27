/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"


static void _fq_zech_vec_scalar_div_fq_zech(
    fq_zech_struct * A,
    const fq_zech_struct * B,
    slong len,
    const fq_zech_t c,
    const fq_zech_ctx_t fqctx)
{
    slong i;
    fq_zech_t cinv;

    fq_zech_init(cinv, fqctx);
    fq_zech_inv(cinv, c, fqctx);

    for (i = 0; i < len; i++)
        fq_zech_mul(A + i, B + i, cinv, fqctx);

    fq_zech_clear(cinv, fqctx);
}

int fq_zech_mpoly_gcd_cofactors(
    fq_zech_mpoly_t G,
    fq_zech_mpoly_t Abar,
    fq_zech_mpoly_t Bbar,
    const fq_zech_mpoly_t A,
    const fq_zech_mpoly_t B,
    const fq_zech_mpoly_ctx_t ctx)
{
    int success;
    fq_nmod_mpoly_ctx_t ctx2;
    fq_nmod_mpoly_t A2, B2, G2, Abar2, Bbar2;

    if (A->length == 0)
    {
        if (B->length == 0)
        {
            fq_zech_mpoly_zero(G, ctx);
            fq_zech_mpoly_zero(Abar, ctx);
            fq_zech_mpoly_zero(Bbar, ctx);
            return 1;
        }
        fq_zech_mpoly_set(G, B, ctx);
        fq_zech_mpoly_zero(Abar, ctx);
        fq_zech_mpoly_one(Bbar, ctx);
        if (!fq_zech_is_one(G->coeffs + 0, ctx->fqctx))
        {
            _fq_zech_vec_scalar_mul_fq_zech(Bbar->coeffs, Bbar->coeffs,
                                      Bbar->length, G->coeffs + 0, ctx->fqctx);
            _fq_zech_vec_scalar_div_fq_zech(G->coeffs, G->coeffs,
                                         G->length, G->coeffs + 0, ctx->fqctx);
        }
        return 1;
    }

    if (B->length == 0)
    {
        fq_zech_mpoly_set(G, A, ctx);
        fq_zech_mpoly_zero(Bbar, ctx);
        fq_zech_mpoly_one(Abar, ctx);
        if (!fq_zech_is_one(G->coeffs + 0, ctx->fqctx))
        {
            _fq_zech_vec_scalar_mul_fq_zech(Abar->coeffs, Abar->coeffs,
                                      Abar->length, G->coeffs + 0, ctx->fqctx);
            _fq_zech_vec_scalar_div_fq_zech(G->coeffs, G->coeffs,
                                         G->length, G->coeffs + 0, ctx->fqctx);
        }
        return 1;
    }

    *ctx2->minfo = *ctx->minfo;
    *ctx2->fqctx = *ctx->fqctx->fq_nmod_ctx;

    fq_nmod_mpoly_init(A2, ctx2);
    fq_nmod_mpoly_init(B2, ctx2);
    fq_nmod_mpoly_init(G2, ctx2);
    fq_nmod_mpoly_init(Abar2, ctx2);
    fq_nmod_mpoly_init(Bbar2, ctx2);

    _fq_zech_mpoly_get_fq_nmod_mpoly(A2, ctx2, A, ctx);
    _fq_zech_mpoly_get_fq_nmod_mpoly(B2, ctx2, B, ctx);
    success = fq_nmod_mpoly_gcd_cofactors(G2, Abar2, Bbar2, A2, B2, ctx2);
    if (success)
    {
        _fq_zech_mpoly_set_fq_nmod_mpoly(G, ctx, G2, ctx2);
        _fq_zech_mpoly_set_fq_nmod_mpoly(Abar, ctx, Abar2, ctx2);
        _fq_zech_mpoly_set_fq_nmod_mpoly(Bbar, ctx, Bbar2, ctx2);
    }

    fq_nmod_mpoly_clear(A2, ctx2);
    fq_nmod_mpoly_clear(B2, ctx2);
    fq_nmod_mpoly_clear(G2, ctx2);
    fq_nmod_mpoly_clear(Abar2, ctx2);
    fq_nmod_mpoly_clear(Bbar2, ctx2);

    return success;
}
