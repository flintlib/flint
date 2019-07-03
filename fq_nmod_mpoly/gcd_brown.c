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
*/
int fq_nmod_mpoly_gcd_brown(
    fq_nmod_mpoly_t G,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong * perm;
    ulong * shift, * stride;
    slong i;
    flint_bitcnt_t new_bits;
    fq_nmod_mpoly_ctx_t uctx;
    fq_nmod_mpolyun_t An, Bn, Gn, Abarn, Bbarn;

    if (fq_nmod_mpoly_is_zero(A, ctx))
    {
        if (fq_nmod_mpoly_is_zero(B, ctx))
        {
            fq_nmod_mpoly_zero(G, ctx);
        }
        else
        {
            fq_nmod_mpoly_make_monic(G, B, ctx);
        }
        return 1;
    }

    if (fq_nmod_mpoly_is_zero(B, ctx))
    {
        fq_nmod_mpoly_make_monic(G, A, ctx);
        return 1;
    }

    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
    {
        return 0;
    }

    perm = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    shift = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    stride = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        perm[i] = i;
        shift[i] = 0;
        stride[i] = 1;
    }

    if (ctx->minfo->nvars == 1)
    {
        fq_nmod_poly_t a, b, g;
        fq_nmod_poly_init(a, ctx->fqctx);
        fq_nmod_poly_init(b, ctx->fqctx);
        fq_nmod_poly_init(g, ctx->fqctx);
        _fq_nmod_mpoly_to_fq_nmod_poly_deflate(a, A, 0, shift, stride, ctx);
        _fq_nmod_mpoly_to_fq_nmod_poly_deflate(b, B, 0, shift, stride, ctx);
        fq_nmod_poly_gcd(g, a, b, ctx->fqctx);
        _fq_nmod_mpoly_from_fq_nmod_poly_inflate(G, A->bits, g, 0, shift, stride, ctx);
        fq_nmod_poly_clear(a, ctx->fqctx);
        fq_nmod_poly_clear(b, ctx->fqctx);
        fq_nmod_poly_clear(g, ctx->fqctx);
        success = 1;
        goto cleanup1;
    }

    new_bits = FLINT_MAX(A->bits, B->bits);

    fq_nmod_mpoly_ctx_init(uctx, ctx->minfo->nvars - 1, ORD_LEX, ctx->fqctx);
    fq_nmod_mpolyun_init(An, new_bits, uctx);
    fq_nmod_mpolyun_init(Bn, new_bits, uctx);
    fq_nmod_mpolyun_init(Gn, new_bits, uctx);
    fq_nmod_mpolyun_init(Abarn, new_bits, uctx);
    fq_nmod_mpolyun_init(Bbarn, new_bits, uctx);

    _fq_nmod_mpoly_to_mpolyun_perm_deflate(An, uctx, A, ctx, perm, shift, stride);
    _fq_nmod_mpoly_to_mpolyun_perm_deflate(Bn, uctx, B, ctx, perm, shift, stride);
    success = fq_nmod_mpolyun_gcd_brown_smprime(Gn, Abarn, Bbarn, An, Bn, uctx->minfo->nvars - 1, uctx);
    if (!success)
    {
        _fq_nmod_mpoly_to_mpolyun_perm_deflate(An, uctx, A, ctx, perm, shift, stride);
        _fq_nmod_mpoly_to_mpolyun_perm_deflate(Bn, uctx, B, ctx, perm, shift, stride);
        success = fq_nmod_mpolyun_gcd_brown_lgprime(Gn, Abarn, Bbarn, An, Bn, uctx->minfo->nvars - 1, uctx);
    }

    if (success)
    {
        _fq_nmod_mpoly_from_mpolyun_perm_inflate(G, new_bits, ctx, Gn, uctx, perm, shift, stride);
        fq_nmod_mpoly_make_monic(G, G, ctx);        
    }

    fq_nmod_mpolyun_clear(An, uctx);
    fq_nmod_mpolyun_clear(Bn, uctx);
    fq_nmod_mpolyun_clear(Gn, uctx);
    fq_nmod_mpolyun_clear(Abarn, uctx);
    fq_nmod_mpolyun_clear(Bbarn, uctx);
    fq_nmod_mpoly_ctx_clear(uctx);

cleanup1:

    flint_free(perm);
    flint_free(shift);
    flint_free(stride);

    return success;
}
