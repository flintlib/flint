/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

/*
    It is assumed that A is not zero, and B is a monomial.
    This function is always successful and always returns 1.
    G is packed into bits = Gbits
*/
int _fmpz_mpoly_gcd_monomial(fmpz_mpoly_t G, flint_bitcnt_t Gbits, 
                                const fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_t g;
    fmpz * minAfields, * minAdegs, * minBdegs;
    TMP_INIT;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length == 1);

    TMP_START;

    /* get the field-wise minimum of A */
    minAfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
        fmpz_init(minAfields + i);
    mpoly_min_fields_fmpz(minAfields, A->exps, A->length, A->bits, ctx->minfo);

    /* unpack to get the min degrees of each variable in A */
    minAdegs = (fmpz *) TMP_ALLOC(ctx->minfo->nvars*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_init(minAdegs + i);
    mpoly_get_monomial_ffmpz_unpacked_ffmpz(minAdegs, minAfields, ctx->minfo);

    /* get the degree of each variable in B */
    minBdegs = (fmpz *) TMP_ALLOC(ctx->minfo->nvars*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_init(minBdegs + i);
    mpoly_get_monomial_ffmpz(minBdegs, B->exps, B->bits, ctx->minfo);

    /* compute the degree of each variable in G */
    _fmpz_vec_min_inplace(minBdegs, minAdegs, ctx->minfo->nvars);

    fmpz_mpoly_fit_length(G, 1, ctx);
    fmpz_mpoly_fit_bits(G, Gbits, ctx);
    G->bits = Gbits;
    mpoly_set_monomial_ffmpz(G->exps, minBdegs, Gbits, ctx->minfo);

    fmpz_init_set(g, B->coeffs + 0);
    _fmpz_vec_content_chained(g, A->coeffs, A->length);
    fmpz_swap(G->coeffs + 0, g);
    fmpz_clear(g);

    _fmpz_mpoly_set_length(G, 1, ctx);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(minAfields + i);
    }
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        fmpz_clear(minAdegs + i);
        fmpz_clear(minBdegs + i);
    }

    TMP_END;
    return 1;
}
