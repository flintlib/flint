/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "nmod_mpoly.h"


void nmod_mpoly_pow(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                                           slong k, const nmod_mpoly_ctx_t ctx)
{
    slong i, bits, exp_bits, N, len1 = 0;
    ulong * max_degs2;
    ulong max = 0;
    ulong * exp2 = poly2->exps;
    int free2 = 0;

    TMP_INIT;

    if (k == 0)
    {
        nmod_mpoly_set_ui(poly1, 1, ctx);
        return;
    }

    if (poly2->length == 0)
    {
        nmod_mpoly_zero(poly1, ctx);
        return;
    }

    if (k == 1)
    {
        nmod_mpoly_set(poly1, poly2, ctx);
        return;
    }

    if (k == 2)
    {
        nmod_mpoly_mul_johnson(poly1, poly2, poly2, ctx);
        return;
    }

    if (poly2->length > 1)
    {
        nmod_mpoly_pow_rmul(poly1, poly2, k, ctx);
        return;
    }

    /* code for powering a monomial */

    TMP_START;

    max_degs2 = (ulong *) TMP_ALLOC(ctx->n*sizeof(ulong));
    mpoly_max_degrees(max_degs2, poly2->exps, poly2->length, poly2->bits, ctx->n);

    for (i = 0; i < ctx->n; i++)
    {
        if (max_degs2[i] > max)
            max = max_degs2[i];
    }

    if (FLINT_BIT_COUNT(max) + FLINT_BIT_COUNT(k) > FLINT_BITS || 0 > (slong) (k*max))
        flint_throw(FLINT_EXPOF, "Exponent overflow in nmod_mpoly_pow");

    bits = FLINT_BIT_COUNT(k*max);
    exp_bits = FLINT_MAX(WORD(8), bits + 1);
    exp_bits = FLINT_MAX(exp_bits, poly2->bits);
    exp_bits = mpoly_optimize_bits(exp_bits, ctx->n);

    N = words_per_exp(ctx->n, exp_bits);

    if (exp_bits > poly2->bits)
    {
        free2 = 1;
        exp2 = (ulong *) flint_malloc(N*poly2->length*sizeof(ulong));
        mpoly_unpack_monomials(exp2, exp_bits, poly2->exps, poly2->bits,
                                                        poly2->length, ctx->n);
    }

    nmod_mpoly_fit_length(poly1, 1, ctx);
    nmod_mpoly_fit_bits(poly1, exp_bits, ctx);
    poly1->bits = exp_bits;

    mpoly_monomial_mul_si(poly1->exps + N*0, poly2->exps + N*0, N, k);
    poly1->coeffs[0] = nmod_pow_ui(poly2->coeffs[0], k, ctx->ffinfo->mod);
    len1 = (poly1->coeffs[0] != 0);

    if (free2)
        flint_free(exp2);

    _nmod_mpoly_set_length(poly1, len1, ctx);

    TMP_END;
}
