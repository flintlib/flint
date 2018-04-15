/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void fmpq_mpoly_get_coeff_monomial(fmpq_t c, const fmpq_mpoly_t poly1,
                         const fmpq_mpoly_t poly2, const fmpq_mpoly_ctx_t ctx)
{


    slong i, index, N, nvars = ctx->zctx->minfo->nvars;
    ulong * cmpmask, * pexp;
    int exists;
    TMP_INIT;

    if (poly2->zpoly->length != WORD(1))
    {
        flint_throw(FLINT_ERROR, "poly2 not monomial in fmpq_mpoly_get_coeff_monomial");
    }

    TMP_START;

    N = mpoly_words_per_exp(poly1->zpoly->bits, ctx->zctx->minfo);

    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    pexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    mpoly_get_cmpmask(cmpmask, N, poly1->zpoly->bits, ctx->zctx->minfo);
    if (poly2->zpoly->bits == poly1->zpoly->bits)
    {
        mpoly_monomial_set(pexp, poly2->zpoly->exps + N*0, N);
    } else
    {
        slong exp_bits;
        fmpz * exps;

        exps = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
        for (i = 0; i < nvars; i++)
            fmpz_init(exps + i);

        mpoly_get_monomial_ffmpz(exps, poly2->zpoly->exps + 0,
                                         poly2->zpoly->bits, ctx->zctx->minfo);

        exp_bits = mpoly_exp_bits_required_ffmpz(exps, ctx->zctx->minfo);

        if (exp_bits > poly1->zpoly->bits) /* exponent too large to be poly1 exponent */
        {
            fmpq_zero(c);
            for (i = 0; i < nvars; i++)
                fmpz_clear(exps + i);
            goto clean_up;
        }

        mpoly_set_monomial_ffmpz(pexp, exps, poly1->zpoly->bits, ctx->zctx->minfo);

        for (i = 0; i < nvars; i++)
            fmpz_clear(exps + i);
    }

    exists = mpoly_monomial_exists(&index, poly1->zpoly->exps,
                                       pexp, poly1->zpoly->length, N, cmpmask);

    if (!exists)
    {
        fmpq_zero(c);
    } else
    {
        fmpq_mul_fmpz(c, poly1->content, poly1->zpoly->coeffs + index);
    }

clean_up:
    TMP_END;
    return;
}
