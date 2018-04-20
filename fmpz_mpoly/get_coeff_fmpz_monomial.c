/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_get_coeff_fmpz_monomial(fmpz_t c, const fmpz_mpoly_t poly1,
                         const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)
{
    slong i, index, N, nvars = ctx->minfo->nvars;
    ulong * cmpmask, * pexp;
    int exists;
    TMP_INIT;

    if (poly2->length != WORD(1))
    {
        flint_throw(FLINT_ERROR, "poly2 not monomial in fmpz_mpoly_get_coeff_fmpz_monomial");
    }

    TMP_START;

    N = mpoly_words_per_exp(poly1->bits, ctx->minfo);

    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    pexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    mpoly_get_cmpmask(cmpmask, N, poly1->bits, ctx->minfo);
    if (poly2->bits == poly1->bits)
    {
        mpoly_monomial_set(pexp, poly2->exps + N*0, N);
    } else
    {
        slong exp_bits;
        fmpz * texps;

        texps = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
        for (i = 0; i < nvars; i++)
            fmpz_init(texps + i);

        mpoly_get_monomial_ffmpz(texps, poly2->exps + 0, poly2->bits, ctx->minfo);

        exp_bits = mpoly_exp_bits_required_ffmpz(texps, ctx->minfo);

        if (exp_bits > poly1->bits) /* exponent too large to be poly1 exponent */
        {
            fmpz_zero(c);
            for (i = 0; i < nvars; i++)
                fmpz_clear(texps + i);
            goto clean_up;
        }

        mpoly_set_monomial_ffmpz(pexp, texps, poly1->bits, ctx->minfo);

        for (i = 0; i < nvars; i++)
            fmpz_clear(texps + i);
    }

    exists = mpoly_monomial_exists(&index, poly1->exps,
                                       pexp, poly1->length, N, cmpmask);

    if (!exists)
        fmpz_zero(c);
    else
        fmpz_set(c, poly1->coeffs + index);

clean_up:
    TMP_END;
    return;
}
