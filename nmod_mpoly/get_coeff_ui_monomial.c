/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

ulong nmod_mpoly_get_coeff_ui_monomial(const nmod_mpoly_t A,
                             const nmod_mpoly_t M, const nmod_mpoly_ctx_t ctx)
{
    slong i, index, N, nvars = ctx->minfo->nvars;
    ulong * cmpmask, * pexp;
    ulong c;
    int exists;
    TMP_INIT;

    if (M->length != WORD(1))
    {
        flint_throw(FLINT_ERROR, "M not monomial in nmod_mpoly_get_coeff_fmpz_monomial");
    }

    TMP_START;

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    pexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    mpoly_get_cmpmask(cmpmask, N, A->bits, ctx->minfo);
    if (M->bits == A->bits)
    {
        mpoly_monomial_set(pexp, M->exps + N*0, N);
    }
    else
    {
        mp_bitcnt_t exp_bits;
        fmpz * texps;

        texps = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
        for (i = 0; i < nvars; i++)
            fmpz_init(texps + i);

        mpoly_get_monomial_ffmpz(texps, M->exps + 0, M->bits, ctx->minfo);

        exp_bits = mpoly_exp_bits_required_ffmpz(texps, ctx->minfo);

        if (exp_bits > A->bits) /* exponent too large to be A exponent */
        {
            c = UWORD(0);
            for (i = 0; i < nvars; i++)
                fmpz_clear(texps + i);
            goto clean_up;
        }

        mpoly_set_monomial_ffmpz(pexp, texps, A->bits, ctx->minfo);

        for (i = 0; i < nvars; i++)
            fmpz_clear(texps + i);
    }

    exists = mpoly_monomial_exists(&index, A->exps, pexp, A->length, N, cmpmask);

    if (!exists)
        c = UWORD(0);
    else
        c = A->coeffs[index];

clean_up:
    TMP_END;
    return c;
}
