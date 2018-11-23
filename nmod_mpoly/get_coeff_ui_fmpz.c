/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

ulong _nmod_mpoly_get_coeff_ui_fmpz(const nmod_mpoly_t A,
                                  const fmpz * exp, const nmod_mpoly_ctx_t ctx)
{
    slong N, index;
    mp_bitcnt_t exp_bits;
    ulong * cmpmask, * packed_exp;
    ulong c;
    int exists;
    TMP_INIT;

    exp_bits = mpoly_exp_bits_required_ffmpz(exp, ctx->minfo);

    if (exp_bits > A->bits) /* exponent too large to be poly exponent */
    {
        return UWORD(0);
    }

    TMP_START;
   
    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, A->bits, ctx->minfo);

    packed_exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_set_monomial_ffmpz(packed_exp, exp, A->bits, ctx->minfo);

    exists = mpoly_monomial_exists(&index, A->exps,
                                  packed_exp, A->length, N, cmpmask);

    if (!exists)
        c = UWORD(0);
    else
        c = A->coeffs[index];

    TMP_END;
    return c;
}

ulong nmod_mpoly_get_coeff_ui_fmpz(const nmod_mpoly_t A,
                                fmpz * const * exp, const nmod_mpoly_ctx_t ctx)
{
    slong i, nvars = ctx->minfo->nvars;
    fmpz * newexp;
    ulong c;
    TMP_INIT;

    TMP_START;
    newexp = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (i = 0; i < nvars; i++)
    {
        fmpz_init(newexp + i);
        fmpz_set(newexp + i, exp[i]);
    }

    c = _nmod_mpoly_get_coeff_ui_fmpz(A, newexp, ctx);

    for (i = 0; i < nvars; i++)
        fmpz_clear(newexp + i);

    TMP_END;
    return c;
}
