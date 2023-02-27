/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"


void mpoly_monomials_shift_right_ui(ulong * Aexps, flint_bitcnt_t Abits,
                slong Alength, const ulong * user_exps, const mpoly_ctx_t mctx)
{
    slong i;
    slong N = mpoly_words_per_exp(Abits, mctx);
    ulong * texps;
    TMP_INIT;

    TMP_START;

    texps = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    mpoly_set_monomial_ui(texps, user_exps, Abits, mctx);

    if (Abits <= FLINT_BITS)
    {
#if FLINT_WANT_ASSERT
        ulong mask = mpoly_overflow_mask_sp(Abits);
#endif
        for (i = 0; i < Alength; i++)
        {
            mpoly_monomial_sub(Aexps + N*i, Aexps + N*i, texps, N);
            FLINT_ASSERT(!mpoly_monomial_overflows(Aexps + N*i, N, mask));
        }
    }
    else
    {
        for (i = 0; i < Alength; i++)
        {
            mpoly_monomial_sub_mp(Aexps + N*i, Aexps + N*i, texps, N);
            FLINT_ASSERT(!mpoly_monomial_overflows_mp(Aexps + N*i, N, Abits));
        }
    }

    TMP_END;
}


void mpoly_monomials_shift_right_ffmpz(ulong * Aexps, flint_bitcnt_t Abits,
                slong Alength, const fmpz * user_exps, const mpoly_ctx_t mctx)
{
    slong i;
    slong N = mpoly_words_per_exp(Abits, mctx);
    ulong * texps;
    TMP_INIT;

    TMP_START;

    texps = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    mpoly_set_monomial_ffmpz(texps, user_exps, Abits, mctx);

    if (Abits <= FLINT_BITS)
    {
#if FLINT_WANT_ASSERT
        ulong mask = mpoly_overflow_mask_sp(Abits);
#endif
        for (i = 0; i < Alength; i++)
        {
            mpoly_monomial_sub(Aexps + N*i, Aexps + N*i, texps, N);
            FLINT_ASSERT(!mpoly_monomial_overflows(Aexps + N*i, N, mask));
        }
    }
    else
    {
        for (i = 0; i < Alength; i++)
        {
            mpoly_monomial_sub_mp(Aexps + N*i, Aexps + N*i, texps, N);
            FLINT_ASSERT(!mpoly_monomial_overflows_mp(Aexps + N*i, N, Abits));
        }
    }

    TMP_END;
}


