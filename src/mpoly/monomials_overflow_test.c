/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

int mpoly_monomials_overflow_test(ulong * exps, slong len, flint_bitcnt_t bits,
                                                        const mpoly_ctx_t mctx)
{
    slong N, i;
    N = mpoly_words_per_exp(bits, mctx);

    if (bits <= FLINT_BITS)
    {
        ulong mask = mpoly_overflow_mask_sp(bits);

        for (i = 0; i < len; i++)
            if (mpoly_monomial_overflows(exps + i*N, N, mask))
                return 1;
    }
    else
    {
        for (i = 0; i < len; i++)
            if (mpoly_monomial_overflows_mp(exps + i*N, N, bits))
                return 1;
    }
    return 0;
}
