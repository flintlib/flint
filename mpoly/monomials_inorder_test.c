/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"


int mpoly_monomials_inorder_test(ulong * exps, slong len,
                                   flint_bitcnt_t bits, const mpoly_ctx_t mctx)
{
    slong N, i;
    ulong * cmpmask;

    N = mpoly_words_per_exp(bits, mctx);
    cmpmask = flint_malloc((N + 1)*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, bits, mctx);

    for (i = 0; i + 1 < len; i++)
    {
        if (!mpoly_monomial_gt(exps + i*N, exps + (i + 1)*N, N, cmpmask))
        {
            flint_free(cmpmask);
            return 0;
        }
    }
    flint_free(cmpmask);
    return 1;
}
