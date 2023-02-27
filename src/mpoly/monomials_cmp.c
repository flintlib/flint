/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/* this file does not need to change with new orderings */

static int _mpoly_monomials_cmp_repack_bits(
    const ulong * Aexps, flint_bitcnt_t Abits,
    const ulong * Bexps, flint_bitcnt_t Bbits,
    slong length,
    const mpoly_ctx_t mctx)
{
    int cmp = 0;
    ulong * newAexps, * cmpmask;
    slong NA = mpoly_words_per_exp(Abits, mctx);
    slong NB = mpoly_words_per_exp(Bbits, mctx);
    const slong max_limit = 32;
    slong i, j, limit;
    TMP_INIT;

    FLINT_ASSERT(Abits <= Bbits);
    FLINT_ASSERT(NA <= NB);

    TMP_START;

    /* repack first monomial, then the next two, then the next four, ... */

    cmpmask = (ulong *) TMP_ALLOC(NB*sizeof(ulong));
    newAexps = (ulong *) TMP_ALLOC(max_limit*NB*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, NB, Bbits, mctx);

    i = 0;
    limit = 1;

try_again:

    if (i + limit <= length)
    {
        FLINT_ASSERT(limit <= max_limit);
        mpoly_repack_monomials(newAexps, Bbits, Aexps + NA*i, Abits, limit, mctx);
        for (j = 0; j < limit; j++, i++)
        {
            cmp = mpoly_monomial_cmp(newAexps + NB*j, Bexps + NB*i, NB, cmpmask);
            if (cmp != 0)
                goto cleanup;
        }

        limit = FLINT_MIN(2*limit, max_limit);
        goto try_again;
    }

    FLINT_ASSERT(length - i <= max_limit);
    mpoly_repack_monomials(newAexps, Bbits, Aexps + NA*i, Abits, length - i, mctx);
    for (j = 0; i < length; j++, i++)
    {
        cmp = mpoly_monomial_cmp(newAexps + NB*j, Bexps + NB*i, NB, cmpmask);
        if (cmp != 0)
            goto cleanup;
    }

cleanup:

    TMP_END;

    return cmp;
}


/* cmp two exponent matrices of the same length */
int mpoly_monomials_cmp(
    const ulong * Aexps, flint_bitcnt_t Abits,
    const ulong * Bexps, flint_bitcnt_t Bbits,
    slong length,
    const mpoly_ctx_t mctx)
{
    int cmp = 0;
    slong i, N;
    ulong * cmpmask, cmpmask1;
    TMP_INIT;

    FLINT_ASSERT(length >= 0);

    if (Abits < Bbits)
    {
        return _mpoly_monomials_cmp_repack_bits(Aexps, Abits, Bexps, Bbits,
                                                                 length, mctx);
    }
    else if (Abits > Bbits)
    {
        return -_mpoly_monomials_cmp_repack_bits(Bexps, Bbits, Aexps, Abits,
                                                                 length, mctx);
    }

    N = mpoly_words_per_exp(Abits, mctx);

    if (N == 1)
    {
        mpoly_get_cmpmask(&cmpmask1, 1, Abits, mctx);
        for (i = 0; i < length && cmp == 0; i++)
            cmp = mpoly_monomial_cmp1(Aexps[i], Bexps[i], cmpmask1);
    }
    else
    {
        TMP_START;
        cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
        mpoly_get_cmpmask(cmpmask, N, Abits, mctx);
        for (i = 0; i < length && cmp == 0; i++)
            cmp = mpoly_monomial_cmp(Aexps + N*i, Bexps + N*i, N, cmpmask);
        TMP_END;
    }

    return cmp;
}
