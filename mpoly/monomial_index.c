/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/* this file does not need to change with new orderings */


slong mpoly_monomial_index1_nomask(ulong * Aexps, slong Alen, ulong e)
{
    slong start = 0, i, stop = Alen;

again:

    if (stop - start < 8)
    {
        for (i = start; i < stop; i++)
        {
            if (Aexps[i] == e)
                return i;
        }
        return -1;
    }

    i = (start + stop)/2;

    FLINT_ASSERT(Aexps[start] > Aexps[i]);
    FLINT_ASSERT(stop >= Alen || Aexps[stop] < Aexps[i]);

    if (Aexps[i] < e)
    {
        stop = i;
        goto again;
    }
    else if (Aexps[i] > e)
    {
        start = i;
        goto again;
    }

    return i;
}


slong mpoly_monomial_index_ui(const ulong * Aexps, flint_bitcnt_t Abits,
                      slong Alength, const ulong * exp, const mpoly_ctx_t mctx)
{
    slong N, index;
    ulong * cmpmask, * packed_exp;
    flint_bitcnt_t exp_bits;
    int exists;
    TMP_INIT;

    exp_bits = mpoly_exp_bits_required_ui(exp, mctx);
    if (exp_bits > Abits)
        return -WORD(1);

    TMP_START;

    N = mpoly_words_per_exp(Abits, mctx);

    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, Abits, mctx);

    packed_exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_set_monomial_ui(packed_exp, exp, Abits, mctx);

    exists = mpoly_monomial_exists(&index, Aexps, packed_exp, Alength, N, cmpmask);

    TMP_END;

    if (!exists)
        return -WORD(1);
    else
        return index;
}

slong mpoly_monomial_index_pfmpz(const ulong * Aexps, flint_bitcnt_t Abits,
                     slong Alength, fmpz * const * exp, const mpoly_ctx_t mctx)
{
    slong N, index;
    ulong * cmpmask, * packed_exp;
    flint_bitcnt_t exp_bits;
    int exists;
    TMP_INIT;

    exp_bits = mpoly_exp_bits_required_pfmpz(exp, mctx);
    if (exp_bits > Abits)
        return -WORD(1);

    TMP_START;

    N = mpoly_words_per_exp(Abits, mctx);

    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, Abits, mctx);

    packed_exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_set_monomial_pfmpz(packed_exp, exp, Abits, mctx);

    exists = mpoly_monomial_exists(&index, Aexps, packed_exp, Alength, N, cmpmask);

    TMP_END;

    if (!exists)
        return -WORD(1);
    else
        return index;
}

slong mpoly_monomial_index_monomial(const ulong * Aexps, flint_bitcnt_t Abits,
                        slong Alength, const ulong * Mexp, flint_bitcnt_t Mbits,
                                                        const mpoly_ctx_t mctx)
{
    slong N, index;
    ulong * cmpmask, * pexp;
    int exists, could_repack;
    TMP_INIT;

    TMP_START;

    N = mpoly_words_per_exp(Abits, mctx);

    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, Abits, mctx);

    if (Mbits == Abits)
    {
        exists = mpoly_monomial_exists(&index, Aexps, Mexp, Alength, N, cmpmask);
    }
    else
    {
        pexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
        could_repack = mpoly_repack_monomials(pexp, Abits, Mexp, Mbits, 1, mctx);
        if (!could_repack)
        {
            FLINT_ASSERT(Mbits > Abits);
            exists = 0;
            index = -WORD(1);
            goto clean_up;
        }

        exists = mpoly_monomial_exists(&index, Aexps, pexp, Alength, N, cmpmask);
    }

clean_up:

    TMP_END;

    if (!exists)
        return -WORD(1);
    else
        return index;

}
