/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/* this file does not need to change with new orderings */

/*
    Convert B to A using the variable permutation vector perm.
    A must be constructed with bits = Abits.

    operation on each term:

        for 0 <= l < n
            Aexp[l] = shift[l]

        for 0 <= k < m + 2
            l = perm[k]
            Aexp[l] += scale[l]*Bexp[k]
*/
void mpoly_from_mpolyl_perm_inflate(
    ulong * Bexps,
    flint_bitcnt_t Bbits,
    const mpoly_ctx_t Bctx,
    ulong * Aexps,
    flint_bitcnt_t Abits,
    const mpoly_ctx_t Actx,
    slong length,
    const slong * perm,
    const ulong * shift,
    const ulong * stride)
{
    slong i, k, l;
    slong n = Bctx->nvars;
    slong m = Actx->nvars;
    slong NB = mpoly_words_per_exp_sp(Bbits, Bctx);
    slong NA = mpoly_words_per_exp_sp(Abits, Actx);
    ulong * aexps, * bexps;
    TMP_INIT;

    FLINT_ASSERT(Actx->ord == ORD_LEX);
    FLINT_ASSERT(Bbits <= FLINT_BITS);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(m <= n);

    TMP_START;

    aexps = TMP_ARRAY_ALLOC(m + n, ulong);
    bexps = aexps + m;

    for (i = 0; i < length; i++)
    {
	    mpoly_get_monomial_ui(aexps, Aexps + NA*i, Abits, Actx);

        for (l = 0; l < n; l++)
        {
            bexps[l] = shift[l];
        }
        for (k = 0; k < m; k++)
        {
            l = perm[k];
            bexps[l] += stride[l]*aexps[k];
        }

        mpoly_set_monomial_ui(Bexps + NB*i, bexps, Bbits, Bctx);
    }

    TMP_END;
}



void mpoly_to_mpolyl_perm_deflate(
    ulong * Aexps,
    flint_bitcnt_t Abits,
    const mpoly_ctx_t Actx,
    ulong * Bexps,
    flint_bitcnt_t Bbits,
    const mpoly_ctx_t Bctx,
    slong length,
    const slong * perm,
    const ulong * shift,
    const ulong * stride)
{
    slong j, k, l;
    slong m = Actx->nvars;
    slong n = Bctx->nvars;
    slong NA = mpoly_words_per_exp_sp(Abits, Actx);
    slong NB = mpoly_words_per_exp_sp(Bbits, Bctx);
    ulong * aexps, * bexps;
    TMP_INIT;

    FLINT_ASSERT(Actx->ord == ORD_LEX);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(Bbits <= FLINT_BITS);
    FLINT_ASSERT(m <= n);

    TMP_START;

    aexps = TMP_ARRAY_ALLOC(m + n, ulong);
    bexps = aexps + m;

    for (j = 0; j < length; j++)
    {
        mpoly_get_monomial_ui(bexps, Bexps + NB*j, Bbits, Bctx);
        for (k = 0; k < m; k++)
        {
            l = perm[k];
            if (stride[l] == 1)
            {
                aexps[k] = (bexps[l] - shift[l]);
            }
            else
            {
                FLINT_ASSERT(stride[l] != 0);
                FLINT_ASSERT(((bexps[l] - shift[l]) % stride[l]) == 0);
                aexps[k] = (bexps[l] - shift[l]) / stride[l];
            }
        }

        mpoly_set_monomial_ui(Aexps + NA*j, aexps, Abits, Actx);
    }

    TMP_END;
}

