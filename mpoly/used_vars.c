/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/* this file DOES need to change with new orderings !!! */

static void mpoly_used_vars_or_sp(
    int * used,
    const ulong * Aexps, slong Alen, flint_bitcnt_t Abits,
    const mpoly_ctx_t mctx)
{
    slong Ai, i, j, m, dir;
    slong N = mpoly_words_per_exp(Abits, mctx);
    slong nvars = mctx->nvars;
    ulong * t;
    slong Aimod, Amodulus = n_sqrt(Alen);
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - Abits);
    ulong u;
    flint_bitcnt_t shift;
    const ulong * exp2;
    TMP_INIT;

    TMP_START;
    t = TMP_ARRAY_ALLOC(N, ulong);
    mpoly_monomial_zero(t, N);

    m = 0;

    for (Aimod = 0; Aimod < Amodulus; Aimod++)
    {
        while (m < nvars && used[m])
            m++;

        /* all variables in [0, m) are used */
        if (m >= nvars)
            goto cleanup;

        for (Ai = Aimod; Ai < Alen; Ai += Amodulus)
        {
            for (j = 0; j < N; j++)
                t[j] |= Aexps[N*Ai + j];
        }

        j = mctx->rev ? 0 : nvars - 1;
        dir = mctx->rev ? 1 : -1;

        /* unpack vec pt */

        exp2 = t;

        i = 0;
        u = *exp2++;
        shift = 0;

        FLINT_ASSERT(0 <= j && j < nvars);
        used[j] |= (u & mask) != 0;
        j += dir;

        u = u >> Abits;      /* number of bits to encode 0th field */
        shift += Abits;      /* number of bits to encode 0th field */
        while (++i < nvars)
        {
            if (shift + Abits > FLINT_BITS)
            {
                u = *exp2++;
                shift = 0;
            }

            FLINT_ASSERT(0 <= j && j < nvars);
            used[j] |= (u & mask) != 0;
            j += dir;

            u = u >> Abits;      /* number of bits to encode ith field */
            shift += Abits;      /* number of bits to encode ith field */
        }
    }

cleanup:

    TMP_END;
}

static void mpoly_used_vars_or_mp(
    int * used,
    const ulong * Aexps, slong Alen, flint_bitcnt_t Abits,
    const mpoly_ctx_t mctx)
{
    slong Ai, i, j, m;
    slong N = mpoly_words_per_exp(Abits, mctx);
    slong wpf = Abits/FLINT_BITS;
    slong nvars = mctx->nvars;
    slong Aimod, Amodulus = n_sqrt(Alen);

    FLINT_ASSERT(Abits%FLINT_BITS == 0);

    m = 0;

    for (Aimod = 0; Aimod < Amodulus; Aimod++)
    {
        while (m < nvars && used[m])
            m++;

        /* all variables in [0, m) are used */
        if (m >= nvars)
            return;

        if (mctx->rev)
        {
            for (Ai = Aimod; Ai < Alen; Ai += Amodulus)
            {
                for (j = m; j < nvars; j++)
                    for (i = 0; i < wpf && !used[j]; i++)
                        used[j] |= Aexps[N*Ai + wpf*j + i] != 0;
            }
        }
        else
        {
            for (Ai = Aimod; Ai < Alen; Ai += Amodulus)
            {
                for (j = m; j < nvars; j++)
                    for (i = wpf - 1; i >= 0 && !used[j]; i--)
                        used[j] |= Aexps[N*Ai + wpf*(nvars - 1 - j) + i] != 0;
            }
        }
    }
}


/* update used by or'ing it with the used variables in A */
void mpoly_used_vars_or(
    int * used,
    const ulong * Aexps, slong Alen, flint_bitcnt_t Abits,
    const mpoly_ctx_t mctx)
{
    if (Abits <= FLINT_BITS)
        mpoly_used_vars_or_sp(used, Aexps, Alen, Abits, mctx);
    else
        mpoly_used_vars_or_mp(used, Aexps, Alen, Abits, mctx);
}
