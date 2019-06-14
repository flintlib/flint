/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"
#include "ulong_extras.h"

/*
    Scan A and fill in the min and max exponents of each variable along
    with the count of terms attached to each.
*/
void mpoly_gcd_info_limits(ulong * Amax_exp, ulong * Amin_exp,
                       slong * Amax_exp_count, slong * Amin_exp_count,
                       const ulong * Aexps, flint_bitcnt_t Abits, slong Alength,
                                                        const mpoly_ctx_t mctx)
{
    ulong * exps;
    slong i, j, N;
    slong nvars = mctx->nvars;
    TMP_INIT;

    FLINT_ASSERT(Alength > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);

    TMP_START;

    exps = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));

    N = mpoly_words_per_exp(Abits, mctx);

    i = 0;
    mpoly_get_monomial_ui(exps, Aexps + N*i, Abits, mctx);
    for (j = 0; j < nvars; j++)
    {
        Amin_exp[j] = exps[j];
        Amax_exp[j] = exps[j];
        Amin_exp_count[j] = 1;
        Amax_exp_count[j] = 1;
    }
    for (i = 1; i < Alength; i++)
    {
        mpoly_get_monomial_ui(exps, Aexps + N*i, Abits, mctx);

        for (j = 0; j < nvars; j++)
        {
            if (Amin_exp[j] > exps[j])
            {
                Amin_exp[j] = exps[j];
                Amin_exp_count[j] = 1;            
            }
            else if (Amin_exp[j] == exps[j])
            {
                Amin_exp_count[j] += 1;
            }

            if (Amax_exp[j] < exps[j])
            {
                Amax_exp[j] = exps[j];
                Amax_exp_count[j] = 1;            
            }
            else if (Amax_exp[j] == exps[j])
            {
                Amax_exp_count[j] += 1;
            }
        }
    }

    TMP_END;
}


/*
    For each variable v, let SA[v] be the set of exponents of variable v in A.
    Ditto for SB[v]. The function computes
        strides[v] = GCD(SA[v] - min(SA[v]), SB[v] - min(SB[v]))
    It is assumed that {A|B}{max|min}_exp are correct.
*/
void mpoly_gcd_info_stride(ulong * strides,
          const ulong * Aexps, flint_bitcnt_t Abits, slong Alength,
                             const ulong * Amax_exp, const ulong * Amin_exp,
          const ulong * Bexps, flint_bitcnt_t Bbits, slong Blength,
                             const ulong * Bmax_exp, const ulong * Bmin_exp,
                                                        const mpoly_ctx_t mctx)
{
    slong i, j, NA, NB;
    slong nvars = mctx->nvars;
    ulong mask;
    ulong * exps;
    TMP_INIT;

    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(Bbits <= FLINT_BITS);

    for (j = 0; j < nvars; j++)
    {
        strides[j] = n_gcd(Amax_exp[j] - Amin_exp[j],
                           Bmax_exp[j] - Bmin_exp[j]);
    }

    TMP_START;
    exps = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));

    NA = mpoly_words_per_exp(Abits, mctx);

    for (i = 0; i < Alength; i++)
    {
        mpoly_get_monomial_ui(exps, Aexps + NA*i, Abits, mctx);
        mask = 0;
        for (j = 0; j < nvars; j++)
        {
            strides[j] = n_gcd(strides[j], exps[j] - Amin_exp[j]);
            mask |= strides[j];
        }
        if (mask < UWORD(2))
        {
            goto cleanup;
        }
    }

    NB = mpoly_words_per_exp(Bbits, mctx);

    for (i = 0; i < Blength; i++)
    {
        mpoly_get_monomial_ui(exps, Bexps + NB*i, Bbits, mctx);
        mask = 0;
        for (j = 0; j < nvars; j++)
        {
            strides[j] = n_gcd(strides[j], exps[j] - Bmin_exp[j]);
            mask |= strides[j];
        }
        if (mask < UWORD(2))
        {
            goto cleanup;
        }
    }

cleanup:

    TMP_END;

    return;
}
