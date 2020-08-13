/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/* !!! this file DOES need to change with new orderings */

int mpoly_monomial_cmp_general(ulong * Aexp, flint_bitcnt_t Abits,
                       ulong * Bexp, flint_bitcnt_t Bbits, const mpoly_ctx_t mctx)
{
    slong N;

    if (Abits == Bbits)
    {
        /* common case */
        N = mpoly_words_per_exp(Abits, mctx);

        if (!mctx->rev)
        {
            /* ORD_DEGREVLEX and ORD_DEG_LEX */
            return mpoly_monomial_cmp_nomask(Aexp, Bexp, N);
        }
        else
        {
            /* ORD_DEGREVLEX */
            slong i = N - 1;
            if (Abits <= FLINT_BITS)
            {
                /* compare the highest word */
                ulong fpw = FLINT_BITS/Abits;
                ulong himask = (UWORD(1) << (mctx->nvars%fpw*Abits)) - UWORD(1);
                if (Aexp[i] != Bexp[i])
                {
                    if ((Aexp[i]^himask) > (Bexp[i]^himask))
                        return 1;
                    else
                        return -1;
                }
                i--;
            }
            else
            {
                /* compare the degree field with usual comparison */
                ulong wpf = Abits/FLINT_BITS;
                do
                {
                    if (Aexp[i] != Bexp[i])
                    {
                        if (Aexp[i] > Bexp[i])
                            return 1;
                        else
                            return -1;
                    }
                    i--;
                } while (--wpf != 0);
            }
            /* compare the remaining fields with reversed comparisons */
            for (; i >= 0; i--)
            {
                if (Aexp[i] != Bexp[i])
                {
                    if (Aexp[i] < Bexp[i])
                        return 1;
                    else
                        return -1;
                }
            }
            return 0;
        }
    }
    else
    {
        int cmp;
        flint_bitcnt_t newbits;
        ulong * newAexp, * newBexp, * cmpmask;
        TMP_INIT;

        TMP_START;

        if (Abits > Bbits)
        {
            newbits = Abits;
            N = mpoly_words_per_exp(newbits, mctx);
            newAexp = Aexp;
            newBexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
            mpoly_repack_monomials(newBexp, newbits, Bexp, Bbits, 1, mctx);
        }
        else
        {
            FLINT_ASSERT(Abits < Bbits);
            newbits = Bbits;
            N = mpoly_words_per_exp(newbits, mctx);
            newBexp = Bexp;
            newAexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
            mpoly_repack_monomials(newAexp, newbits, Aexp, Abits, 1, mctx);            
        }

        cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
        mpoly_get_cmpmask(cmpmask, N, newbits, mctx);
        cmp = mpoly_monomial_cmp(newAexp, newBexp, N, cmpmask);

        TMP_END;
        return cmp;
    }
}
