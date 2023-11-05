/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

/*
    sort terms in [left, right) by exponent
    assuming that bits in position >= pos are already sorted
    and assuming exponent vectors fit into one word
    and assuming that all bit positions that need to be sorted are in totalmask
*/
void _fmpz_mod_mpoly_radix_sort1(
    fmpz * Acoeffs,
    ulong * Aexps,
    slong left, slong right,
    flint_bitcnt_t pos,
    ulong cmpmask,
    ulong totalmask)
{
    ulong mask, cmp;
    slong mid, cur;

    while (pos > 0)
    {
        pos--;

        FLINT_ASSERT(left <= right);
        FLINT_ASSERT(pos < FLINT_BITS);

        mask = UWORD(1) << pos;
        cmp = cmpmask & mask;

        /* insertion base case */
        if (right - left < 10)
        {
            slong i, j;

            for (i = left + 1; i < right; i++)
            {
                for (j = i; j > left && mpoly_monomial_gt1(Aexps[j],
                                                   Aexps[j - 1], cmpmask); j--)
                {
                    fmpz_swap(Acoeffs + j, Acoeffs + j - 1);
                    FLINT_SWAP(ulong, Aexps[j], Aexps[j - 1]);
                }
            }

            return;
        }

        /* return if there is no information to sort on this bit */
        if ((totalmask & mask) == 0)
            continue;

        /* find first 'zero' */
        mid = left;
        while (mid < right && (Aexps[mid] & mask) != cmp)
            mid++;

        /* make sure [left,mid)  doesn't match cmpmask in position pos 'one'
                     [mid,right)    does match cmpmask in position pos 'zero' */
        cur = mid;
        while (++cur < right)
        {
            if ((Aexps[cur] & mask) != cmp)
            {
                fmpz_swap(Acoeffs + cur, Acoeffs + mid);
                FLINT_SWAP(ulong, Aexps[cur], Aexps[mid]);
                mid++;
            }
        }

        if (mid - left < right - mid)
        {
            _fmpz_mod_mpoly_radix_sort1(Acoeffs, Aexps, left, mid,
                                                      pos, cmpmask, totalmask);
            left = mid;
        }
        else
        {
            _fmpz_mod_mpoly_radix_sort1(Acoeffs, Aexps, mid, right,
                                                      pos, cmpmask, totalmask);
            right = mid;
        }
    }
}


/*
    sort terms in [left, right) by exponent
    assuming that bits in position >= pos are already sorted
*/
void _fmpz_mod_mpoly_radix_sort(
    fmpz * Acoeffs,
    ulong * Aexps,
    slong left, slong right,
    flint_bitcnt_t pos,
    slong N,
    ulong * cmpmask)
{
    ulong off, bit, mask, cmp;
    slong mid, check;

    while (pos > 0)
    {
        pos--;

        FLINT_ASSERT(left <= right);
        FLINT_ASSERT(pos < N*FLINT_BITS);

        off = pos/FLINT_BITS;
        bit = pos%FLINT_BITS;
        mask = UWORD(1) << bit;
        cmp = cmpmask[off] & mask;

        /* insertion base case */
        if (right - left < 20)
        {
            slong i, j;

            for (i = left + 1; i < right; i++)
            {
                for (j = i; j > left && mpoly_monomial_gt(Aexps + N*j,
                                         Aexps + N*(j - 1), N, cmpmask); j--)
                {
                    fmpz_swap(Acoeffs + j, Acoeffs + j - 1);
                    mpoly_monomial_swap(Aexps + N*j, Aexps + N*(j - 1), N);
                }
            }

            return;
        }

        /* find first 'zero' */
        mid = left;
        while (mid < right && ((Aexps+N*mid)[off] & mask) != cmp)
            mid++;

        /* make sure [left,mid)  doesn't match cmpmask in position pos 'one'
                     [mid,right)    does match cmpmask in position pos 'zero' */
        check = mid;
        while (++check < right)
        {
            if (((Aexps + N*check)[off] & mask) != cmp)
            {
                fmpz_swap(Acoeffs + check, Acoeffs + mid);
                mpoly_monomial_swap(Aexps + N*check, Aexps + N*mid, N);
                mid++;
            }
        }

        FLINT_ASSERT(left <= mid && mid <= right);

        if (mid - left < right - mid)
        {
            _fmpz_mod_mpoly_radix_sort(Acoeffs, Aexps, left, mid,
                                                              pos, N, cmpmask);
            left = mid;
        }
        else
        {
            _fmpz_mod_mpoly_radix_sort(Acoeffs, Aexps, mid, right,
                                                              pos, N, cmpmask);
            right = mid;
        }
    }
}


/*
    sort the terms in A by exponent
    assuming that the exponents are valid (other than being in order)
*/
void fmpz_mod_mpoly_sort_terms(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, N;
    flint_bitcnt_t pos;
    fmpz * Acoeffs = A->coeffs;
    ulong * Aexps = A->exps;
    ulong himask, * ptempexp;
    TMP_INIT;

    TMP_START;
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    ptempexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(ptempexp, N, A->bits, ctx->minfo);

    himask = 0;
    for (i = 0; i < A->length; i++)
        himask |= (Aexps + N*i)[N - 1];

    pos = FLINT_BIT_COUNT(himask);
    if (N == 1)
        _fmpz_mod_mpoly_radix_sort1(Acoeffs, Aexps, 0, A->length,
                                                     pos, ptempexp[0], himask);
    else
        _fmpz_mod_mpoly_radix_sort(Acoeffs, Aexps, 0, A->length,
                                        (N - 1)*FLINT_BITS + pos, N, ptempexp);

    TMP_END;
}
