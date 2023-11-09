/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech.h"
#include "fq_zech_mpoly.h"

/*
    sort terms in [left, right) by exponent
    assuming that bits in position >= pos are already sorted
    and assuming exponent vectors fit into one word
    and assuming that all bit positions that need to be sorted are in totalmask
*/
void _fq_zech_mpoly_radix_sort1(fq_zech_mpoly_t A, slong left, slong right,
                            flint_bitcnt_t pos, ulong cmpmask, ulong totalmask)
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
        if (right - left < 20)
        {
            slong i, j;

            for (i = left + 1; i < right; i++)
            {
                for (j = i; j > left && mpoly_monomial_gt1(A->exps[j],
                                                 A->exps[j - 1], cmpmask); j--)
                {
                    fq_zech_swap(A->coeffs + j, A->coeffs + j - 1, NULL);
                    FLINT_SWAP(ulong, A->exps[j], A->exps[j - 1]);
                }
            }

            return;
        }

        /* return if there is no information to sort on this bit */
        if ((totalmask & mask) == 0)
            continue;

        /* find first 'zero' */
        mid = left;
        while (mid < right && (A->exps[mid] & mask) != cmp)
            mid++;

        /* make sure [left,mid)  doesn't match cmpmask in position pos 'one'
                     [mid,right)    does match cmpmask in position pos 'zero' */
        cur = mid;
        while (++cur < right)
        {
            if ((A->exps[cur] & mask) != cmp)
            {
                fq_zech_swap(A->coeffs + cur, A->coeffs + mid, NULL);
                FLINT_SWAP(ulong, A->exps[cur], A->exps[mid]);
                mid++;
            }
        }

        if (mid - left < right - mid)
        {
            _fq_zech_mpoly_radix_sort1(A, left, mid, pos, cmpmask, totalmask);
            left = mid;
        }
        else
        {
            _fq_zech_mpoly_radix_sort1(A, mid, right, pos, cmpmask, totalmask);
            right = mid;
        }
    }
}


/*
    sort terms in [left, right) by exponent
    assuming that bits in position > pos are already sorted
*/
void _fq_zech_mpoly_radix_sort(fq_zech_mpoly_t A, slong left, slong right,
                                  flint_bitcnt_t pos, slong N, ulong * cmpmask)
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
        if (right - left < 10)
        {
            slong i, j;

            for (i = left + 1; i < right; i++)
            {
                for (j = i; j > left && mpoly_monomial_gt(A->exps + N*j,
                                         A->exps + N*(j - 1), N, cmpmask); j--)
                {
                    fq_zech_swap(A->coeffs + j, A->coeffs + j - 1, NULL);
                    mpoly_monomial_swap(A->exps + N*j, A->exps + N*(j - 1), N);
                }
            }

            return;
        }

        /* find first 'zero' */
        mid = left;
        while (mid < right && ((A->exps+N*mid)[off] & mask) != cmp)
            mid++;

        /* make sure [left,mid)  doesn't match cmpmask in position pos 'one'
                     [mid,right)    does match cmpmask in position pos 'zero' */
        check = mid;
        while (++check < right)
        {
            if (((A->exps + N*check)[off] & mask) != cmp)
            {
                fq_zech_swap(A->coeffs + check, A->coeffs + mid, NULL);
                mpoly_monomial_swap(A->exps + N*check, A->exps + N*mid, N);
                mid++;
            }
        }

        FLINT_ASSERT(left <= mid && mid <= right);

        if (mid - left < right - mid)
        {
            _fq_zech_mpoly_radix_sort(A, left, mid, pos, N, cmpmask);
            left = mid;
        }
        else
        {
            _fq_zech_mpoly_radix_sort(A, mid, right, pos, N, cmpmask);
            right = mid;
        }
    }
}


/*
    sort the terms in A by exponent
    assuming that the exponents are valid (other than being in order)
*/
void fq_zech_mpoly_sort_terms(fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
{
    slong i, N;
    flint_bitcnt_t pos;
    ulong himask, * ptempexp;
    TMP_INIT;

    TMP_START;
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    ptempexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(ptempexp, N, A->bits, ctx->minfo);

    himask = 0;
    for (i = 0; i < A->length; i++)
        himask |= (A->exps + N*i)[N - 1];

    pos = FLINT_BIT_COUNT(himask);

    if (N == 1)
        _fq_zech_mpoly_radix_sort1(A, 0, A->length, pos, ptempexp[0], himask);
    else
        _fq_zech_mpoly_radix_sort(A, 0, A->length,
                                        (N - 1)*FLINT_BITS + pos, N, ptempexp);

    TMP_END;
}
