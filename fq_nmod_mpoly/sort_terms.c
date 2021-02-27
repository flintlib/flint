/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"


/*
    sort terms in [left, right) by exponent
    assuming that bits in position > pos are already sorted
    and assuming exponent vectors fit into one word
    and assuming that all bit positions that need to be sorted are in totalmask
*/
void _fq_nmod_mpoly_radix_sort1(
    fq_nmod_mpoly_t A,
    slong left,
    slong right,
    flint_bitcnt_t pos,
    ulong cmpmask,
    ulong totalmask,
    slong d)
{
    ulong mask = UWORD(1) << pos;
    ulong cmp = cmpmask & mask;
    slong mid, cur;

    FLINT_ASSERT(left <= right);
    FLINT_ASSERT(pos < FLINT_BITS);

    /* do nothing on lists of 0 or 1 elements */
    if (left + 1 >= right)
    {
        return;
    }

    /* return if there is no information to sort on this bit */
    if ((totalmask & mask) == WORD(0))
    {
        if (pos < 1)
            return;

        pos--;
        _fq_nmod_mpoly_radix_sort1(A, left, right, pos, cmpmask, totalmask, d);
        return;
    }

    /* find first 'zero' */
    mid = left;
    while (mid < right && ((A->exps + 1*mid)[0] & mask) != cmp)
    {
        mid++;
    }

    /* make sure [left,mid)  doesn't match cmpmask in position pos 'one'
                 [mid,right)    does match cmpmask in position pos 'zero' */
    cur = mid;
    while (++cur < right)
    {
        if (((A->exps + 1*cur)[0] & mask) != cmp)
        {
            _n_fq_swap(A->coeffs + d*cur, A->coeffs + d*mid, d);
            mpoly_monomial_swap(A->exps + 1*cur, A->exps + 1*mid, 1);
            mid++;
        }
    }

    if (pos < 1)
        return;

    pos--;
    _fq_nmod_mpoly_radix_sort1(A, left,  mid, pos, cmpmask, totalmask, d);
    _fq_nmod_mpoly_radix_sort1(A, mid, right, pos, cmpmask, totalmask, d);
}


/*
    sort terms in [left, right) by exponent
    assuming that bits in position > pos are already sorted

    TODO: Stack depth is proportional to N*FLINT_BITS
            Might turn into iterative version
            Low priority
*/
void _fq_nmod_mpoly_radix_sort(
    fq_nmod_mpoly_t A,
    slong left,
    slong right,
    flint_bitcnt_t pos,
    slong N,
    ulong * cmpmask,
    slong d)
{
    ulong off = pos/FLINT_BITS;
    ulong bit = pos%FLINT_BITS;
    ulong mask = UWORD(1) << bit;
    ulong cmp = cmpmask[off] & mask;
    slong mid, check;

    FLINT_ASSERT(left <= right);
    FLINT_ASSERT(pos < N*FLINT_BITS);

    /* do nothing on lists of 0 or 1 elements */
    if (left + 1 >= right)
        return;

    /* find first 'zero' */
    mid = left;
    while (mid < right && ((A->exps+N*mid)[off] & mask) != cmp)
    {
        mid++;
    }

    /* make sure [left,mid)  doesn't match cmpmask in position pos 'one'
                 [mid,right)    does match cmpmask in position pos 'zero' */
    check = mid;
    while (++check < right)
    {
        if (((A->exps + N*check)[off] & mask) != cmp)
        {
            _n_fq_swap(A->coeffs + d*check, A->coeffs + d*mid, d);
            mpoly_monomial_swap(A->exps + N*check, A->exps + N*mid, N);
            mid++;
        }
    }

    if (pos < 1)
        return;

    pos--;
    _fq_nmod_mpoly_radix_sort(A, left,  mid, pos, N, cmpmask, d);
    _fq_nmod_mpoly_radix_sort(A, mid, right, pos, N, cmpmask, d);
}


/*
    sort the terms in A by exponent
    assuming that the exponents are valid (other than being in order)
*/
void fq_nmod_mpoly_sort_terms(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    slong d;
    slong i, msb, N;
    ulong himask, * ptempexp;
    TMP_INIT;

    TMP_START;
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    ptempexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(ptempexp, N, A->bits, ctx->minfo);

    himask = 0;
    for (i = 0; i < A->length; i++)
    {
        himask |= (A->exps + N*i)[N - 1];
    }

    if (himask != 0)
    {
        count_leading_zeros(msb, himask);
        msb = (FLINT_BITS - 1)^msb;
    }
    else
    {
        msb = -WORD(1);
    }

    d = fq_nmod_ctx_degree(ctx->fqctx);

    if (N == 1)
    {
        if (msb >= 0)
        {
            _fq_nmod_mpoly_radix_sort1(A, 0, A->length, msb, ptempexp[0], himask, d);
        }
    }
    else
    {
        _fq_nmod_mpoly_radix_sort(A, 0, A->length, (N - 1)*FLINT_BITS + msb, N, ptempexp, d);
    }

    TMP_END;
}
