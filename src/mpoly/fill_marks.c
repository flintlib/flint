/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/*
    Fill in the indices that partition A into univariate chunks.
    D will be allocated to at least one more than its length so that
        D->coeffs[D->length] = Alen
*/
void mpoly1_fill_marks(
    ulong ** Dcoeffs,
    slong * Dlength,
    slong * Dalloc,
    const ulong * Aexps,
    slong Alen,
    flint_bitcnt_t Abits,
    const mpoly_ctx_t mctx)
{
    slong off0, shift0;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - Abits);
    slong i, N = mpoly_words_per_exp_sp(Abits, mctx);
    ulong e0;

    FLINT_ASSERT(mctx->ord == ORD_LEX);
    FLINT_ASSERT(mctx->nvars >= 1);

    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, Abits, mctx);

    *Dlength = 0;
    i = 0;
    while (i < Alen)
    {
        if (*Dalloc < *Dlength + 1)
        {
            *Dalloc = FLINT_MAX(*Dlength + 1, *Dalloc + *Dalloc/2);
            *Dcoeffs = FLINT_ARRAY_REALLOC(*Dcoeffs, *Dalloc, ulong);
        }
        (*Dcoeffs)[*Dlength] = i;
        *Dlength = *Dlength + 1;

        e0 = (Aexps[N*i + off0] >> shift0) & mask;
        while (1)
        {
            i++;
            if (i >= Alen)
                break;
            if (((Aexps[N*i + off0] >> shift0) & mask) != e0)
                break;
        }
    }

    if (*Dalloc < *Dlength + 1)
    {
        *Dalloc = FLINT_MAX(*Dlength + 1, *Dalloc + *Dalloc/2);
        *Dcoeffs = FLINT_ARRAY_REALLOC(*Dcoeffs, *Dalloc, ulong);
    }
    (*Dcoeffs)[*Dlength] = Alen;
}

/*
    Fill in the indices that partition A into bivariate chunks.
    D will be allocated to at least one more than its length so that
        D->coeffs[D->length] = Alen
*/
void mpoly2_fill_marks(
    ulong ** Dcoeffs,
    slong * Dlength,
    slong * Dalloc,
    const ulong * Aexps,
    slong Alen,
    flint_bitcnt_t Abits,
    const mpoly_ctx_t mctx)
{
    slong off0, off1, shift0, shift1;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - Abits);
    slong i, N = mpoly_words_per_exp_sp(Abits, mctx);
    ulong e0, e1;

    FLINT_ASSERT(mctx->ord == ORD_LEX);
    FLINT_ASSERT(mctx->nvars >= 2);

    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, Abits, mctx);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, Abits, mctx);

    *Dlength = 0;
    i = 0;
    while (i < Alen)
    {
        if (*Dalloc < *Dlength + 1)
        {
            *Dalloc = FLINT_MAX(*Dlength + 1, *Dalloc + *Dalloc/2);
            *Dcoeffs = FLINT_ARRAY_REALLOC(*Dcoeffs, *Dalloc, ulong);
        }
        (*Dcoeffs)[*Dlength] = i;
        *Dlength = *Dlength + 1;

        e0 = (Aexps[N*i + off0] >> shift0) & mask;
        e1 = (Aexps[N*i + off1] >> shift1) & mask;
        while (1)
        {
            i++;
            if (i >= Alen)
                break;
            if (((Aexps[N*i + off0] >> shift0) & mask) != e0)
                break;
            if (((Aexps[N*i + off1] >> shift1) & mask) != e1)
                break;
        }
    }

    if (*Dalloc < *Dlength + 1)
    {
        *Dalloc = FLINT_MAX(*Dlength + 1, *Dalloc + *Dalloc/2);
        *Dcoeffs = FLINT_ARRAY_REALLOC(*Dcoeffs, *Dalloc, ulong);
    }
    (*Dcoeffs)[*Dlength] = Alen;
}
