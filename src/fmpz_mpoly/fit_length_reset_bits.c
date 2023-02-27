/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

/*
    Ensure that A has space for "len" terms with bits "bits" and set A->bits
    The value of A is destroyed.

    Since the coefficient alloc and the exponent alloc are coupled, we must
    have the assumption that A is allocated wrt ctx upon entry.

    TODO: decouple the two allocations as in nmod_mpoly/fq_nmod_mpoly
*/
void fmpz_mpoly_fit_length_reset_bits(
    fmpz_mpoly_t A,
    slong len,
    flint_bitcnt_t bits,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong oldN = mpoly_words_per_exp(A->bits, ctx->minfo);
    slong newN = mpoly_words_per_exp(bits, ctx->minfo);

    if (len > A->alloc)
    {
        len = FLINT_MAX(len, 2*A->alloc);

        A->exps = (ulong *) flint_realloc(A->exps, newN*len*sizeof(ulong));
        A->coeffs = (fmpz *) flint_realloc(A->coeffs, len*sizeof(fmpz));

        for (i = A->alloc; i < len; i++)
            fmpz_init(A->coeffs + i);

        A->alloc = len;
    }
    else if (newN > oldN && A->alloc > 0)
    {
        A->exps = (ulong *) flint_realloc(A->exps, newN*A->alloc*sizeof(ulong));
    }

    A->bits = bits;
}
