/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"


void fmpz_mod_mpoly_fit_length_fit_bits(
    fmpz_mod_mpoly_t A,
    slong len,
    flint_bitcnt_t bits,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, N = mpoly_words_per_exp(A->bits, ctx->minfo);

    if (len > A->coeffs_alloc)
    {
        slong old_alloc = A->coeffs_alloc;
        slong new_alloc = FLINT_MAX(len, 2*old_alloc);
        A->coeffs_alloc = new_alloc;
        A->coeffs = (fmpz*) flint_realloc(A->coeffs, new_alloc*sizeof(fmpz));
        for (i = old_alloc; i < new_alloc; i++)
            fmpz_init(A->coeffs + i);
    }

    if (bits > A->bits)
    {
        slong newN = mpoly_words_per_exp(bits, ctx->minfo);
        slong new_exps_alloc = newN*len;
        ulong * t;

        if (len < 1)
        {
            A->bits = bits;
            return;
        }

        t = (ulong *) flint_malloc(new_exps_alloc*sizeof(ulong));

        if (A->length > 0)
            mpoly_repack_monomials(t, bits, A->exps, A->bits, A->length, ctx->minfo);

        if (A->exps_alloc > 0)
            flint_free(A->exps);

        A->exps = t;
        A->exps_alloc = new_exps_alloc;
        A->bits = bits;
    }
    else
    {
        if (N*len > A->exps_alloc)
        {
            A->exps_alloc = FLINT_MAX(N*len, 2*A->exps_alloc);
            A->exps = flint_realloc(A->exps, A->exps_alloc*sizeof(ulong));
        }
    }
}
