/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "fq_zech_mpoly.h"

void fq_zech_mpoly_fit_bits(fq_zech_mpoly_t A, slong bits, const fq_zech_mpoly_ctx_t ctx)
{
    if (A->bits < (ulong) bits)
    {
        if (A->alloc != 0)
        {
            slong N = mpoly_words_per_exp(bits, ctx->minfo);
            ulong * t = (ulong *) flint_malloc(N*A->alloc*sizeof(ulong));
            mpoly_repack_monomials(t, bits, A->exps, A->bits, A->length,
                                                                   ctx->minfo);
            flint_free(A->exps);
            A->exps = t;
        }

        A->bits = bits;
    }
}
