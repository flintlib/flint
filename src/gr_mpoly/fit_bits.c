/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "gr_mpoly.h"

void
gr_mpoly_fit_bits(gr_mpoly_t A, flint_bitcnt_t bits, const mpoly_ctx_t mctx, gr_ctx_t FLINT_UNUSED(cctx))
{
    if (A->bits < bits)
    {
        if (A->exps_alloc != 0)
        {
            slong N = mpoly_words_per_exp(bits, mctx);
            ulong * t = (ulong *) flint_malloc(N*A->exps_alloc*sizeof(ulong));
            mpoly_repack_monomials(t, bits, A->exps, A->bits, A->length, mctx);
            flint_free(A->exps);
            A->exps = t;
            A->exps_alloc = N*A->exps_alloc;
        }

        A->bits = bits;
    }
}
