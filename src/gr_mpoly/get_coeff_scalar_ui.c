/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mpoly.h"

/* todo: should have an
   mpoly_monomial_index_pfmpz function */

int gr_mpoly_get_coeff_scalar_ui(gr_ptr c, const gr_mpoly_t A,
                            const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
{

    slong index;
    index = mpoly_monomial_index_ui(A->exps, A->bits, A->length, exp, mctx);

    if (index < 0)
    {
        return gr_zero(c, cctx);
    }
    else
    {
        FLINT_ASSERT(index < A->length);
        return gr_set(c, GR_ENTRY(A->coeffs, index, cctx->sizeof_elem), cctx);
    }
}
