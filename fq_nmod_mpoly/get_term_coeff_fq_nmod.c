/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_get_term_coeff_fq_nmod(fq_nmod_t c, const fq_nmod_mpoly_t A,
                                        slong i, const fq_nmod_mpoly_ctx_t ctx)
{
    if ((ulong) i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR, "index out of range in fmpz_mpoly_get_term_coeff_fmpz");
    }

    fq_nmod_set(c, A->coeffs + i, ctx->fqctx);
}
