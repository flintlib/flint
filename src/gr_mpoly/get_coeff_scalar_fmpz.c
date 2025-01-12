/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "mpoly.h"
#include "gr_mpoly.h"

/* todo: should have an
   mpoly_monomial_index_pfmpz function */

int gr_mpoly_get_coeff_scalar_fmpz(gr_ptr c, const gr_mpoly_t A,
                            const fmpz * exp, gr_mpoly_ctx_t ctx)
{
    int status;
    slong index;
    slong nvars = GR_MPOLY_NVARS(ctx);
    fmpz ** exp_ptr;
    slong i;
    TMP_INIT;

    for (i = 0; i < nvars; i++)
    {
        if (fmpz_sgn(exp + i) < 0)
            return GR_DOMAIN;
    }

    TMP_START;

    exp_ptr = (fmpz **) TMP_ALLOC(sizeof(fmpz *) * nvars);
    for (i = 0; i < nvars; i++)
        exp_ptr[i] = (fmpz *) (exp + i);

    index = mpoly_monomial_index_pfmpz(A->exps, A->bits, A->length, exp_ptr, GR_MPOLY_MCTX(ctx));
    if (index < 0)
    {
        status = gr_zero(c, GR_MPOLY_CCTX(ctx));
    }
    else
    {
        FLINT_ASSERT(index < A->length);
        status = gr_set(c, GR_ENTRY(A->coeffs, index, GR_MPOLY_CCTX(ctx)->sizeof_elem), GR_MPOLY_CCTX(ctx));
    }

    TMP_END;

    return status;
}
