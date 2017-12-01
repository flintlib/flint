/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"


void nmod_mpoly_make_monic(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                                                    const nmod_mpoly_ctx_t ctx)
{
    if (poly2->length == 0)
    {
        flint_printf("Exception (nmod_mpoly_make_monic). Division by zero.\n");
        flint_abort();
    }

    nmod_mpoly_scalar_mul_ui(poly1, poly2,
                            nmod_inv(poly2->coeffs[0], ctx->ffinfo->mod), ctx);
}

