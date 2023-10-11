/*
    Copyright (C) 2017 Luca De Feo
    Copyright (C) 2020, 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "nmod_mat.h"
#include "fmpz.h"
#include "fq_nmod.h"

int fq_nmod_get_fmpz(fmpz_t a, const fq_nmod_t b, const fq_nmod_ctx_t ctx)
{
    if (b->length > 1)
        return 0;

    if (b->length == 1)
        fmpz_set_ui(a, b->coeffs[0]);
    else
        fmpz_zero(a);

    return 1;
}

void fq_nmod_get_nmod_poly(nmod_poly_t a, const fq_nmod_t b,
                                                       const fq_nmod_ctx_t ctx)
{
    FLINT_ASSERT(b->mod.n == ctx->modulus->mod.n);
    a->mod = ctx->modulus->mod;
    nmod_poly_set(a, b);
}

void fq_nmod_get_nmod_mat(nmod_mat_t col,
                          const fq_nmod_t a,
                          const fq_nmod_ctx_t ctx)
{
    slong i, n = fq_nmod_ctx_degree(ctx);
    for (i = 0; i < a->length; i++)
        nmod_mat_set_entry(col, i, 0, a->coeffs[i]);
    for ( ; i < n; i++)
        nmod_mat_entry(col, i, 0) = 0;
}
