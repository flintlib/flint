/*
    Copyright (C) 2017 Luca De Feo
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2023 Albin Ahlbäck

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

void fq_nmod_set_fmpz(fq_nmod_t rop, const fmpz_t x, const fq_nmod_ctx_t ctx)
{
    fmpz_t rx;
    fmpz_init(rx);
    fmpz_mod(rx, x, fq_nmod_ctx_prime(ctx));
    nmod_poly_zero(rop);
    nmod_poly_set_coeff_ui(rop, 0, fmpz_get_ui(rx));
    fmpz_clear(rx);
}

void fq_nmod_set_nmod_poly(fq_nmod_t a, const nmod_poly_t b,
                                                       const fq_nmod_ctx_t ctx)
{
    FLINT_ASSERT(a->mod.n == b->mod.n);
    FLINT_ASSERT(a->mod.n == ctx->modulus->mod.n);

    if (b->length <= 2*(ctx->modulus->length - 1))
    {
        nmod_poly_set(a, b);
        fq_nmod_reduce(a, ctx);
    }
    else
    {
        nmod_poly_rem(a, b, ctx->modulus);
    }
}

void fq_nmod_set_nmod_mat(fq_nmod_t a,
                          const nmod_mat_t col,
                          const fq_nmod_ctx_t ctx)
{
    slong i, n = fq_nmod_ctx_degree(ctx);
    nmod_poly_fit_length(a, n);
    a->length = n;
    for (i = 0; i < n; i++)
        a->coeffs[i] = nmod_mat_entry(col, i, 0);
    _nmod_poly_normalise(a);
}
