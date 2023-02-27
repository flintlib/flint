/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq.h"


void
fq_rand(fq_t rop, flint_rand_t state, const fq_ctx_t ctx)
{
    const slong d = fq_ctx_degree(ctx);
    slong i;

    fmpz_poly_fit_length(rop, d);

    for (i = 0; i < d; i++)
    {
        fmpz_randm(rop->coeffs + i, state, fq_ctx_prime(ctx));
    }

    _fmpz_poly_set_length(rop, d);
    _fmpz_poly_normalise(rop);
}

void fq_rand_not_zero(fq_t rop, flint_rand_t state, const fq_ctx_t ctx)
{
    int tries = 3;

    do {
        fq_rand(rop, state, ctx);
        if (!fq_is_zero(rop, ctx))
            return;
    } while (--tries >= 0);

    fq_one(rop, ctx);
}

