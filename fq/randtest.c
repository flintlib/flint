/*
    Copyright (C) 2011, 2012 Sebastian Pancratz 
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq.h"

void
fq_randtest(fq_t rop, flint_rand_t state, const fq_ctx_t ctx)
{
    const slong d = fq_ctx_degree(ctx);
    slong i, sparse;

    fmpz_poly_fit_length(rop, d);

    if (n_randint(state, 2))
    {
        for (i = 0; i < d; i++)
            fmpz_randm(rop->coeffs + i, state, fq_ctx_prime(ctx));
    }
    else
    {
        sparse = 1 + n_randint(state, FLINT_MAX(2, d));

        for (i = 0; i < d; i++)
            if (n_randint(state, sparse))
                fmpz_zero(rop->coeffs + i);
            else
                fmpz_randm(rop->coeffs + i, state, fq_ctx_prime(ctx));
    }

    _fmpz_poly_set_length(rop, d);
    _fmpz_poly_normalise(rop);
}

void
fq_randtest_dense(fq_t rop, flint_rand_t state, const fq_ctx_t ctx)
{
    const slong d = fq_ctx_degree(ctx);
    slong i;

    fmpz_poly_fit_length(rop, d);

    for (i = 0; i < d - 1; i++)
    {
        fmpz_randm(rop->coeffs + i, state, fq_ctx_prime(ctx));
    }

    fmpz_one(rop->coeffs + d - 1);

    _fmpz_poly_set_length(rop, d);
    _fmpz_poly_normalise(rop);
}

void
fq_randtest_not_zero(fq_t rop, flint_rand_t state, const fq_ctx_t ctx)
{
    slong i;

    fq_randtest(rop, state, ctx);
    for (i = 0; fq_is_zero(rop, ctx) && (i < 10); i++)
        fq_randtest(rop, state, ctx);

    if (fq_is_zero(rop, ctx))
        fq_one(rop, ctx);
}
