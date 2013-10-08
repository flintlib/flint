/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011, 2012 Sebastian Pancratz 
    Copyright (C) 2012 Andres Goens

******************************************************************************/

#include "fq.h"

void fq_randtest(fq_t rop, flint_rand_t state, const fq_ctx_t ctx)
{
    const long d = fq_ctx_degree(ctx);
    long i, sparse;

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

void fq_randtest_not_zero(fq_t rop, flint_rand_t state, const fq_ctx_t ctx)
{
    long i;

    fq_randtest(rop, state, ctx);
    for (i = 0; fq_is_zero(rop) && (i < 10); i++)
        fq_randtest(rop, state, ctx);

    if (fq_is_zero(rop))
        fq_one(rop);
}

