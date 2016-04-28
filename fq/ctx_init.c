/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <string.h>

#include "fq.h"
#include "fq_poly.h"

void
fq_ctx_init(fq_ctx_t ctx, const fmpz_t p, slong d, const char *var)
{
    flint_rand_t state;
    fmpz_mod_poly_t poly;

    if (_fq_ctx_init_conway(ctx, p, d, var))
    {
        return;
    }

    flint_randinit(state);

    fmpz_mod_poly_init2(poly, p, d + 1);
    fmpz_mod_poly_randtest_sparse_irreducible(poly, state, d + 1);

    fq_ctx_init_modulus(ctx, poly, var);

    fmpz_mod_poly_clear(poly);
    flint_randclear(state);
}
