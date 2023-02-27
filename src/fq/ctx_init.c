/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
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
    fmpz_mod_ctx_t ctxp;

    if (_fq_ctx_init_conway(ctx, p, d, var))
    {
        ctx->is_conway = 1;
	    return;
    }
    else
    {
	    ctx->is_conway = 0;
    }

    flint_randinit(state);

    fmpz_mod_ctx_init(ctxp, p);
    fmpz_mod_poly_init2(poly, d + 1, ctxp);
    fmpz_mod_poly_randtest_sparse_irreducible(poly, state, d + 1, ctxp);

    fq_ctx_init_modulus(ctx, poly, ctxp, var);

    fmpz_mod_poly_clear(poly, ctxp);
    fmpz_mod_ctx_clear(ctxp);
    flint_randclear(state);
}
