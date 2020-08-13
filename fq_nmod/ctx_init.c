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

#include "fq_nmod.h"
#include "nmod_poly.h"

void fq_nmod_ctx_init(fq_nmod_ctx_t ctx, const fmpz_t p, slong d, const char *var)
{
    flint_rand_t state;
    nmod_poly_t poly;

    if (_fq_nmod_ctx_init_conway(ctx, p, d, var))
    {
        ctx->is_conway = 1;
	return;
    } else
        ctx->is_conway = 0;
    
    flint_randinit(state);

    nmod_poly_init2(poly, fmpz_get_ui(p), d + 1);
    nmod_poly_randtest_sparse_irreducible(poly, state, d + 1);

    fq_nmod_ctx_init_modulus(ctx, poly, var);

    nmod_poly_clear(poly);
    flint_randclear(state);    
}
