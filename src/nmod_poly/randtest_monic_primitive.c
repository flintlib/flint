/*
    Copyright (C) 2018 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "fq_nmod.h"

void
nmod_poly_randtest_monic_primitive(nmod_poly_t poly, flint_rand_t state, slong len)
{
    fq_nmod_ctx_t ctx;
    fq_nmod_t X;
    int primitive = 0;
    
    while (!primitive)
    {
        nmod_poly_randtest_monic_irreducible(poly, state, len);
        fq_nmod_ctx_init_modulus(ctx, poly, "X");
        fq_nmod_init(X, ctx);
        fq_nmod_gen(X, ctx);
        primitive = fq_nmod_is_primitive(X, ctx);
        fq_nmod_clear(X, ctx);
        fq_nmod_ctx_clear(ctx);
    }
}
