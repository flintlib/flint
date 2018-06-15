/*
    Copyright (C) 2018 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "fq.h"

void
fmpz_mod_poly_randtest_monic_primitive(fmpz_mod_poly_t f,
                                       flint_rand_t state, slong len)
{
    fq_ctx_t ctx;
    fq_t X;
    int primitive = 0;
    
    while (!primitive)
    {
        fmpz_mod_poly_randtest_monic_irreducible(f, state, len);
        fq_ctx_init_modulus(ctx, f, "X");
        fq_init(X, ctx);
        fq_gen(X, ctx);
        primitive = fq_is_primitive(X, ctx);
        fq_clear(X, ctx);
        fq_ctx_clear(ctx);
    }
}
