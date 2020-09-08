/*
    Copyright (C) 2018 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "fq.h"

void
fmpz_mod_poly_randtest_monic_primitive(fmpz_mod_poly_t f,
                       flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx)
{
    fq_ctx_t fqctx;
    fq_t X;
    int primitive = 0;
    
    while (!primitive)
    {
        fmpz_mod_poly_randtest_monic_irreducible(f, state, len, ctx);
        fq_ctx_init_modulus(fqctx, f, ctx, "X");
        fq_init(X, fqctx);
        fq_gen(X, fqctx);
        primitive = fq_is_primitive(X, fqctx);
        fq_clear(X, fqctx);
        fq_ctx_clear(fqctx);
    }
}
