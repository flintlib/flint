/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"

void
fq_nmod_ctx_randtest(fq_nmod_ctx_t ctx, flint_rand_t state)
{
    nmod_poly_t modulus;
    mp_limb_t x;
    fmpz_t p;
    slong d;

    fmpz_init(p);
    fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 6), 1));
    d = n_randint(state, 10) + 1;
    fq_nmod_ctx_init_conway(ctx, p, d, "a");
    fmpz_clear(p);

    /* Test non-monic modulus */
    if (n_randint(state, 2))
    {
        nmod_poly_init(modulus, ctx->mod.n);
        nmod_poly_set(modulus, ctx->modulus);
        x = n_randint(state, ctx->mod.n - 1) + 1;
        nmod_poly_scalar_mul_nmod(modulus, modulus, x);
        fq_nmod_ctx_clear(ctx);
        fq_nmod_ctx_init_modulus(ctx, modulus, "a");
        nmod_poly_clear(modulus);
    }
}
