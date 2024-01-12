/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "fq_nmod.h"

void
fq_nmod_ctx_randtest(fq_nmod_ctx_t ctx, flint_rand_t state)
{
    struct _prime_degree_struct pd;

    pd = _nmod_poly_conway_rand(state);
    fq_nmod_ctx_init_conway_ui(ctx, pd.prime, pd.degree, "a");

    /* Test non-monic modulus */
    if (n_randint(state, 2))
    {
        nmod_poly_t modulus;
        mp_limb_t x;

        nmod_poly_init(modulus, ctx->mod.n);
        nmod_poly_set(modulus, ctx->modulus);
        x = n_randint(state, ctx->mod.n - 1) + 1;
        nmod_poly_scalar_mul_nmod(modulus, modulus, x);
        fq_nmod_ctx_clear(ctx);
        fq_nmod_ctx_init_modulus(ctx, modulus, "a");
        nmod_poly_clear(modulus);
    }
}
