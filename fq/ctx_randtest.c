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

#include "fq.h"

void
fq_ctx_randtest(fq_ctx_t ctx, flint_rand_t state)
{
    fmpz_mod_poly_t modulus;
    fmpz_t p, x;
    slong d;

    fmpz_init(p);
    fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 6), 1));
    d = n_randint(state, 10) + 1;
    fq_ctx_init_conway(ctx, p, d, "a");
    fmpz_clear(p);

    /* Test non-monic modulus */
    if (n_randint(state, 2))
    {
        fmpz_mod_ctx_t ctxp;
        fmpz_mod_ctx_init(ctxp, p);
        fmpz_init_set(x, p);
        fmpz_sub_ui(x, x, 1);
        fmpz_mod_poly_init(modulus, ctxp);
        fmpz_mod_poly_set(modulus, ctx->modulus, ctxp);
        fmpz_randm(x, state, x);
        fmpz_add_ui(x, x, 1);
        fmpz_mod_poly_scalar_mul_fmpz(modulus, modulus, x, ctxp);
        fq_ctx_clear(ctx);
        fq_ctx_init_modulus(ctx, modulus, ctxp, "a");
        fmpz_mod_poly_clear(modulus, ctxp);
        fmpz_mod_ctx_clear(ctxp);
        fmpz_clear(x);
    }
}
