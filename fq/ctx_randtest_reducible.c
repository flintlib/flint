/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq.h"

void
fq_ctx_randtest_reducible(fq_ctx_t ctx, flint_rand_t state)
{
    fmpz_mod_ctx_t ctxp;
    fmpz_mod_poly_t mod;
    fmpz_t p;
    slong d;

    fmpz_init(p);
    fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 6), 1));
    d = n_randint(state, 10) + 1;

    fmpz_mod_ctx_init(ctxp, p);
    fmpz_mod_poly_init(mod, ctxp);
    fmpz_mod_poly_randtest_monic(mod, state, d + 1, ctxp);
    fq_ctx_init_modulus(ctx, mod, ctxp, "a");

    fmpz_mod_poly_clear(mod, ctxp);
    fmpz_mod_ctx_clear(ctxp);
    fmpz_clear(p);
}
