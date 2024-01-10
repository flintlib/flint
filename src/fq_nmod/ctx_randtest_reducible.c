/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_poly.h"
#include "fmpz.h"
#include "fq_nmod.h"

void
fq_nmod_ctx_randtest_reducible(fq_nmod_ctx_t ctx, flint_rand_t state)
{
    ulong p;
    slong d;
    nmod_poly_t mod;

    p = n_randprime(state, 2 + n_randint(state, 6), 1);
    d = n_randint(state, 10) + 1;
    nmod_poly_init(mod, p);
    nmod_poly_randtest_monic(mod, state, d + 1);
    fq_nmod_ctx_init_modulus(ctx, mod, "a");
    nmod_poly_clear(mod);
}
