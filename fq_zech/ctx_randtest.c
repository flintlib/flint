/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech.h"
#include "fq_nmod.h"
#include "ulong_extras.h"
#include "long_extras.h"
#include <math.h>

void
fq_zech_ctx_randtest(fq_zech_ctx_t ctx, flint_rand_t state)
{
    fmpz_t p;
    slong max_d, d;

    fmpz_init(p);
    fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 4), 1));
    if (n_randlimb(state) % 16 == 0)   /* slow */
        max_d = floor(log(n_pow(2, 15)) / log(fmpz_get_ui(p)));
    else
        max_d = floor(log(n_pow(2, 11)) / log(fmpz_get_ui(p)));
    d = n_randint(state, max_d - 1) + 2;
    fq_zech_ctx_init_random(ctx, p, d, "a");
    fmpz_clear(p);

    ctx->owns_fq_nmod_ctx = 1;
}
