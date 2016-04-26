/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
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

    fq_nmod_ctx_struct * fq_nmod_ctx;

    fq_nmod_ctx = flint_malloc(sizeof(fq_nmod_ctx_struct));

    fmpz_init(p);
    fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 4), 1));
    max_d = floor(log(n_pow(2, 16)) / log(fmpz_get_ui(p)));
    d = n_randint(state, max_d - 1) + 2;
    fq_nmod_ctx_init(fq_nmod_ctx, p, d, "a");
    fq_zech_ctx_init_fq_nmod_ctx(ctx, fq_nmod_ctx);
    fmpz_clear(p);

    ctx->owns_fq_nmod_ctx = 1;
}
