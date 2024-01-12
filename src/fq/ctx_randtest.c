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

#include "fq_nmod.h"
#include "fq.h"

void fq_ctx_init_set_clear_small_fq_nmod_ctx(fq_ctx_t, fq_nmod_ctx_t);

/* FIXME: Test big primes? */
void
fq_ctx_randtest(fq_ctx_t ctx, flint_rand_t state)
{
    fq_nmod_ctx_t nmod_ctx;

    fq_nmod_ctx_randtest(nmod_ctx, state);

    fq_ctx_init_set_clear_small_fq_nmod_ctx(ctx, nmod_ctx);
}
