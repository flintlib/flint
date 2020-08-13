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

void fq_nmod_ctx_clear(fq_nmod_ctx_t ctx)
{
    nmod_poly_clear(ctx->modulus);
    nmod_poly_clear(ctx->inv);
    fmpz_clear(fq_nmod_ctx_prime(ctx));
    _nmod_vec_clear(ctx->a);
    flint_free(ctx->j);
    flint_free(ctx->var);
}
