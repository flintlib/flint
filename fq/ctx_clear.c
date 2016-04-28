/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq.h"

void
fq_ctx_clear(fq_ctx_t ctx)
{
    fmpz_mod_poly_clear(ctx->modulus);
    fmpz_mod_poly_clear(ctx->inv);
    fmpz_clear(fq_ctx_prime(ctx));
    flint_free(ctx->var);
    _fmpz_vec_clear(ctx->a, ctx->len);
    flint_free(ctx->j);
}
