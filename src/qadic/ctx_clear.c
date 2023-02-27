/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "padic.h"
#include "qadic.h"

void qadic_ctx_clear(qadic_ctx_t ctx)
{
    padic_ctx_clear(&ctx->pctx);
    _fmpz_vec_clear(ctx->a, ctx->len);
    flint_free(ctx->j);
    flint_free(ctx->var);
}

