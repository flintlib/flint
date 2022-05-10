/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "flint-impl.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "fq.h"

int fq_ctx_fprint(FILE * file, const fq_ctx_t ctx)
{
    int r;

    r = fprintf(file, "p = ");
    if (r <= 0)
        return r;

    r = fmpz_fprint(file, fq_ctx_prime(ctx));
    if (r <= 0)
        return r;

    r = fprintf(file, "\nd = " WORD_FMT "d\n", fq_ctx_degree(ctx));
    if (r <= 0)
        return r;

    r = fprintf(file, "f(X) = ");
    if (r <= 0)
        return r;

    r = fmpz_mod_poly_fprint_pretty(file, ctx->modulus, "X", ctx->ctxp);
    if (r <= 0)
        return r;

    r = fprintf(file, "\n");

    return r;
}
