/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"
#include "fq.h"

/* printing *******************************************************************/

int fq_ctx_fprint(FILE * file, const fq_ctx_t ctx)
{
    int r;

    r = flint_fprintf(file, "p = ");
    if (r <= 0)
        return r;

    r = fmpz_fprint(file, fq_ctx_prime(ctx));
    if (r <= 0)
        return r;

    r = flint_fprintf(file, "\nd = %wd\n", fq_ctx_degree(ctx));
    if (r <= 0)
        return r;

    r = flint_fprintf(file, "f(X) = ");
    if (r <= 0)
        return r;

    r = fmpz_mod_poly_fprint_pretty(file, ctx->modulus, "X", ctx->ctxp);
    if (r <= 0)
        return r;

    r = flint_fprintf(file, "\n");

    return r;
}

void fq_ctx_print(const fq_ctx_t ctx) { fq_ctx_fprint(stdout, ctx); }

int fq_fprint(FILE * file, const fq_t op, const fq_ctx_t ctx) { return fmpz_poly_fprint(file, op); }
void fq_print(const fq_t op, const fq_ctx_t ctx) { fmpz_poly_print(op); }
int fq_fprint_pretty(FILE * file, const fq_t op, const fq_ctx_t ctx) { return fmpz_poly_fprint_pretty(file, op, ctx->var); }
int fq_print_pretty(const fq_t op, const fq_ctx_t ctx) { return fmpz_poly_print_pretty(op, ctx->var); }
