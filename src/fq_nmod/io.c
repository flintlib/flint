/*
    Copyright (C) 2023, 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "nmod_poly.h"
#include "fmpz.h"
#include "fq_nmod.h"

/* printing *******************************************************************/

/* TODO */
int fq_nmod_ctx_fprint(FILE * file, const fq_nmod_ctx_t ctx)
{
    int r;
    slong i, k;

    r = flint_fprintf(file, "p = %wu", fq_nmod_ctx_prime(ctx));
    if (r <= 0)
        return r;

    r = flint_fprintf(file, "\nd = %wd\nf(X) = ", ctx->j[ctx->len - 1]);
    if (r <= 0)
        return r;

    r = flint_fprintf(file, "%wu", ctx->a[0]);
    if (r <= 0)
        return r;

    for (k = 1; k < ctx->len; k++)
    {
        i = ctx->j[k];
        r = flint_fprintf(file, " + ");
        if (r <= 0)
            return r;

        if (ctx->a[k] == UWORD(1))
        {
            if (i == 1)
                r = flint_fprintf(file, "X");
            else
                r = flint_fprintf(file, "X^%wd", i);
            if (r <= 0)
                return r;
        }
        else
        {
            r = flint_fprintf(file, "%wu", ctx->a[k]);
            if (r <= 0)
                return r;

            if (i == 1)
                r = flint_fprintf(file, "*X");
            else
                r = flint_fprintf(file, "*X^%wd", i);
            if (r <= 0)
                return r;
        }
    }
    r = flint_fprintf(file, "\n");
    return r;
}

void fq_nmod_ctx_print(const fq_nmod_ctx_t ctx) { fq_nmod_ctx_fprint(stdout, ctx); }

int fq_nmod_fprint(FILE * file, const fq_nmod_t op, const fq_nmod_ctx_t ctx) { return nmod_poly_fprint(file, op); }
void fq_nmod_print(const fq_nmod_t op, const fq_nmod_ctx_t ctx) { nmod_poly_print(op); }
int fq_nmod_fprint_pretty(FILE * file, const fq_nmod_t op, const fq_nmod_ctx_t ctx) { return nmod_poly_fprint_pretty(file, op, ctx->var); }
void fq_nmod_print_pretty(const fq_nmod_t op, const fq_nmod_ctx_t ctx) { nmod_poly_print_pretty(op, ctx->var); }
