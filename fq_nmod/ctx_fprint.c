/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fmpz.h"
#include "flint-impl.h"
#include "fq_nmod.h"

int fq_nmod_ctx_fprint(FILE * file, const fq_nmod_ctx_t ctx)
{
    int r;
    slong i, k;

    r = fprintf(file, "p = ");
    if (r <= 0)
        return r;

    r = fmpz_fprint(file, fq_nmod_ctx_prime(ctx));
    if (r <= 0)
        return r;

    r = fprintf(file, "\nd = " WORD_FMT "d\nf(X) = ", ctx->j[ctx->len - 1]);
    if (r <= 0)
        return r;

    r = fprintf(file, WORD_FMT "u", ctx->a[0]);
    if (r <= 0)
        return r;

    for (k = 1; k < ctx->len; k++)
    {
        i = ctx->j[k];
        r = fprintf(file, " + ");
        if (r <= 0)
            return r;

        if (ctx->a[k] == UWORD(1))
        {
            if (i == 1)
                r = fprintf(file, "X");
            else
                r = fprintf(file, "X^" WORD_FMT "d", i);
            if (r <= 0)
                return r;
        }
        else
        {
            r = fprintf(file, WORD_FMT "u", ctx->a[k]);
            if (r <= 0)
                return r;

            if (i == 1)
                r = fprintf(file, "*X");
            else
                r = fprintf(file, "*X^" WORD_FMT "d", i);
            if (r <= 0)
                return r;
        }
    }
    r = fprintf(file, "\n");
    return r;
}
