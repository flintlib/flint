/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fq_nmod.h"
#include "fq_zech.h"

/* printing *******************************************************************/

int
fq_zech_ctx_fprint(FILE * file, const fq_zech_ctx_t ctx)
{
    int r;
    r = flint_fprintf(file, "Zech Representation:\n");
    if (r <= 0)
        return r;
    return fq_nmod_ctx_fprint(file, ctx->fq_nmod_ctx);
}

void fq_zech_ctx_print(const fq_zech_ctx_t ctx) { fq_zech_ctx_fprint(stdout, ctx); }

int fq_zech_fprint(FILE * file, const fq_zech_t op, const fq_zech_ctx_t FLINT_UNUSED(ctx)) { return flint_fprintf(file, "%wd", op->value); }
void fq_zech_print(const fq_zech_t op, const fq_zech_ctx_t ctx) { fq_zech_fprint(stdout, op, ctx); }
int fq_zech_fprint_pretty(FILE * file, const fq_zech_t op, const fq_zech_ctx_t ctx) { return flint_fprintf(file, "%s^%wd", ctx->fq_nmod_ctx->var, op->value); }
void fq_zech_print_pretty(const fq_zech_t op, const fq_zech_ctx_t ctx) { fq_zech_fprint_pretty(stdout, op, ctx); }
