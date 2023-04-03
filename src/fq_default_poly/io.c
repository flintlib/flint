/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fq_default_poly.h"

/* printing *******************************************************************/

int fq_default_poly_fprint(FILE * file, const fq_default_poly_t poly, const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_poly_fprint(file, poly->fq_zech, ctx->ctx.fq_zech);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_poly_fprint(file, poly->fq_nmod, ctx->ctx.fq_nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        return nmod_poly_fprint(file, poly->nmod);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_mod_poly_fprint(file, poly->fmpz_mod, ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        return fq_poly_fprint(file, poly->fq, ctx->ctx.fq);
    }
}

int fq_default_poly_fprint_pretty(FILE * file, const fq_default_poly_t poly, const char *x, const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_poly_fprint_pretty(file,
                                           poly->fq_zech, x, ctx->ctx.fq_zech);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_poly_fprint_pretty(file,
                                           poly->fq_nmod, x, ctx->ctx.fq_nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        return nmod_poly_fprint_pretty(file, poly->nmod, x);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_mod_poly_fprint_pretty(file, poly->fmpz_mod, x,
                                                        ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        return fq_poly_fprint_pretty(file, poly->fq, x, ctx->ctx.fq);
    }
}

int fq_default_poly_print(const fq_default_poly_t poly, const fq_default_ctx_t ctx) { return fq_default_poly_fprint(stdout, poly, ctx); }
int fq_default_poly_print_pretty(const fq_default_poly_t poly, const char *x, const fq_default_ctx_t ctx) { return fq_default_poly_fprint_pretty(stdout, poly, x, ctx); }
