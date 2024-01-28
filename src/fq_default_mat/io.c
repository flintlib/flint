/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fq_default_mat.h"

int fq_default_mat_fprint(FILE * file, const fq_default_mat_t mat, const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_mat_fprint(file, mat->fq_zech, ctx->ctx.fq_zech);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_mat_fprint(file, mat->fq_nmod, ctx->ctx.fq_nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        return nmod_mat_fprint(file, mat->nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        return fmpz_mod_mat_fprint(file, mat->fmpz_mod, ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        return fq_mat_fprint(file, mat->fq, ctx->ctx.fq);
    }
}

int fq_default_mat_fprint_pretty(FILE * file, const fq_default_mat_t mat, const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_mat_fprint_pretty(file, mat->fq_zech, ctx->ctx.fq_zech);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_mat_fprint_pretty(file, mat->fq_nmod, ctx->ctx.fq_nmod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        return nmod_mat_fprint_pretty(file, mat->nmod);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_mod_mat_fprint_pretty(file, mat->fmpz_mod, ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        return fq_mat_fprint_pretty(file, mat->fq, ctx->ctx.fq);
    }
}

int fq_default_mat_print(const fq_default_mat_t mat, const fq_default_ctx_t ctx) { return fq_default_mat_fprint(stdout, mat, ctx); }

int fq_default_mat_print_pretty(const fq_default_mat_t mat, const fq_default_ctx_t ctx) { return fq_default_mat_fprint_pretty(stdout, mat, ctx); }
