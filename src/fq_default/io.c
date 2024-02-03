/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fq_default.h"

#include "gr.h"

/* printing *******************************************************************/

int fq_default_ctx_fprint(FILE * file, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_ctx_fprint(file, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_ctx_fprint(file, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        return flint_fprintf(file, "p = %wu\n", FQ_DEFAULT_CTX_NMOD(ctx).n);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        int z = flint_fprintf(file, "p = ");
        if (z <= 0)
            return z;
        z = fmpz_fprint(file, fmpz_mod_ctx_modulus(FQ_DEFAULT_CTX_FMPZ_MOD(ctx)));
        if (z <= 0)
            return z;
        return flint_fprintf(file, "\n");
    }
    else
    {
        return fq_ctx_fprint(file, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

int fq_default_fprint(FILE * file, const fq_default_t op, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_fprint(file, op->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_fprint(file, op->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        return flint_fprintf(file, "%wu", op->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_fprint(file, op->fmpz_mod);
    }
    else
    {
        return fq_fprint(file, op->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

int fq_default_fprint_pretty(FILE * file, const fq_default_t op, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_fprint_pretty(file, op->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_fprint_pretty(file, op->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        return flint_fprintf(file, "%wu", op->nmod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_fprint(file, op->fmpz_mod);
    }
    else
    {
        return fq_fprint_pretty(file, op->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

void fq_default_ctx_print(const fq_default_ctx_t ctx) { fq_default_ctx_fprint(stdout, ctx); }
void fq_default_print(const fq_default_t op, const fq_default_ctx_t ctx) { fq_default_fprint(stdout, op, ctx); }
void fq_default_print_pretty(const fq_default_t op, const fq_default_ctx_t ctx) { fq_default_fprint_pretty(stdout, op, ctx); }
