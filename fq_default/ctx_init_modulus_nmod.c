/*
    Copyright (C) 2021 William Hart
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_default.h"

void fq_default_ctx_init_modulus_nmod_type(fq_default_ctx_t ctx,
                         const nmod_poly_t modulus, const char * var, int type)
{
    ulong p = modulus->mod.n;
    int bits = FLINT_BIT_COUNT(p);
    int d = nmod_poly_degree(modulus);

    if (type == FQ_DEFAULT_FQ_ZECH || (type == 0 && d > 1 && bits*d <= 16))
    {
        fq_nmod_ctx_struct * fq_nmod_ctx =
                                     flint_malloc(sizeof(fq_nmod_ctx_struct));
        ctx->type = FQ_DEFAULT_FQ_ZECH;
        fq_nmod_ctx_init_modulus(fq_nmod_ctx, modulus, var);
        if (fq_zech_ctx_init_fq_nmod_ctx_check(ctx->ctx.fq_zech, fq_nmod_ctx))
        {
            ctx->ctx.fq_zech->owns_fq_nmod_ctx = 1;
        }
        else
        {
            *ctx->ctx.fq_nmod = *fq_nmod_ctx;
            flint_free(fq_nmod_ctx);
            ctx->type = FQ_DEFAULT_FQ_NMOD;
        }
    }
    else if (type == FQ_DEFAULT_FQ_NMOD || (type == 0 && d > 1))
    {
        ctx->type = FQ_DEFAULT_FQ_NMOD;
        fq_nmod_ctx_init_modulus(ctx->ctx.fq_nmod, modulus, var);
    }
    else if (type == FQ_DEFAULT_NMOD || (type == 0 && d == 1))
    {
        ctx->type = FQ_DEFAULT_NMOD;
        nmod_init(&ctx->ctx.nmod.mod, p);
        ctx->ctx.nmod.a = 0;
    }
    else if (type == FQ_DEFAULT_FMPZ_MOD || (type == 0 && d == 1))
    {
        ctx->type = FQ_DEFAULT_FMPZ_MOD;
        fmpz_mod_ctx_init_ui(ctx->ctx.fmpz_mod.mod, p);
        fmpz_init_set_ui(ctx->ctx.fmpz_mod.a, 0);
    }
    else
    {
        fmpz_mod_ctx_t fmod_ctx;
        fmpz_mod_poly_t fmod;
        fmpz_t p;
        ctx->type = FQ_DEFAULT_FQ;
        fmpz_init_set_ui(p, modulus->mod.n);
        fmpz_mod_ctx_init(fmod_ctx, p);
        fmpz_mod_poly_init(fmod, fmod_ctx);
        fmpz_mod_poly_set_nmod_poly(fmod, modulus);
        fq_ctx_init_modulus(ctx->ctx.fq, fmod, fmod_ctx, var);
        fmpz_mod_poly_clear(fmod, fmod_ctx);
        fmpz_mod_ctx_clear(fmod_ctx);
        fmpz_clear(p);
    }
}

void fq_default_ctx_init_modulus_nmod(fq_default_ctx_t ctx,
                                   const nmod_poly_t modulus, const char * var)
{
    fq_default_ctx_init_modulus_nmod_type(ctx, modulus, var, 0);
}

