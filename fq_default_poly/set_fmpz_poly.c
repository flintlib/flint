/*
    Copyright (C) 2021 William Hart
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_default_poly.h"

void fq_default_poly_set_fmpz_poly(fq_default_poly_t rop,
                              const fmpz_poly_t op, const fq_default_ctx_t ctx)
{
    fmpz_mod_ctx_t mod;
    fmpz_mod_poly_t mod_poly;
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        fmpz_mod_ctx_init(mod, fq_zech_ctx_prime(ctx->ctx.fq_zech));
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        fmpz_mod_ctx_init(mod, fq_nmod_ctx_prime(ctx->ctx.fq_nmod));
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        fmpz_mod_ctx_init_ui(mod, ctx->ctx.nmod.mod.n);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_init(mod_poly, ctx->ctx.fmpz_mod.mod);
        fmpz_mod_poly_set_fmpz_poly(mod_poly, op, ctx->ctx.fmpz_mod.mod);
        fq_default_poly_set_fmpz_mod_poly(rop, mod_poly, ctx);
        fmpz_mod_poly_clear(mod_poly, ctx->ctx.fmpz_mod.mod);
        return;
    }
    else
    {
        fmpz_mod_ctx_init(mod, fq_ctx_prime(ctx->ctx.fq));
    }
    fmpz_mod_poly_init(mod_poly, mod);
    fmpz_mod_poly_set_fmpz_poly(mod_poly, op, mod);
    fq_default_poly_set_fmpz_mod_poly(rop, mod_poly, ctx);
    fmpz_mod_poly_clear(mod_poly, mod);
    fmpz_mod_ctx_clear(mod);
}

