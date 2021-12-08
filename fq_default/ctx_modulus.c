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

void fq_default_ctx_modulus(fmpz_mod_poly_t p, const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        nmod_poly_struct const * mod = fq_zech_ctx_modulus(ctx->ctx.fq_zech);
        fmpz_mod_poly_set_nmod_poly(p, mod);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        nmod_poly_struct const * mod = fq_nmod_ctx_modulus(ctx->ctx.fq_nmod);
        fmpz_mod_poly_set_nmod_poly(p, mod);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        _fmpz_mod_poly_fit_length(p, 2);
        _fmpz_mod_poly_set_length(p, 2);
        fmpz_set_ui(p->coeffs + 0, nmod_neg(ctx->ctx.nmod.a, ctx->ctx.nmod.mod));
        fmpz_one(p->coeffs + 1);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        _fmpz_mod_poly_fit_length(p, 2);
        _fmpz_mod_poly_set_length(p, 2);
        fmpz_mod_neg(p->coeffs + 0, ctx->ctx.fmpz_mod.a, ctx->ctx.fmpz_mod.mod);
        fmpz_one(p->coeffs + 1);
    }
    else
    {
        fmpz_mod_ctx_struct const * mod = ctx->ctx.fq->ctxp;
        fmpz_mod_poly_set(p, ctx->ctx.fq->modulus, mod);
    }
}

