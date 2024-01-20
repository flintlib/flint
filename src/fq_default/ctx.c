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

void fq_default_ctx_init_modulus_type(fq_default_ctx_t ctx,
                const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx,
                                                    const char * var, int type)
{
    fmpz const * p = fmpz_mod_ctx_modulus(mod_ctx);
    int bits = fmpz_bits(p);
    int d = fmpz_mod_poly_degree(modulus, mod_ctx);

    if (type == FQ_DEFAULT_FQ_ZECH || (type == 0 && d > 1 && bits*d <= 16))
    {
        nmod_poly_t nmodulus;
        fq_nmod_ctx_struct * fq_nmod_ctx;
        ctx->type = FQ_DEFAULT_FQ_ZECH;
        nmod_poly_init(nmodulus, fmpz_get_ui(p));
        fmpz_mod_poly_get_nmod_poly(nmodulus, modulus);
        fq_nmod_ctx = flint_malloc(sizeof(fq_nmod_ctx_struct));
        fq_nmod_ctx_init_modulus(fq_nmod_ctx, nmodulus, var);
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
        nmod_poly_clear(nmodulus);
    }
    else if (type == FQ_DEFAULT_FQ_NMOD || (type == 0 && d > 1 && fmpz_abs_fits_ui(p)))
    {
        nmod_poly_t nmodulus;
        ctx->type = FQ_DEFAULT_FQ_NMOD;
        nmod_poly_init(nmodulus, fmpz_get_ui(p));
        fmpz_mod_poly_get_nmod_poly(nmodulus, modulus);
        fq_nmod_ctx_init_modulus(ctx->ctx.fq_nmod, nmodulus, var);
        nmod_poly_clear(nmodulus);
    }
    else if (type == FQ_DEFAULT_NMOD || (type == 0 && d == 1 && fmpz_abs_fits_ui(p)))
    {
        mp_limb_t c0, c1;
        ctx->type = FQ_DEFAULT_NMOD;
        nmod_init(&ctx->ctx.nmod.mod, fmpz_get_ui(p));
        c0 = fmpz_get_ui(modulus->coeffs + 0);
        c1 = fmpz_get_ui(modulus->coeffs + 1);
        c0 = nmod_neg(c0, ctx->ctx.nmod.mod);
        ctx->ctx.nmod.a = nmod_div(c0, c1, ctx->ctx.nmod.mod);
    }
    else if (type == FQ_DEFAULT_FMPZ_MOD || (type == 0 && d == 1))
    {
        ctx->type = FQ_DEFAULT_FMPZ_MOD;
        fmpz_mod_ctx_init(ctx->ctx.fmpz_mod.mod, p);
        fmpz_init(ctx->ctx.fmpz_mod.a);
        fmpz_mod_divides(ctx->ctx.fmpz_mod.a, modulus->coeffs + 0,
                                   modulus->coeffs + 1, ctx->ctx.fmpz_mod.mod);
        fmpz_mod_neg(ctx->ctx.fmpz_mod.a, ctx->ctx.fmpz_mod.a,
                                                        ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        ctx->type = FQ_DEFAULT_FQ;
        fq_ctx_init_modulus(ctx->ctx.fq, modulus, mod_ctx, var);
    }
}

void fq_default_ctx_init_modulus(fq_default_ctx_t ctx,
       const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx, const char * var)
{
    fq_default_ctx_init_modulus_type(ctx, modulus, mod_ctx, var, 0);
}

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
        mp_limb_t c0, c1;
        ctx->type = FQ_DEFAULT_NMOD;
        nmod_init(&ctx->ctx.nmod.mod, p);
        c0 = modulus->coeffs[0];
        c1 = modulus->coeffs[1];
        c0 = nmod_neg(c0, ctx->ctx.nmod.mod);
        ctx->ctx.nmod.a = nmod_div(c0, c1, ctx->ctx.nmod.mod);
    }
    else if (type == FQ_DEFAULT_FMPZ_MOD || (type == 0 && d == 1))
    {
        mp_limb_t c0, c1;
        ctx->type = FQ_DEFAULT_FMPZ_MOD;
        fmpz_mod_ctx_init_ui(ctx->ctx.fmpz_mod.mod, p);
        fmpz_init_set_ui(ctx->ctx.fmpz_mod.a, 0);
        c0 = modulus->coeffs[0];
        c1 = modulus->coeffs[1];
        c0 = nmod_neg(c0, modulus->mod);
        fmpz_set_ui(ctx->ctx.fmpz_mod.a, nmod_div(c0, c1, modulus->mod));
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
