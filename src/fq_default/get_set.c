/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "fq_default.h"

void fq_default_get_fmpz_mod_poly(fmpz_mod_poly_t poly,
                             const fq_default_t op, const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        nmod_poly_t p;
        ulong mod = fq_zech_ctx_prime(ctx->ctx.fq_zech);
        nmod_poly_init(p, mod);
        fq_zech_get_nmod_poly(p, op->fq_zech, ctx->ctx.fq_zech);
        fmpz_mod_poly_set_nmod_poly(poly, p);
        nmod_poly_clear(p);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        nmod_poly_t p;
        ulong mod = fq_nmod_ctx_prime(ctx->ctx.fq_nmod);
        nmod_poly_init(p, mod);
        fq_nmod_get_nmod_poly(p, op->fq_nmod, ctx->ctx.fq_nmod);
        fmpz_mod_poly_set_nmod_poly(poly, p);
        nmod_poly_clear(p);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        _fmpz_mod_poly_fit_length(poly, 1);
        fmpz_set_ui(poly->coeffs + 0, op->nmod);
        _fmpz_mod_poly_set_length(poly, op->nmod != 0);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        _fmpz_mod_poly_fit_length(poly, 1);
        fmpz_set(poly->coeffs + 0, op->fmpz_mod);
        _fmpz_mod_poly_set_length(poly, !fmpz_is_zero(op->fmpz_mod));
    }
    else
    {
        fq_get_fmpz_mod_poly(poly, op->fq, ctx->ctx.fq);
    }
}

void fq_default_get_fmpz_poly(fmpz_poly_t poly,
                             const fq_default_t op, const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        nmod_poly_t p;
        ulong mod = fq_zech_ctx_prime(ctx->ctx.fq_zech);
        nmod_poly_init(p, mod);
        fq_zech_get_nmod_poly(p, op->fq_zech, ctx->ctx.fq_zech);
        fmpz_poly_set_nmod_poly(poly, p);
        nmod_poly_clear(p);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        nmod_poly_t p;
        ulong mod = fq_nmod_ctx_prime(ctx->ctx.fq_nmod);
        nmod_poly_init(p, mod);
        fq_nmod_get_nmod_poly(p, op->fq_nmod, ctx->ctx.fq_nmod);
        fmpz_poly_set_nmod_poly(poly, p);
        nmod_poly_clear(p);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        fmpz_poly_fit_length(poly, 1);
        fmpz_set_ui(poly->coeffs + 0, op->nmod);
        _fmpz_poly_set_length(poly, op->nmod != 0);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_poly_fit_length(poly, 1);
        fmpz_set(poly->coeffs + 0, op->fmpz_mod);
        _fmpz_poly_set_length(poly, !fmpz_is_zero(op->fmpz_mod));
    }
    else
    {
        fq_get_fmpz_poly(poly, op->fq, ctx->ctx.fq);
    }
}

void fq_default_set_fmpz_mod_poly(fq_default_t op,
                        const fmpz_mod_poly_t poly, const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        nmod_poly_t p;
        ulong mod = fq_zech_ctx_prime(ctx->ctx.fq_zech);
        nmod_poly_init(p, mod);
        fmpz_mod_poly_get_nmod_poly(p, poly);
        fq_zech_set_nmod_poly(op->fq_zech, p, ctx->ctx.fq_zech);
        nmod_poly_clear(p);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        nmod_poly_t p;
        ulong mod = fq_nmod_ctx_prime(ctx->ctx.fq_nmod);
        nmod_poly_init(p, mod);
        fmpz_mod_poly_get_nmod_poly(p, poly);
        fq_nmod_set_nmod_poly(op->fq_nmod, p, ctx->ctx.fq_nmod);
        nmod_poly_clear(p);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        nmod_poly_t p;
        nmod_poly_init_mod(p, ctx->ctx.nmod.mod);
        fmpz_mod_poly_get_nmod_poly(p, poly);
        op->nmod = nmod_poly_evaluate_nmod(p, ctx->ctx.nmod.a);
        nmod_poly_clear(p);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_evaluate_fmpz(op->fmpz_mod, poly, ctx->ctx.fmpz_mod.a,
                                                        ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        fq_set_fmpz_mod_poly(op->fq, poly, ctx->ctx.fq);
    }
}

void fq_default_set_fmpz_poly(fq_default_t op,
                            const fmpz_poly_t poly, const fq_default_ctx_t ctx)
{
    if (ctx->type == FQ_DEFAULT_FQ_ZECH)
    {
        nmod_poly_t p;
        ulong mod = fq_zech_ctx_prime(ctx->ctx.fq_zech);
        nmod_poly_init(p, mod);
        fmpz_poly_get_nmod_poly(p, poly);
        fq_zech_set_nmod_poly(op->fq_zech, p, ctx->ctx.fq_zech);
        nmod_poly_clear(p);
    }
    else if (ctx->type == FQ_DEFAULT_FQ_NMOD)
    {
        nmod_poly_t p;
        ulong mod = fq_nmod_ctx_prime(ctx->ctx.fq_nmod);
        nmod_poly_init(p, mod);
        fmpz_poly_get_nmod_poly(p, poly);
        fq_nmod_set_nmod_poly(op->fq_nmod, p, ctx->ctx.fq_nmod);
        nmod_poly_clear(p);
    }
    else if (ctx->type == FQ_DEFAULT_NMOD)
    {
        nmod_poly_t p;
        nmod_poly_init_mod(p, ctx->ctx.nmod.mod);
        fmpz_poly_get_nmod_poly(p, poly);
        op->nmod = nmod_poly_evaluate_nmod(p, ctx->ctx.nmod.a);
        nmod_poly_clear(p);
    }
    else if (ctx->type == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_poly_t p;
        fmpz_mod_poly_init(p, ctx->ctx.fmpz_mod.mod);
        fmpz_mod_poly_set_fmpz_poly(p, poly, ctx->ctx.fmpz_mod.mod);
        fmpz_mod_poly_evaluate_fmpz(op->fmpz_mod, p, ctx->ctx.fmpz_mod.a,
                                                        ctx->ctx.fmpz_mod.mod);
        fmpz_mod_poly_clear(p, ctx->ctx.fmpz_mod.mod);
    }
    else
    {
        fq_set_fmpz_poly(op->fq, poly, ctx->ctx.fq);
    }
}
