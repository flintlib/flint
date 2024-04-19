/*
    Copyright (C) 2021 William Hart
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_default.h"

void fq_default_ctx_init_type(fq_default_ctx_t ctx,
                            const fmpz_t p, slong d, const char *var, int type)
{
    int bits = fmpz_bits(p);

    if (type == FQ_DEFAULT_FQ_ZECH || (type == 0 && d > 1 && bits*d <= 16))
    {
        gr_ctx_init_fq_zech(FQ_DEFAULT_GR_CTX(ctx), *p, d, var);
    }
    else if (type == FQ_DEFAULT_FQ_NMOD || (type == 0 && d > 1 && fmpz_abs_fits_ui(p)))
    {
        gr_ctx_init_fq_nmod(FQ_DEFAULT_GR_CTX(ctx), *p, d, var);
    }
    else if (type == FQ_DEFAULT_NMOD || (type == 0 && d == 1 && fmpz_abs_fits_ui(p)))
    {
        gr_ctx_init_nmod(FQ_DEFAULT_GR_CTX(ctx), fmpz_get_ui(p));
        NMOD_CTX_A(FQ_DEFAULT_GR_CTX(ctx))[0] = 0;
        GR_MUST_SUCCEED(gr_ctx_set_is_field(FQ_DEFAULT_GR_CTX(ctx), T_TRUE));
    }
    else if (type == FQ_DEFAULT_FMPZ_MOD || (type == 0 && d == 1))
    {
        gr_ctx_init_fmpz_mod(FQ_DEFAULT_GR_CTX(ctx), p);
        GR_MUST_SUCCEED(gr_ctx_set_is_field(FQ_DEFAULT_GR_CTX(ctx), T_TRUE));
    }
    else
    {
        gr_ctx_init_fq(FQ_DEFAULT_GR_CTX(ctx), p, d, var);
    }
}

int gr_ctx_init_fq_zech_modulus_nmod_poly(gr_ctx_t ctx, const nmod_poly_t modulus, const char * var);
int gr_ctx_init_fq_zech_modulus_fmpz_mod_poly(gr_ctx_t ctx, const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx, const char * var);

int gr_ctx_init_fq_nmod_modulus_nmod_poly(gr_ctx_t ctx, const nmod_poly_t modulus, const char * var);
int gr_ctx_init_fq_nmod_modulus_fmpz_mod_poly(gr_ctx_t ctx, const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx, const char * var);

int gr_ctx_init_fq_modulus_nmod_poly(gr_ctx_t ctx, const nmod_poly_t modulus, const char * var);
int gr_ctx_init_fq_modulus_fmpz_mod_poly(gr_ctx_t ctx, const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx, const char * var);


void fq_default_ctx_init_modulus_type(fq_default_ctx_t ctx,
                const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx,
                                                    const char * var, int type)
{
    fmpz const * p = fmpz_mod_ctx_modulus(mod_ctx);
    int bits = fmpz_bits(p);
    int d = fmpz_mod_poly_degree(modulus, mod_ctx);

    if (type == FQ_DEFAULT_FQ_ZECH || (type == 0 && d > 1 && bits*d <= 16))
    {
        if (gr_ctx_init_fq_zech_modulus_fmpz_mod_poly(FQ_DEFAULT_GR_CTX(ctx), modulus, mod_ctx, var) != GR_SUCCESS)
            fq_default_ctx_init_modulus_type(ctx, modulus, mod_ctx, var, _FQ_DEFAULT_FQ_NMOD);
    }
    else if (type == FQ_DEFAULT_FQ_NMOD || (type == 0 && d > 1 && fmpz_abs_fits_ui(p)))
    {
        GR_MUST_SUCCEED(gr_ctx_init_fq_nmod_modulus_fmpz_mod_poly(FQ_DEFAULT_GR_CTX(ctx), modulus, mod_ctx, var));
    }
    else if (type == FQ_DEFAULT_NMOD || (type == 0 && d == 1 && fmpz_abs_fits_ui(p)))
    {
        mp_limb_t c0, c1, a;
        nmod_t mod;

        nmod_init(&mod, fmpz_get_ui(p));
        c0 = fmpz_get_ui(modulus->coeffs + 0);
        c1 = fmpz_get_ui(modulus->coeffs + 1);
        c0 = nmod_neg(c0, mod);
        a = nmod_div(c0, c1, mod);

        _gr_ctx_init_nmod(FQ_DEFAULT_GR_CTX(ctx), &mod);
        NMOD_CTX_A(FQ_DEFAULT_GR_CTX(ctx))[0] = a;
        GR_MUST_SUCCEED(gr_ctx_set_is_field(FQ_DEFAULT_GR_CTX(ctx), T_TRUE));
    }
    else if (type == FQ_DEFAULT_FMPZ_MOD || (type == 0 && d == 1))
    {
        gr_ctx_init_fmpz_mod(FQ_DEFAULT_GR_CTX(ctx), p);
        GR_MUST_SUCCEED(gr_ctx_set_is_field(FQ_DEFAULT_GR_CTX(ctx), T_TRUE));

        fmpz_mod_divides(FQ_DEFAULT_CTX_FMPZ_MOD_A(ctx), modulus->coeffs + 0,
            modulus->coeffs + 1, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
        fmpz_mod_neg(FQ_DEFAULT_CTX_FMPZ_MOD_A(ctx), FQ_DEFAULT_CTX_FMPZ_MOD_A(ctx),
                                                        FQ_DEFAULT_CTX_FMPZ_MOD(ctx));

        fmpz_set(FMPZ_MOD_CTX_A(FQ_DEFAULT_GR_CTX(ctx)), FQ_DEFAULT_CTX_FMPZ_MOD_A(ctx));
    }
    else
    {
        GR_MUST_SUCCEED(gr_ctx_init_fq_modulus_fmpz_mod_poly(FQ_DEFAULT_GR_CTX(ctx), modulus, mod_ctx, var));
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
        if (gr_ctx_init_fq_zech_modulus_nmod_poly(FQ_DEFAULT_GR_CTX(ctx), modulus, var) != GR_SUCCESS)
            fq_default_ctx_init_modulus_nmod_type(ctx, modulus, var, _FQ_DEFAULT_FQ_NMOD);
    }
    else if (type == FQ_DEFAULT_FQ_NMOD || (type == 0 && d > 1))
    {
        GR_MUST_SUCCEED(gr_ctx_init_fq_nmod_modulus_nmod_poly(FQ_DEFAULT_GR_CTX(ctx), modulus, var));
    }
    else if (type == FQ_DEFAULT_NMOD || (type == 0 && d == 1))
    {
        mp_limb_t c0, c1, a;
        nmod_t mod;

        nmod_init(&mod, p);
        c0 = modulus->coeffs[0];
        c1 = modulus->coeffs[1];
        c0 = nmod_neg(c0, mod);
        a = nmod_div(c0, c1, mod);

        _gr_ctx_init_nmod(FQ_DEFAULT_GR_CTX(ctx), &mod);
        NMOD_CTX_A(FQ_DEFAULT_GR_CTX(ctx))[0] = a;
        GR_MUST_SUCCEED(gr_ctx_set_is_field(FQ_DEFAULT_GR_CTX(ctx), T_TRUE));
    }
    else if (type == FQ_DEFAULT_FMPZ_MOD || (type == 0 && d == 1))
    {
        fmpz_t pp;
        mp_limb_t c0, c1, a;

        c0 = modulus->coeffs[0];
        c1 = modulus->coeffs[1];
        c0 = nmod_neg(c0, modulus->mod);
        a = nmod_div(c0, c1, modulus->mod);

        fmpz_init_set_ui(pp, p);
        gr_ctx_init_fmpz_mod(FQ_DEFAULT_GR_CTX(ctx), pp);
        fmpz_clear(pp);
        GR_MUST_SUCCEED(gr_ctx_set_is_field(FQ_DEFAULT_GR_CTX(ctx), T_TRUE));
        fmpz_set_ui(FMPZ_MOD_CTX_A(FQ_DEFAULT_GR_CTX(ctx)), a);
    }
    else
    {
        GR_MUST_SUCCEED(gr_ctx_init_fq_modulus_nmod_poly(FQ_DEFAULT_GR_CTX(ctx), modulus, var));
    }
}

void fq_default_ctx_init_modulus_nmod(fq_default_ctx_t ctx,
                                   const nmod_poly_t modulus, const char * var)
{
    fq_default_ctx_init_modulus_nmod_type(ctx, modulus, var, 0);
}

void fq_default_ctx_modulus(fmpz_mod_poly_t p, const fq_default_ctx_t ctx)
{
    if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_ZECH)
    {
        nmod_poly_struct const * mod = fq_zech_ctx_modulus(FQ_DEFAULT_CTX_FQ_ZECH(ctx));
        fmpz_mod_poly_set_nmod_poly(p, mod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FQ_NMOD)
    {
        nmod_poly_struct const * mod = fq_nmod_ctx_modulus(FQ_DEFAULT_CTX_FQ_NMOD(ctx));
        fmpz_mod_poly_set_nmod_poly(p, mod);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_NMOD)
    {
        _fmpz_mod_poly_fit_length(p, 2);
        _fmpz_mod_poly_set_length(p, 2);
        fmpz_set_ui(p->coeffs + 0, nmod_neg(FQ_DEFAULT_CTX_NMOD_A(ctx)[0], FQ_DEFAULT_CTX_NMOD(ctx)));
        fmpz_one(p->coeffs + 1);
    }
    else if (_FQ_DEFAULT_TYPE(ctx) == _FQ_DEFAULT_FMPZ_MOD)
    {
        _fmpz_mod_poly_fit_length(p, 2);
        _fmpz_mod_poly_set_length(p, 2);
        fmpz_mod_neg(p->coeffs + 0, FQ_DEFAULT_CTX_FMPZ_MOD_A(ctx), FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
        fmpz_one(p->coeffs + 1);
    }
    else
    {
        fmpz_mod_ctx_struct const * mod = FQ_DEFAULT_CTX_FQ(ctx)->ctxp;
        fmpz_mod_poly_set(p, FQ_DEFAULT_CTX_FQ(ctx)->modulus, mod);
    }
}
