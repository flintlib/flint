/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_poly.h"
#include "fq_nmod.h"

void fq_nmod_set(fq_nmod_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
{
    nmod_poly_set(rop, op);
}

void fq_nmod_set_si(fq_nmod_t rop, const slong x, const fq_nmod_ctx_t ctx)
{
    mp_limb_t rx = x < 0 ? -x : x;
    rx =  n_mod2_preinv(rx, ctx->mod.n, ctx->mod.ninv);
    if (x < 0)
        rx = ctx->mod.n - rx;

    nmod_poly_zero(rop);
    nmod_poly_set_coeff_ui(rop, 0, rx);
}

void fq_nmod_set_ui(fq_nmod_t rop, const ulong x, const fq_nmod_ctx_t ctx)
{
    nmod_poly_zero(rop);
    nmod_poly_set_coeff_ui(rop, 0, n_mod2_preinv(x, ctx->mod.n, ctx->mod.ninv));
}

void fq_nmod_swap(fq_nmod_t op1, fq_nmod_t op2, const fq_nmod_ctx_t ctx)
{
    nmod_poly_swap(op1, op2);
}

void fq_nmod_zero(fq_nmod_t rop,  const fq_nmod_ctx_t ctx)
{
    nmod_poly_zero(rop);
}

void fq_nmod_one(fq_nmod_t rop,  const fq_nmod_ctx_t ctx)
{
    nmod_poly_one(rop);
}

void fq_nmod_gen(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
{
    if (ctx->modulus->length == 2)
    {
        nmod_poly_set_coeff_ui(rop, 0,
              nmod_neg(nmod_div(ctx->modulus->coeffs[0],
              ctx->modulus->coeffs[1], ctx->mod), ctx->mod));
    }
    else
    {
        nmod_poly_zero(rop);
        nmod_poly_set_coeff_ui(rop, 0, 0);
        nmod_poly_set_coeff_ui(rop, 1, 1);
    }
}
