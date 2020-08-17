/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq.h"

void
fq_gcdinv(fq_t rop, fq_t inv, const fq_t op, const fq_ctx_t ctx)
{
    fmpz *g, *s;
    slong lenG;
    fmpz_t f;

    if (fq_is_zero(op, ctx))
    {
        fq_zero(rop, ctx);
        return;
    }

    fmpz_init(f);

    if (rop == op)
    {
        g = _fmpz_vec_init(op->length);
    }
    else
    {
        fmpz_poly_fit_length(rop, op->length);
        g = rop->coeffs;
    }

    if (inv == op)
    {
        s = _fmpz_vec_init(ctx->modulus->length - 1);
    }
    else
    {
        fmpz_poly_fit_length(inv, ctx->modulus->length - 1);
        s = inv->coeffs;
    }

    lenG = _fmpz_mod_poly_gcdinv_f(f, g, s, op->coeffs, op->length,
                                 ctx->modulus->coeffs, ctx->modulus->length,
                                 fq_ctx_prime(ctx));
    if (rop == op)
    {
        _fmpz_vec_clear(rop->coeffs, rop->alloc);
        rop->coeffs = g;
        rop->alloc = op->length;
    }
    if (inv == op)
    {
        _fmpz_vec_clear(inv->coeffs, inv->alloc);
        inv->coeffs = s;
        inv->alloc = ctx->modulus->length - 1;
    }

    if (!fmpz_is_one(f))
    {
        _fmpz_poly_set_length(inv, 0);
        fq_zero(rop, ctx);
        goto cleanup;
    }

    _fmpz_poly_set_length(rop, lenG);
    _fmpz_poly_set_length(inv, ctx->modulus->length - lenG);
    _fmpz_poly_normalise(inv);

    if (fmpz_is_one(fmpz_poly_lead(rop)))
        goto cleanup;

    if (!fmpz_invmod(f, fmpz_poly_lead(rop), fq_ctx_prime(ctx)))
    {
        fq_zero(rop, ctx);
        goto cleanup;
    }
    _fmpz_mod_poly_scalar_mul_fmpz(rop->coeffs, rop->coeffs, rop->length,
                                   f, fq_ctx_prime(ctx));
    _fmpz_mod_poly_scalar_mul_fmpz(inv->coeffs, inv->coeffs, inv->length,
                                   f, fq_ctx_prime(ctx));
cleanup:

    fmpz_clear(f);
    return;
}

