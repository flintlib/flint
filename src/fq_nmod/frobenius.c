/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fq_nmod.h"

/*
    Sets (rop, 2d-1) to the image of (op, len) under the Frobenius operator 
    raised to the e-th power, assuming that neither op nor e are zero.
 */

void _fq_nmod_frobenius(mp_limb_t *rop, const mp_limb_t *op, slong len, slong e, 
                        const fq_nmod_ctx_t ctx)
{
    const slong d = ctx->j[ctx->len - 1];

    if (len == 1)  /* op is in Fp, not just Fq */
    {
        _nmod_vec_set(rop, op, len);
        _nmod_vec_zero(rop + len, (2*d - 1)  - len);
    }
    else
    {
        fmpz_t t;

        fmpz_init(t);
        fmpz_pow_ui(t, fq_nmod_ctx_prime(ctx), e);
        _fq_nmod_pow(rop, op, len, t, ctx);
        fmpz_clear(t);
    }
}

void fq_nmod_frobenius(fq_nmod_t rop, const fq_nmod_t op, slong e, const fq_nmod_ctx_t ctx)
{
    const slong d = fq_nmod_ctx_degree(ctx);

    e = e % d;
    if (e < 0)
        e += d;

    if (fq_nmod_is_zero(op, ctx))
    {
        fq_nmod_zero(rop, ctx);
    }
    else if (e == 0)
    {
        fq_nmod_set(rop, op, ctx);
    }
    else
    {
        mp_limb_t *t;

        if (rop == op)
        {
            t = _nmod_vec_init(2 * d - 1);
        }
        else
        {
            nmod_poly_fit_length(rop, 2 * d - 1);
            t = rop->coeffs;
        }

        _fq_nmod_frobenius(t, op->coeffs, op->length, e, ctx);

        if (rop == op)
        {
            _nmod_vec_clear(rop->coeffs);
            rop->coeffs = t;
            rop->alloc  = 2 * d - 1;
            rop->length = d;
        }
        else
        {
            _nmod_poly_set_length(rop, d);
        }
        _nmod_poly_normalise(rop);
    }
}

