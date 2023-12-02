/*
    Copyright (C) 2011, 2012 Sebastian Pancratz
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"
#include "fq.h"

void
_fq_inv(fmpz * rop, const fmpz * op, slong len, const fq_ctx_t ctx)
{
    const slong d = fq_ctx_degree(ctx);

    if (len == 1)
    {
        fmpz_invmod(rop, op, fq_ctx_prime(ctx));
        _fmpz_vec_zero(rop + 1, d - 1);
    }
    else
    {
        _fmpz_mod_poly_invmod(rop, op, len, ctx->modulus->coeffs, d + 1, ctx->ctxp);
    }
}

void
fq_inv(fq_t rop, const fq_t op, const fq_ctx_t ctx)
{
    if (fq_is_zero(op, ctx))
    {
        flint_throw(FLINT_ERROR, "Exception (fq_inv).  Zero is not invertible.\n");
    }
    else
    {
        const slong d = fq_ctx_degree(ctx);
        fmpz *t;

        if (rop == op)
        {
            t = _fmpz_vec_init(d);
        }
        else
        {
            fmpz_poly_fit_length(rop, d);
            t = rop->coeffs;
        }

        _fq_inv(t, op->coeffs, op->length, ctx);

        if (rop == op)
        {
            _fmpz_vec_clear(rop->coeffs, rop->alloc);
            rop->coeffs = t;
            rop->alloc = d;
            rop->length = d;
        }
        else
        {
            _fmpz_poly_set_length(rop, d);
        }
        _fmpz_poly_normalise(rop);
    }
}
