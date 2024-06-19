/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2017 Claus Fieker

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "fq.h"

/*
    Computes the norm on $\mathbf{F}_q$.
 */

void
_fq_norm(fmpz_t rop, const fmpz * op, slong len, const fq_ctx_t ctx)
{
    const slong d = fq_ctx_degree(ctx);

    if (d == 1)
    {
        fmpz_set(rop, op + 0);
    }
    else if (len == 1)
    {
        fmpz_mod_pow_ui(rop, op + 0, d, ctx->ctxp);
    }
    else
    {
        _fmpz_mod_poly_resultant(rop, ctx->modulus->coeffs, ctx->modulus->length,
                        op, len, ctx->ctxp);
        /*
           XXX:  This part of the code is currently untested as the Conway
           polynomials used for the extension Fq/Fp are monic.

           TODO: make polynomial monic!!!
           reading the source, a[] is monic, modulus is not. Why????
         */
        if (!fmpz_is_one(ctx->modulus->coeffs + d))
        {
            fmpz_t f;

            fmpz_init(f);
            fmpz_mod_pow_ui(f, ctx->modulus->coeffs + d, len - 1, ctx->ctxp);
            fmpz_mod_inv(f, f, ctx->ctxp);
            fmpz_mod_mul(rop, f, rop, ctx->ctxp);
            fmpz_clear(f);
        }
    }
}

void
fq_norm(fmpz_t rop, const fq_t op, const fq_ctx_t ctx)
{
    if (fq_is_zero(op, ctx))
    {
        fmpz_zero(rop);
    }
    else
    {
        _fq_norm(rop, op->coeffs, op->length, ctx);
    }
}
