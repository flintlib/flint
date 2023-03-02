/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2017 Claus Fieker


    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"
#include "nmod_mat.h"

/*
    This computes the norm on $\mathbf{F}_q$.
 */

void _fq_nmod_norm(fmpz_t rop2, const mp_limb_t *op, slong len, 
                        const fq_nmod_ctx_t ctx)
{
    const slong d = fq_nmod_ctx_degree(ctx);

    mp_limb_t rop;

    if (d == 1)
    {
        rop = op[0];
    } 
    else if (len == 1) /* element scalar */
    {
        rop = n_powmod2_ui_preinv(op[0], d, ctx->mod.n, ctx->mod.ninv);
    }
    else 
    {
        rop = _nmod_poly_resultant(ctx->modulus->coeffs, ctx->modulus->length,
            op, len, ctx->mod);

        /*
            XXX:  This part of the code is currently untested as the Conway 
            polynomials used for the extension Fq/Fp are monic.

            TODO: make polynomial monic!!!
              reading the source, a[] is monic, modulus is not. Why????
         */
        if (ctx->modulus->coeffs[d] != WORD(1))
        {
            mp_limb_t f;
            f = n_powmod2_ui_preinv(ctx->modulus->coeffs[d], len - 1, ctx->mod.n, ctx->mod.ninv);
            f = n_invmod(f, ctx->mod.n);
            rop = n_mulmod2_preinv(f, rop, ctx->mod.n, ctx->mod.ninv);
        }
    }

    fmpz_set_ui(rop2, rop);

}

void fq_nmod_norm(fmpz_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
{
    if (fq_nmod_is_zero(op, ctx))
    {
        fmpz_zero(rop);
        return;
    }
    else
    {
        _fq_nmod_norm(rop, op->coeffs, op->length, ctx); 
    }
}

