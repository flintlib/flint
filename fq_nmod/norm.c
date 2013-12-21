/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include "fq_nmod.h"
#include "nmod_mat.h"

/*
    Computes the norm on $\mathbf{Q}_q$ to precision $N \geq 1$. 
    When $N = 1$, this computes the norm on $\mathbf{F}_q$.
 */

void _fq_nmod_norm(fmpz_t rop2, const mp_limb_t *op, slong len, 
                        const fq_nmod_ctx_t ctx)
{
    const slong d = fq_nmod_ctx_degree(ctx);

    mp_limb_t rop;

    if (len == 1)
    {
        rop = n_powmod2_ui_preinv(op[0], d, ctx->mod.n, ctx->mod.ninv);
    }
    else 
    {
        /* Construct an ad hoc matrix M and set rop to det(M) */
        {
            const slong n = d + len - 1;
            slong i, k;
            nmod_mat_t M;

            nmod_mat_init(M, n, n, ctx->mod.n);

            for (k = 0; k < len-1; k++)
            {
                for (i = 0; i < ctx->len; i++)
                {
                    M->entries[k * n + k + (d - ctx->j[i])] = ctx->a[i];
                }
            }
            for (k = 0; k < d; k++)
            {
                for (i = 0; i < len; i++)
                {
                    M->entries[(len-1 + k) * n + k + (len-1 - i)] = op[i];
                }
            }

            rop = _nmod_mat_det(M);

            nmod_mat_clear(M);
        }

        /*
            XXX:  This part of the code is currently untested as the Conway 
            polynomials used for the extension Fq/Fp are monic.
         */
        if (ctx->a[ctx->len - 1] != WORD(1))
        {
            mp_limb_t f;
            f = n_powmod2_ui_preinv(ctx->a[ctx->len - 1], len - 1, ctx->mod.n, ctx->mod.ninv);
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

