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

    Copyright (C) 2011, 2012 Sebastian Pancratz
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include "fq_nmod.h"

void _fq_nmod_inv(mp_limb_t *rop, const mp_limb_t *op, slong len,
                  const fq_nmod_ctx_t ctx)
{
    const slong d = fq_nmod_ctx_degree(ctx);

    if (len == 1)
    {
        rop[0] = n_invmod(op[0], ctx->mod.n);
        _nmod_vec_zero(rop + 1, d - 1);
    }
    else
    {
        _nmod_poly_invmod(rop, op, len, ctx->modulus->coeffs, d + 1, ctx->mod);
    }
}

void fq_nmod_inv(fq_nmod_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
{
    if (fq_nmod_is_zero(op, ctx))
    {
        flint_printf("Exception (fq_nmod_inv).  Zero is not invertible.\n");
        abort();
    }
    else
    {
        const slong d = fq_nmod_ctx_degree(ctx);
        mp_limb_t *t;

        if (rop == op)
        {
            t = _nmod_vec_init(d);
        }
        else
        {
            nmod_poly_fit_length(rop, d);
            t = rop->coeffs;
        }

        _fq_nmod_inv(t, op->coeffs, op->length, ctx);

        if (rop == op)
        {
            _nmod_vec_clear(rop->coeffs);
            rop->coeffs = t;
            rop->alloc  = d;
            rop->length = d;
        }
        else
        {
            _nmod_poly_set_length(rop, d);
        }
        _nmod_poly_normalise(rop);
    }
}
