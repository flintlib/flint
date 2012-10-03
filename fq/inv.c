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
 
******************************************************************************/

#include "fq.h"

void _fq_inv(fmpz *rop, const fmpz *op, long len, 
             const fmpz *a, const long *j, long lena, const fmpz_t p)
{
    const long d = j[lena - 1];

    if (len == 1)
    {
        fmpz_invmod(rop, op, p);
        _fmpz_vec_zero(rop + 1, d - 1);
    }
    else
    {
        fmpz *f = _fmpz_vec_init(d + 1);
        long k;

        for (k = 0; k < lena; k++)
            fmpz_set(f + j[k], a + k);

        _fmpz_mod_poly_invmod(rop, op, len, f, d + 1, p);

        _fmpz_vec_clear(f, d + 1);
    }
}

void fq_inv(fq_t rop, const fq_t op, const fq_ctx_t ctx)
{
    if (fq_is_zero(op))
    {


        printf("Exception (fq_inv).  Zero is not invertible.\n");
        abort();
    }
    else
    {
        const long d = fq_ctx_degree(ctx);
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

        _fq_inv(t, op->coeffs, op->length, 
                ctx->a, ctx->j, ctx->len, fq_ctx_prime(ctx));

        if (rop == op)
        {
            _fmpz_vec_clear(rop->coeffs, rop->alloc);
            rop->coeffs = t;
            rop->alloc  = d;
            rop->length = d;
        }
        else
        {
            _fmpz_poly_set_length(rop, d);
        }
        _fmpz_poly_normalise(rop);
    }
}
