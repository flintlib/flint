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
 
******************************************************************************/

#include "fmpz_mod_poly.h"
#include "qadic.h"

void _qadic_exp(fmpz *rop, const fmpz *op, long v, long len, 
                           const fmpz *a, const long *j, long lena, 
                           const fmpz_t p, long N, const fmpz_t pN)
{
    if (N < (1L << 13) / (long) fmpz_bits(p))
    {
        _qadic_exp_rectangular(rop, op, v, len, a, j, lena, p, N, pN);
    }
    else
    {
        const long d = j[lena - 1];

        _qadic_exp_balanced(rop, op, v, len, a, j, lena, p, N, pN);
        _fmpz_vec_zero(rop + d, d - 1);
    }
}

int qadic_exp(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
{
    const long N  = (&ctx->pctx)->N;
    const long v  = op->val;
    const fmpz *p = (&ctx->pctx)->p;

    if (padic_poly_is_zero(op))
    {
        padic_poly_one(rop, &ctx->pctx);
        return 1;
    }

    if ((*p == 2L && v <= 1) || (v <= 0))
    {
        return 0;
    }
    else
    {
        if (v < N)
        {
            const long d = qadic_ctx_degree(ctx);

            fmpz *t;
            fmpz_t pN;
            int alloc;

            alloc = _padic_ctx_pow_ui(pN, N, &ctx->pctx);

            if (rop != op)
            {
                padic_poly_fit_length(rop, 2 * d - 1);
                t = rop->coeffs;
            }
            else
            {
                t = _fmpz_vec_init(2 * d - 1);
            }

            _qadic_exp(t, op->coeffs, v, op->length, 
                       ctx->a, ctx->j, ctx->len, p, N, pN);
            rop->val = 0;

            if (rop == op)
            {
                _fmpz_vec_clear(rop->coeffs, rop->alloc);
                rop->coeffs = t;
                rop->alloc  = 2 * d - 1;
                rop->length = d;
            }
            _padic_poly_set_length(rop, d);
            _padic_poly_normalise(rop);

            if (alloc)
                fmpz_clear(pN);
        }
        else
        {
            padic_poly_one(rop, &ctx->pctx);
        }
        return 1;
    }
}

