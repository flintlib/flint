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

#include "qadic.h"

void _qadic_trace(fmpz_t rop, const fmpz *op, long len, 
                  const fmpz *a, const long *j, long lena, const fmpz_t pN)
{
    const long d = j[lena - 1];

    long i;
    fmpz *t;

    t  = _fmpz_vec_init(2 * d - 1);

    fmpz_set(rop, op);
    for (i = 1; i < d; i++)
    {
        _fmpz_vec_zero(t, i);
        _fmpz_vec_set(t + i, op, len);
        _fmpz_mod_poly_reduce(t, len + i, a, j, lena, pN);
        fmpz_add(rop, rop, t + i);
    }
    fmpz_mod(rop, rop, pN);

    _fmpz_vec_clear(t, 2 * d - 1);
}

void qadic_trace(padic_t rop, const qadic_t op, const qadic_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;

    if (qadic_is_zero(op) || op->val >= N)
    {
        padic_zero(rop);
    }
    else
    {
        fmpz_t pN;
        int alloc;

        alloc = _padic_ctx_pow_ui(pN, N - op->val, &ctx->pctx);

        _qadic_trace(padic_unit(rop), op->coeffs, op->length, 
                     ctx->a, ctx->j, ctx->len, pN);
        padic_val(rop) = op->val;

        _padic_canonicalise(rop, (&ctx->pctx)->p);

        if (alloc)
            fmpz_clear(pN);
    }
}

