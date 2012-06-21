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

    long i, l;
    fmpz *t;

    t = _fmpz_vec_init(d);

    fmpz_set_ui(t + 0, d);
    for (i = 1; i < d; i++)
    {
        for (l = lena - 2; l >= 0 && j[l] >= d - (i - 1); l--)
        {
            fmpz_addmul(t + i, t + (j[l] + i - d), a + l);
        }

        if (l >= 0 && j[l] == d - i)
        {
            fmpz_addmul_ui(t + i, a + l, i);
        }

        fmpz_neg(t + i, t + i);
        fmpz_mod(t + i, t + i, pN);
    }

    fmpz_zero(rop);
    for (i = 0; i < d; i++)
    {
        fmpz_addmul(rop, op + i, t + i);
    }
    fmpz_mod(rop, rop, pN);

    _fmpz_vec_clear(t, d);
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

        _padic_canonicalise(rop, &ctx->pctx);

        if (alloc)
            fmpz_clear(pN);
    }
}

