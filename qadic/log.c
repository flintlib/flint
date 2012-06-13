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

void _qadic_log(fmpz *z, const fmpz *y, long v, long len, 
                const fmpz *a, const long *j, long lena, 
                const fmpz_t p, long N, const fmpz_t pN)
{
    const long d = j[lena - 1];

    if (N < (1L < 10) / (long) fmpz_bits(p))
    {
        _qadic_log_rectangular(z, y, v, len, a, j, lena, p, N, pN);
    }
    else
    {
        _qadic_log_balanced(z, y, len, a, j, lena, p, N, pN);
    }
}

int qadic_log(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
{
    const fmpz *p  = (&ctx->pctx)->p;
    const long d   = qadic_ctx_degree(ctx);
    const long N   = (&ctx->pctx)->N;
    const long len = op->length;

    if (op->val < 0)
    {
        return 0;
    }
    else
    {
        fmpz *x;
        fmpz_t pN;
        int alloc, ans;

        x = _fmpz_vec_init(len + 1);
        alloc = _padic_ctx_pow_ui(pN, N, &ctx->pctx);

        /* Set x := (1 - op) mod p^N */
        fmpz_pow_ui(x + len, p, op->val);
        _fmpz_vec_scalar_mul_fmpz(x, op->coeffs, len, x + len);
        fmpz_sub_ui(x, x, 1);
        _fmpz_vec_neg(x, x, len);
        _fmpz_vec_scalar_mod_fmpz(x, x, len, pN);

        if (_fmpz_vec_is_zero(x, len))
        {
            padic_poly_zero(rop);
            ans = 1;
        }
        else
        {
            const long v = _fmpz_vec_ord_p(x, len, p);

            if (v >= 2 || (*p != 2L && v >= 1))
            {
                if (v >= N)
                {
                    padic_poly_zero(rop);
                }
                else
                {
                    padic_poly_fit_length(rop, d);

                    _qadic_log(rop->coeffs, x, v, len, 
                               ctx->a, ctx->j, ctx->len, p, N, pN);
                    rop->val = 0;

                    _padic_poly_set_length(rop, d);
                    _padic_poly_normalise(rop);
                    padic_poly_canonicalise(rop, p);
                }
                ans = 1;
            }
            else
            {
                ans = 0;
            }
        }

        _fmpz_vec_clear(x, len + 1);
        if (alloc)
            fmpz_clear(pN);
        return ans;
    }
}

