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

/*
    Discussion on the choice of the norm algorithm.

    When the logarithm function does not converge for x, 
    the only choice is the resultant method.

    However, when the logarithm function converges, we 
    can choose between the analytic method and the resultant 
    method.  Roughly speaking, we postulate that the analytic 
    method has runtime A (log N)^2 mu(p,d,N), where mu(p,d,N) 
    is (d log d) M(N log p).  The resultant method has runtime 
    B d^4 M(N log p).  Experimentally, we find that A/B is 
    somewhere around 4.

    TODO:  Repeat the experiments with p=2, which is an 
    important special case.
 */

void _qadic_norm(fmpz_t rop, const fmpz *op, long len, 
                 const fmpz *a, const long *j, long lena, 
                 const fmpz_t p, long N)
{
    const long d = j[lena - 1];

    if (len == 1)
    {
        fmpz_t pN;

        fmpz_init(pN);
        fmpz_pow_ui(pN, p, N);
        fmpz_powm_ui(rop, op + 0, d, pN);
        fmpz_clear(pN);
    }
    else
    {
        fmpz *y;
        long w;

        y = _fmpz_vec_init(len);

        /* (y,len) := 1 - (op,len) */
        _fmpz_vec_neg(y, op, len);
        fmpz_add_ui(y + 0, y + 0, 1);

        w = _fmpz_vec_ord_p(y, len, p);

        if (w >= 2 || (*p != 2L && w >= 1))
        {
            if (4 * FLINT_FLOG2(N) * FLINT_FLOG2(N) * FLINT_FLOG2(d) < d*d*d)
            {
                _qadic_norm_analytic(rop, y, w, len, a, j, lena, p, N);
            }
            else
            {
                _qadic_norm_resultant(rop, op, len, a, j, lena, p, N);
            }
        }
        else
        {
            _qadic_norm_resultant(rop, op, len, a, j, lena, p, N);
        }

        _fmpz_vec_clear(y, len);
    }
}

void qadic_norm(padic_t rop, const qadic_t op, const qadic_ctx_t ctx)
{
    const long N  = (&ctx->pctx)->N;
    const long d  = qadic_ctx_degree(ctx);
    const fmpz *p = (&ctx->pctx)->p;

    /* N(p^v u) = p^{dv} N(u) */

    if (qadic_is_zero(op) || d * op->val >= N)
    {
        padic_zero(rop);
    }
    else
    {
        _qadic_norm(padic_unit(rop), op->coeffs, op->length, 
                    ctx->a, ctx->j, ctx->len, p, N - d * op->val);
        padic_val(rop) = d * op->val;
    }
}

