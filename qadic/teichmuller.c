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

#include <stdlib.h>
#include "fmpz_mod_poly.h"
#include "qadic.h"

void _qadic_teichmuller(fmpz *rop, const fmpz *op, long len, 
                        const fmpz *a, const long *j, long lena, 
                        const fmpz_t p, long N)
{
    const long d = j[lena - 1];

    if (d == 1)
    {
        _padic_teichmuller(rop, op, p, N);
    }
    else if (N == 1)
    {
        _fmpz_vec_scalar_mod_fmpz(rop, op, len, p);
        _fmpz_vec_zero(rop + len, d - len);
    }
    else
    {
        long *e, i, n;
        fmpz *pow, *inv, *r, *s, *t;
        fmpz_t qm1, qm2;

        n = FLINT_CLOG2(N) + 1;

        e = flint_malloc(n * sizeof(long));
        for (e[i = 0] = N; e[i] > 1; i++)
            e[i + 1] = (e[i] + 1) / 2;

        pow = _fmpz_vec_init(n);
        inv = _fmpz_vec_init(2 * d - 1);
        r   = _fmpz_vec_init(2 * d - 1);
        s   = _fmpz_vec_init(2 * d - 1);
        t   = _fmpz_vec_init(2 * d - 1);

        fmpz_init(qm1);
        fmpz_init(qm2);

        fmpz_pow_ui(qm1, p, d);
        fmpz_sub_ui(qm1, qm1, 1);
        fmpz_sub_ui(qm2, qm1, 1);

        /* Compute powers of p */
        {
            fmpz_one(t);
            fmpz_set(pow + i, p);
        }
        for (i--; i >= 1; i--)
        {
            if (e[i] & 1L)
            {
                fmpz_mul(pow + i, t, pow + (i + 1));
                fmpz_mul(t, t, t);
            }
            else
            {
                fmpz_mul(t, t, pow + (i + 1));
                fmpz_mul(pow + i, pow + (i + 1), pow + (i + 1));
            }
        }
        {
            if (e[i] & 1L)
                fmpz_mul(pow + i, t, pow + (i + 1));
            else
                fmpz_mul(pow + i, pow + (i + 1), pow + (i + 1));
        }

        /* Run Newton iteration */
        i = n - 1;
        {
            _fmpz_vec_scalar_mod_fmpz(rop, op, len, pow + i);
            _fmpz_vec_zero(rop + len, d - len);

            /* {(q-1) x^{q-2}}^{-1} */
            _fmpz_mod_poly_neg(inv, rop, d, pow + i);
        }
        for (i--; i >= 0; i--)
        {
            /* Lift rop */
            _qadic_pow(t, rop, d, qm1, a, j, lena, pow + i);
            fmpz_sub_ui(t, t, 1);
            _fmpz_mod_poly_mul(s, t, d, inv, d, pow + i);
            _fmpz_mod_poly_reduce(s, 2*d - 1, a, j, lena, pow + i);
            _fmpz_mod_poly_sub(rop, rop, d, s, d, pow + i);
            
            /* Lift inv */
            if (i > 0)
            {
                _qadic_pow(s, rop, d, qm2, a, j, lena, pow + i);
                _fmpz_mod_poly_scalar_mul_fmpz(t, inv, d, qm1, pow + i);
                _fmpz_mod_poly_mul(r, s, d, t, d, pow + i);
                _fmpz_mod_poly_reduce(r, 2*d - 1, a, j, lena, pow + i);
                fmpz_sub_ui(r, r, 2);
                _fmpz_mod_poly_neg(r, r, d, pow + i);
                _fmpz_mod_poly_mul(s, inv, d, r, d, pow + i);
                _fmpz_mod_poly_reduce(s, 2*d - 1, a, j, lena, pow + i);
                {
                    fmpz *__t;
                    __t = s; s = inv; inv = __t;
                }
            }
        }

        _fmpz_vec_clear(pow, n);
        _fmpz_vec_clear(inv, 2 * d - 1);
        _fmpz_vec_clear(r, 2 * d - 1);
        _fmpz_vec_clear(s, 2 * d - 1);
        _fmpz_vec_clear(t, 2 * d - 1);
        fmpz_clear(qm1);
        fmpz_clear(qm2);
        flint_free(e);
    }
}

void qadic_teichmuller(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;

    if (op->val < 0)
    {
        printf("Exception (qadic_teichmuller).  val(op) is negative.\n");
        abort();
    }

    if (qadic_is_zero(op) || op->val > 0 || N <= 0)
    {
        qadic_zero(rop);
    }
    else
    {
        const long d = qadic_ctx_degree(ctx);

        padic_poly_fit_length(rop, d);

        _qadic_teichmuller(rop->coeffs, op->coeffs, op->length, 
                           ctx->a, ctx->j, ctx->len, (&ctx->pctx)->p, N);
        rop->val = 0;
        _padic_poly_set_length(rop, d);
        _padic_poly_normalise(rop);
    }
}

