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

/*
    Uses Hensel lifting along the polynomial $X^q - X$, which yields 
    the formula $z' = z - (z^q - z) / (q z^{q-1} - 1)$.

    We observe that the denominator is $p$-adically close to $q - 1$ 
    since $z^{q-1}$ is close to~$1$.  This allows us to use the formula 
    $z = z + (1-q)^{-1} (z^q - z)$, during the iteration where the 
    $p$-adic inverse of $(1-q)$ is updated at each step, too.

    Supports aliasing between \code{rop} and \code{op}.
 */

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
        fmpz *pow, *u, *t, *w;
        fmpz_t inv, q, qm1;

        n = FLINT_CLOG2(N) + 1;

        e = flint_malloc(n * sizeof(long));
        for (e[i = 0] = N; e[i] > 1; i++)
            e[i + 1] = (e[i] + 1) / 2;

        w   = _fmpz_vec_init(n + n + 2 * d - 1);
        pow = w;
        u   = w + n;
        t   = w + 2 * n;

        fmpz_init(inv);
        fmpz_init(q);
        fmpz_init(qm1);

        fmpz_pow_ui(q, p, d);
        fmpz_sub_ui(qm1, q, 1);
        fmpz_neg(qm1, qm1);

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

        /* Compute reduced units for (1-q) */
        {
            fmpz_mod(u + 0, qm1, pow + 0);
        }
        for (i = 1; i < n; i++)
        {
            fmpz_mod(u + i, u + (i - 1), pow + i);
        }

        /* Run Newton iteration */
        i = n - 1;
        {
            fmpz_one(inv);

            _fmpz_vec_scalar_mod_fmpz(rop, op, len, pow + i);
            _fmpz_vec_zero(rop + len, d - len);
        }
        for (i--; i >= 0; i--)
        {
            /* Lift u := (1-q)^{-1} */
            fmpz_mul(t, u + i, inv);
            fmpz_sub_ui(t, t, 1);
            fmpz_mul(t + 1, t, inv);
            fmpz_sub(inv, inv, t + 1);
            fmpz_mod(inv, inv, pow + i);

            /* Lift rop */
            _qadic_pow(t, rop, d, q, a, j, lena, pow + i);
            _fmpz_mod_poly_sub(t, t, d, rop, d, pow + i);
            _fmpz_mod_poly_scalar_mul_fmpz(t, t, d, inv, pow + i);
            _fmpz_mod_poly_add(rop, rop, d, t, d, pow + i);
        }

        _fmpz_vec_clear(w, n + n + 2 * d - 1);
        fmpz_clear(inv);
        fmpz_clear(q);
        fmpz_clear(qm1);
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

