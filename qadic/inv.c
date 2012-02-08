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

#include <stdio.h>
#include <string.h>

#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"
#include "padic.h"
#include "qadic.h"

/*
    N >= 1.  Sets (rop, d).  Does not support aliasing.
 */

void _qadic_inv(fmpz *rop, const fmpz *op, long len, 
                const fmpz *a, const long *j, long lena, 
                const fmpz_t p, long N)
{
    const long d = j[lena - 1];

    if (d == 1)
    {
        _padic_inv(rop, op, p, N);
    }
    else if (N == 1)
    {
        fmpz *P = _fmpz_vec_init(d + 1);
        long k;

        for (k = 0; k < lena; k++)
            fmpz_set(P + j[k], a + k);

        _fmpz_mod_poly_invmod(rop, op, len, P, d + 1, p);

        _fmpz_vec_clear(P, d + 1);
    }
    else  /* d, N >= 2 */
    {
        long *e, i, n;
        fmpz *pow, *u;
        fmpz *s, *t;

        n = FLINT_CLOG2(N) + 1;

        /* Compute sequence of exponents */
        e = flint_malloc(n * sizeof(long));
        for (e[i = 0] = N; e[i] > 1; i++)
            e[i + 1] = (e[i] + 1) / 2;

        pow = _fmpz_vec_init(n);
        u   = _fmpz_vec_init(len * n);
        s   = _fmpz_vec_init(2 * d - 1);
        t   = _fmpz_vec_init(2 * d - 1);

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

        /* Compute reduced units */
        {
            _fmpz_vec_scalar_mod_fmpz(u + 0 * len, op, len, pow + 0);
        }
        for (i = 1; i < n; i++)
        {
            _fmpz_vec_scalar_mod_fmpz(u + i * len, u + (i - 1) * len, len, pow + i);
        }

        /* Run Newton iteration */
        i = n - 1;
        {
            fmpz *P = _fmpz_vec_init(d + 1);
            long k;

            for (k = 0; k < lena; k++)
                fmpz_set(P + j[k], a + k);

            _fmpz_mod_poly_invmod(rop, u + i * len, len, P, d + 1, pow + i);

            _fmpz_vec_clear(P, d + 1);
        }
        for (i--; i >= 0; i--)
        {
            _fmpz_mod_poly_mul(s, rop, d, u + i * len, len, pow + i);
            _fmpz_mod_poly_reduce(s, d + len - 1, a, j, lena, pow + i);

            fmpz_sub_ui(s, s, 1);

            _fmpz_mod_poly_mul(t, s, d, rop, d, pow + i);
            _fmpz_mod_poly_reduce(t, 2*d - 1, a, j, lena, pow + i);

            _fmpz_mod_poly_sub(rop, rop, d, t, d, pow + i);
        }

        _fmpz_vec_clear(pow, n);
        _fmpz_vec_clear(u, len * n);
        _fmpz_vec_clear(s, 2 * d - 1);
        _fmpz_vec_clear(t, 2 * d - 1);
        flint_free(e);
    }
}

void qadic_inv(qadic_t x, const qadic_t y, const qadic_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;

    if (qadic_is_zero(y))
    {
        printf("Exception (qadic_inv).  Zero is not invertible.\n");
        abort();
    }

    /*
        If y = u p^v has negative valuation with N <= -v then the 
        exact inverse of y is zero when reduced modulo $p^N$
     */
    if (N + y->val <= 0)
    {
        qadic_zero(x);
    }
    else
    {
        const long d = qadic_ctx_degree(ctx);
        fmpz *t;

        if (x == y)
        {
            t = _fmpz_vec_init(d);
        }
        else
        {
            padic_poly_fit_length(x, d);
            t = x->coeffs;
        }

        _qadic_inv(t, y->coeffs, y->length, 
                   ctx->a, ctx->j, ctx->len, (&ctx->pctx)->p, N + y->val);
        x->val = - y->val;

        if (x == y)
        {
            _fmpz_vec_clear(x->coeffs, x->alloc);
            x->coeffs = t;
            x->alloc  = d;
            x->length = d;
        }
        else
        {
            _padic_poly_set_length(x, d);
        }
        _padic_poly_normalise(x);
    }
}

