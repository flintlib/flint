/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "qadic.h"

/*
    Uses Hensel lifting along the polynomial $X^q - X$, which yields
    the formula $z' = z - (z^q - z) / (q z^{q-1} - 1)$.

    We observe that the denominator is an approximation to $q - 1$,
    which allows us to use the formula $z' = z - (q-1)^{-1} (z^q - z)$
    during the iteration.

    Supports aliasing between \code{rop} and \code{op}.
 */

void _qadic_teichmuller(fmpz *rop, const fmpz *op, slong len,
                        const fmpz *a, const slong *j, slong lena,
                        const fmpz_t p, slong N)
{
    const slong d = j[lena - 1];

    if (len == 1)
    {
        _padic_teichmuller(rop, op, p, N);
        _fmpz_vec_zero(rop + 1, d - 1);
    }
    else if (N == 1)
    {
        _fmpz_vec_scalar_mod_fmpz(rop, op, len, p);
        _fmpz_vec_zero(rop + len, d - len);
    }
    else  /* d, N >= 2 */
    {
        slong *e, i, n;
        fmpz *pow, *u, *t, *w;
        fmpz_t inv, q, qm1;

        n = FLINT_CLOG2(N) + 1;

        e = flint_malloc(n * sizeof(slong));
        for (e[i = 0] = N; e[i] > 1; i++)
            e[i + 1] = (e[i] + 1) / 2;

        w   = _fmpz_vec_init(n + n + (2 * d - 1));
        pow = w;
        u   = w + n;
        t   = w + 2 * n;

        fmpz_init(inv);
        fmpz_init(q);
        fmpz_init(qm1);

        fmpz_pow_ui(q, p, d);
        fmpz_sub_ui(qm1, q, 1);

        /* Compute powers of p */
        {
            fmpz_one(t);
            fmpz_set(pow + i, p);
        }
        for (i--; i >= 1; i--)
        {
            if (e[i] & WORD(1))
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
            if (e[i] & WORD(1))
                fmpz_mul(pow + i, t, pow + (i + 1));
            else
                fmpz_mul(pow + i, pow + (i + 1), pow + (i + 1));
        }

        /* Compute reduced units for (q-1) */
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
            _fmpz_vec_scalar_mod_fmpz(rop, op, len, pow + i);
            _fmpz_vec_zero(rop + len, d - len);

            fmpz_sub_ui(inv, p, 1);
        }
        for (i--; i >= 0; i--)
        {
            /* Lift rop */
            _qadic_pow(t, rop, d, q, a, j, lena, pow + i);
            _fmpz_poly_sub(t, t, d, rop, d);
            _fmpz_vec_scalar_submul_fmpz(rop, t, d, inv);
            _fmpz_vec_scalar_mod_fmpz(rop, rop, d, pow + i);

            /* Lift inv */
            if (i > 0)
            {
                fmpz_mul(t, inv, inv);
                fmpz_mul(t + 1, u + i, t);
                fmpz_mul_2exp(inv, inv, 1);
                fmpz_sub(inv, inv, t + 1);
                fmpz_mod(inv, inv, pow + i);
            }
        }

        _fmpz_vec_clear(w, n + n + (2 * d - 1));
        fmpz_clear(inv);
        fmpz_clear(q);
        fmpz_clear(qm1);
        flint_free(e);
    }
}

void qadic_teichmuller(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
{
    const slong N = qadic_prec(rop);

    if (op->val < 0)
    {
        flint_throw(FLINT_ERROR, "Exception (qadic_teichmuller).  val(op) is negative.\n");
    }

    if (qadic_is_zero(op) || op->val > 0 || N <= 0)
    {
        qadic_zero(rop);
    }
    else
    {
        const slong d = qadic_ctx_degree(ctx);

        padic_poly_fit_length(rop, d);

        _qadic_teichmuller(rop->coeffs, op->coeffs, op->length,
                           ctx->a, ctx->j, ctx->len, (&ctx->pctx)->p, N);
        rop->val = 0;
        _padic_poly_set_length(rop, d);
        _padic_poly_normalise(rop);
    }
}

