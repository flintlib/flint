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

#include "padic.h"

void _padic_teichmuller(fmpz_t rop, const fmpz_t op, const fmpz_t p, long N)
{
    if (*p == 2L)
    {
        fmpz_one(rop);
    }
    else if (N == 1)
    {
        fmpz_mod(rop, op, p);
    }
    else
    {
        long *a, i, n;
        fmpz *pow;
        fmpz_t r, s, t;
        fmpz_t inv;
        fmpz_t pm1, pm2;

        n = FLINT_CLOG2(N) + 1;

        a = malloc(n * sizeof(long));
        for (a[i = 0] = N; a[i] > 1; i++)
            a[i + 1] = (a[i] + 1) / 2;

        pow = _fmpz_vec_init(n);

        fmpz_init(r);
        fmpz_init(s);
        fmpz_init(t);
        fmpz_init(inv);
        fmpz_init(pm1);
        fmpz_init(pm2);

        fmpz_sub_ui(pm1, p, 1);
        fmpz_sub_ui(pm2, p, 2);

        /* Compute powers of p */
        {
            fmpz_one(t);
            fmpz_set(pow + i, p);
        }
        for (i--; i >= 1; i--)
        {
            if (a[i] & 1L)
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
            if (a[i] & 1L)
                fmpz_mul(pow + i, t, pow + (i + 1));
            else
                fmpz_mul(pow + i, pow + (i + 1), pow + (i + 1));
        }

        /* Run Newton iteration */
        i = n - 1;
        {
            fmpz_mod(rop, op, pow + i);

            /* {(p-1) x^{p-2}}^{-1} */
            fmpz_sub(inv, p, rop);
        }
        for (i--; i >= 0; i--)
        {
            /* Lift rop */
            fmpz_powm(t, rop, pm1, pow + i);
            fmpz_sub_ui(t, t, 1);
            fmpz_mul(s, t, inv);
            fmpz_sub(rop, rop, s);
            fmpz_mod(rop, rop, pow + i);

            /* Lift inv */
            if (i > 0)
            {
                fmpz_powm(s, rop, pm2, pow + i);
                fmpz_mul(t, inv, pm1);
                fmpz_mul(r, s, t);
                fmpz_sub_ui(r, r, 2);
                fmpz_neg(r, r);
                fmpz_mul(inv, inv, r);
                fmpz_mod(inv, inv, pow + i);
            }
        }

        _fmpz_vec_clear(pow, n);

        fmpz_clear(r);
        fmpz_clear(s);
        fmpz_clear(t);
        fmpz_clear(inv);
        fmpz_clear(pm1);
        fmpz_clear(pm2);
    }
}

void padic_teichmuller(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (padic_val(op) < 0)
    {
        printf("ERROR (padic_teichmuller).  op is not a p-adic integer.\n");
        abort();
    }

    if (_padic_is_zero(op) || padic_val(op) > 0 || ctx->N <= 0)
    {
        padic_zero(rop);
        return;
    }

    _padic_teichmuller(padic_unit(rop), padic_unit(op), ctx->p, ctx->N);
    padic_val(rop) = 0;
}

