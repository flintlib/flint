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

void _padic_teichmuller(fmpz_t rop, const fmpz_t op, const fmpz_t p, len_t N)
{
    if (fmpz_equal_ui(p, 2))
    {
        fmpz_one(rop);
    }
    else if (N == 1)
    {
        fmpz_mod(rop, op, p);
    }
    else
    {
        len_t *a, i, n;
        fmpz *pow, *u;
        fmpz_t s, t;
        fmpz_t inv;
        fmpz_t pm1;

        a = _padic_lifts_exps(&n, N);

        pow = _fmpz_vec_init(2 * n);
        u   = pow + n;

        _padic_lifts_pows(pow, a, n, p);

        fmpz_init(s);
        fmpz_init(t);
        fmpz_init(inv);
        fmpz_init(pm1);

        fmpz_sub_ui(pm1, p, 1);

        /* Compute reduced units for (p-1) */
        {
            fmpz_mod(u + 0, pm1, pow + 0);
        }
        for (i = 1; i < n; i++)
        {
            fmpz_mod(u + i, u + (i - 1), pow + i);
        }

        /* Run Newton iteration */
        i = n - 1;
        {
            fmpz_mod(rop, op, pow + i);

            fmpz_set(inv, pm1);
        }
        for (i--; i >= 0; i--)
        {
            /* Lift rop */
            fmpz_powm(s, rop, p, pow + i);
            fmpz_sub(s, s, rop);
            fmpz_mul(t, s, inv);
            fmpz_sub(rop, rop, t);
            fmpz_mod(rop, rop, pow + i);

            /* Lift inv */
            if (i > 0)
            {
                fmpz_mul(s, inv, inv);
                fmpz_mul(t, u + i, s);
                fmpz_mul_2exp(inv, inv, 1);
                fmpz_sub(inv, inv, t);
                fmpz_mod(inv, inv, pow + i);
            }
        }

        fmpz_clear(s);
        fmpz_clear(t);
        fmpz_clear(inv);
        fmpz_clear(pm1);
        _fmpz_vec_clear(pow, 2 * n);
        flint_free(a);
    }
}

void padic_teichmuller(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (padic_val(op) < 0)
    {
        printf("Exception (padic_teichmuller).  op is not a p-adic integer.\n");
        abort();
    }

    if (padic_is_zero(op) || padic_val(op) > 0 || padic_prec(rop) <= 0)
    {
        padic_zero(rop);
    }
    else
    {
        _padic_teichmuller(padic_unit(rop), padic_unit(op), ctx->p, padic_prec(rop));
        padic_val(rop) = 0;
    }
}

