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

    Copyright (C) 2011 Jan Tuitman
    Copyright (C) 2011, 2012 Sebastian Pancratz

******************************************************************************/

#include "padic.h"

#define n    (S->n)
#define pow  (S->pow)

void _padic_inv_precompute(padic_inv_t S, const fmpz_t p, len_t N)
{
    len_t *a;

    a = _padic_lifts_exps(&n, N);

    pow = _fmpz_vec_init(2 * n + 2);

    _padic_lifts_pows(pow, a, n, p);

    flint_free(a);
}

void _padic_inv_clear(padic_inv_t S)
{
    _fmpz_vec_clear(pow, 2 * n + 2);
}

void _padic_inv_precomp(fmpz_t rop, const fmpz_t op, padic_inv_t S)
{
    len_t i;
    fmpz *t, *u;

    u = pow + n;
    t = pow + 2 * n;

    /* Compute reduced units */
    {
        fmpz_mod(u + 0, op, pow + 0);
    }
    for (i = 1; i < n; i++)
    {
        fmpz_mod(u + i, u + (i - 1), pow + i);
    }

    /* Run Newton iteration */
    i = n - 1;
    {
        fmpz_invmod(rop, u + i, pow + i);
    }
    for (i--; i >= 0; i--)
    {
        fmpz_mul(t, rop, rop);
        fmpz_mul(t + 1, u + i, t);
        fmpz_mul_2exp(rop, rop, 1);
        fmpz_sub(rop, rop, t + 1);
        fmpz_mod(rop, rop, pow + i);
    }
}

#undef n
#undef pow

void _padic_inv(fmpz_t rop, const fmpz_t op, const fmpz_t p, len_t N)
{
    if (N == 1)
    {
        fmpz_invmod(rop, op, p);
    }
    else
    {
        padic_inv_t S;

        _padic_inv_precompute(S, p, N);
        _padic_inv_precomp(rop, op, S);
        _padic_inv_clear(S);
    }
}

void padic_inv(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (padic_is_zero(op))
    {
        printf("Exception (padic_inv).  Zero is not invertible.\n");
        abort();
    }

    /*
        If x = u p^v has negative valuation with N <= -v then the 
        exact inverse of x is zero when reduced modulo $p^N$
     */
    if (padic_prec(rop) + padic_val(op) <= 0)
    {
        padic_zero(rop);
    }
    else
    {
        _padic_inv(padic_unit(rop), 
                   padic_unit(op), ctx->p, padic_prec(rop) + padic_val(op));

        padic_val(rop) = - padic_val(op);
    }
}

