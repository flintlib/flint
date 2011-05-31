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

    Copyright (C) 2011 Sebastian Pancratz
 
******************************************************************************/

#include "padic.h"

void _padic_inv_naive(fmpz_t rop, const fmpz_t op, const fmpz_t p, long N)
{
    fmpz_t pow;

    fmpz_init(pow);
    fmpz_pow_ui(pow, p, N);
    fmpz_invmod(rop, op, pow);
    fmpz_clear(pow);
}

void padic_inv_naive(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    fmpz_t pow;
    
    if (_padic_is_zero(op))
    {
        printf("Exception (padic_inv_naive).  Zero is not invertible.\n");
        abort();
    }

    /*
        If x = u p^v has negative valuation with N <= -v then its 
        exact inverse is equal to zero when reduced modulo p^N.
     */
    if (ctx->N + padic_val(op) <= 0)
    {
        padic_zero(rop, ctx);
        return;
    }

    _padic_inv_naive(padic_unit(rop), 
                     padic_unit(op), ctx->p, ctx->N + padic_val(op));

    padic_val(rop) = - padic_val(op);
}

