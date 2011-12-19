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

void padic_pow_si(padic_t rop, const padic_t op, long e, const padic_ctx_t ctx)
{
    if (_padic_is_zero(op) || e * padic_val(op) >= ctx->N)
    {
        padic_zero(rop, ctx);
        return;
    }

    if (e)
    {
        fmpz_t pow;
        int alloc;

        if (e > 0)
        {
            padic_val(rop) = e * padic_val(op);

            alloc = _padic_ctx_pow_ui(pow, ctx->N - padic_val(rop), ctx);
            fmpz_powm_ui(padic_unit(rop), padic_unit(op), e, pow);
        }
        else if (e < 0)
        {
            /* u^{-1} to precision ceil((N - v e) / -e) */
            _padic_inv(padic_unit(rop), padic_unit(op), 
                       ctx->p, (ctx->N - padic_val(op) * e + (-e - 1)) / -e);

            padic_val(rop) = e * padic_val(op);

            alloc = _padic_ctx_pow_ui(pow, ctx->N - padic_val(rop), ctx);
            fmpz_powm_ui(padic_unit(rop), padic_unit(rop), -e, pow);
        }

        if (alloc)
            fmpz_clear(pow);
    }
    else
    {
        padic_one(rop, ctx);
    }
}

