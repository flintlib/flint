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

void padic_pow_si(padic_t rop, const padic_t op, len_t e, const padic_ctx_t ctx)
{
    if (e == 0)
    {
        padic_one(rop);
    }
    else if (padic_is_zero(op) || e * padic_val(op) >= padic_prec(rop))
    {
        padic_zero(rop);
    }
    else
    {
        fmpz_t pow;
        int alloc;

        padic_val(rop) = e * padic_val(op);

        alloc = _padic_ctx_pow_ui(pow, padic_prec(rop) - padic_val(rop), ctx);
        if (e > 0)
        {
            fmpz_powm_ui(padic_unit(rop), padic_unit(op), e, pow);
        }
        else  /* e < 0 */
        {
            /* u^{-1} to precision ceil((N - v e) / -e) */
            _padic_inv(padic_unit(rop), padic_unit(op), 
                       ctx->p, (padic_prec(rop) - padic_val(rop) + (-e - 1)) / -e);

            fmpz_powm_ui(padic_unit(rop), padic_unit(rop), -e, pow);
        }
        if (alloc)
            fmpz_clear(pow);
    }
}

