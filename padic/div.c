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

void padic_div(padic_t rop, const padic_t op1, const padic_t op2, 
               const padic_ctx_t ctx)
{
    if (padic_is_zero(op2))
    {
        flint_printf("Exception (padic_div).  op2 is zero.\n");
        abort();
    }

    if (padic_is_zero(op1) || padic_val(op1) - padic_val(op2) >= padic_prec(rop))
    {
        padic_zero(rop);
    }
    else
    {
        padic_t inv;

        padic_init(inv);

        _padic_inv(padic_unit(inv), padic_unit(op2), ctx->p, 
                   padic_prec(rop) - padic_val(op1) + padic_val(op2));
        padic_val(inv) = - padic_val(op2);
        padic_mul(rop, op1, inv, ctx);

        padic_clear(inv);
    }
}

