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

int padic_div_exact_fmpz(padic_t rop, const padic_t op1, const fmpz_t op2,
                                   const padic_ctx_t ctx)
{
    int res = 1;

    if (fmpz_is_zero(op2))
    {
        flint_printf("Exception (padic_div_exact_fmpz).  op2 is zero.\n");
        abort();
    }

    if (fmpz_sgn(op2) < 0)
       return 0;

    if (padic_is_zero(op1))
    {
        padic_zero(rop);
    }
    else if (!fmpz_is_one(op2))
    {
        fmpz_t r;
        slong val;

        fmpz_init(r);
       
        val = fmpz_remove(r, op2, ctx->p);

        fmpz_tdiv_qr(padic_unit(rop), r, padic_unit(op1), r);

        res = fmpz_is_zero(r);

        fmpz_clear(r);

        padic_val(rop) = padic_val(op1) - val;
    }

    return res;
}

