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

void padic_get_fmpz(fmpz_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (padic_val(op) < 0)
    {
        flint_printf("Exception (padic_get_fmpz).  Negative valuation.\n");
        abort();
    }

    if (padic_is_zero(op))
    {
        fmpz_zero(rop);
    }
    else if (padic_val(op) == 0)
    {
        fmpz_set(rop, padic_unit(op));
    }
    else
    {
        fmpz_t pow;
        int alloc;

        alloc = _padic_ctx_pow_ui(pow, padic_val(op), ctx);
        fmpz_mul(rop, padic_unit(op), pow);

        if (alloc)
            fmpz_clear(pow);
    }
}

