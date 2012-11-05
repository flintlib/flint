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

void padic_get_fmpq(fmpq_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (padic_is_zero(op))
    {
        fmpq_zero(rop);
    }
    else
    {
        fmpz_t pow;
        int alloc = 0;

        if (padic_val(op) == 0)
        {
            fmpz_set(fmpq_numref(rop), padic_unit(op));
            fmpz_one(fmpq_denref(rop));
        }
        else if (padic_val(op) > 0)
        {
            alloc = _padic_ctx_pow_ui(pow, padic_val(op), ctx);
            fmpz_mul(fmpq_numref(rop), padic_unit(op), pow);
            fmpz_one(fmpq_denref(rop));
        }
        else  /* padic_val(op) < 0 */
        {
            alloc = _padic_ctx_pow_ui(pow, - padic_val(op), ctx);
            fmpz_set(fmpq_numref(rop), padic_unit(op));
            fmpz_set(fmpq_denref(rop), pow);
        }

        if (alloc)
            fmpz_clear(pow);
    }
}

