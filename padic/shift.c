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

void padic_shift(padic_t rop, const padic_t op, len_t v, const padic_ctx_t ctx)
{
    if (padic_is_zero(op) || (padic_val(op) + v >= padic_prec(rop)))
    {
        padic_zero(rop);
    }
    else
    {
        fmpz_set(padic_unit(rop), padic_unit(op));
        padic_val(rop) = padic_val(op) + v;
        _padic_reduce(rop, ctx);
    }
}

