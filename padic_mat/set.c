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

    Copyright (C) 2011, 2013 Sebastian Pancratz

******************************************************************************/

#include "fmpz_mat.h"
#include "padic_mat.h"

void padic_mat_set(padic_mat_t rop, const padic_mat_t op, const padic_ctx_t ctx)
{
    if (op != rop)
    {
        if (padic_mat_val(op) >= padic_mat_prec(rop))
        {
            padic_mat_zero(rop);
        }
        else if (padic_mat_prec(rop) >= padic_mat_prec(op))
        {
            fmpz_mat_set(padic_mat(rop), padic_mat(op));
            padic_mat_val(rop) = padic_mat_val(op);
        }
        else
        {
            fmpz_mat_set(padic_mat(rop), padic_mat(op));
            padic_mat_val(rop) = padic_mat_val(op);

            _padic_mat_reduce(rop, ctx);
        }
    }
}

