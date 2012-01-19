/*============================================================================

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

    Copyright (C) 2012 Sebastian Pancratz
 
******************************************************************************/

#include "padic_poly.h"

void padic_poly_shift_left(padic_poly_t rop, const padic_poly_t op, long n)
{
    if (n == 0)
    {
        padic_poly_set(rop, op);
    }
    else if (op->length == 0)
    {
        padic_poly_zero(rop);
    }
    else
    {
        padic_poly_fit_length(rop, op->length + n);
        _fmpz_poly_shift_left(rop->coeffs, op->coeffs, op->length, n);
        rop->val = op->val;
        _padic_poly_set_length(rop, op->length + n);
    }
}

