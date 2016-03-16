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

    Copyright (C) 2010, 2011 Sebastian Pancratz
   
******************************************************************************/

#include "fmpz_poly_q.h"

void fmpz_poly_q_inv(fmpz_poly_q_t rop, const fmpz_poly_q_t op)
{
    if (fmpz_poly_is_zero(op->num))
    {
        flint_printf("Exception (fmpz_poly_q_inv). Zero is not invertible.\n");
        flint_abort();
    }
    
    if (rop == op)
    {
        fmpz_poly_swap(rop->num, rop->den);
        if (fmpz_sgn(fmpz_poly_lead(rop->den)) < 0)
        {
            fmpz_poly_neg(rop->num, rop->num);
            fmpz_poly_neg(rop->den, rop->den);
        }
    }
    else
    {
        if (fmpz_sgn(fmpz_poly_lead(op->num)) > 0)
        {
            fmpz_poly_set(rop->num, op->den);
            fmpz_poly_set(rop->den, op->num);
        }
        else
        {
            fmpz_poly_neg(rop->num, op->den);
            fmpz_poly_neg(rop->den, op->num);
        }
    }
}

