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

void fmpz_poly_q_submul(fmpz_poly_q_t rop, 
                        const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)
{
    fmpz_poly_q_t temp;

    fmpz_poly_q_init(temp);
    fmpz_poly_q_mul(temp, op1, op2);
    fmpz_poly_q_sub(rop, rop, temp);
    fmpz_poly_q_clear(temp);
}
