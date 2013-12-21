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

void fmpz_poly_q_scalar_div_mpz(fmpz_poly_q_t rop, 
                                const fmpz_poly_q_t op, const mpz_t x)
{
    fmpz_t y;

    if (mpz_sgn(x) == 0)
    {
        flint_printf("Exception (fmpz_poly_q_scalar_div_mpz). Division by zero.\n");
        abort();
    }

    fmpz_init(y);
    fmpz_set_mpz(y, x);

    fmpz_poly_set(rop->num, op->num);
    fmpz_poly_scalar_mul_fmpz(rop->den, op->den, y);
    fmpz_poly_q_canonicalise(rop);

    fmpz_clear(y);
}
