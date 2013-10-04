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

void fmpz_poly_q_scalar_div_mpq(fmpz_poly_q_t rop, 
                                const fmpz_poly_q_t op, const mpq_t x)
{
    fmpz_t num, den;

    if (mpz_sgn(mpq_numref(x)) == 0)
    {
        flint_printf("Exception (fmpz_poly_q_scalar_div_mpq). Division by zero.\n");
        abort();
    }

    fmpz_init(num);
    fmpz_init(den);
    fmpz_set_mpz(num, mpq_numref(x));
    fmpz_set_mpz(den, mpq_denref(x));

    fmpz_poly_scalar_mul_fmpz(rop->num, op->num, den);
    fmpz_poly_scalar_mul_fmpz(rop->den, op->den, num);
    fmpz_poly_q_canonicalise(rop);

    fmpz_clear(num);
    fmpz_clear(den);
}
