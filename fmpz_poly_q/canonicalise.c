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

void fmpz_poly_q_canonicalise(fmpz_poly_q_t rop)
{
    fmpz_poly_t gcd;

    if (fmpz_poly_is_zero(rop->den))
    {
        flint_printf("Exception (fmpz_poly_q_canonicalise). Denominator is zero.\n");
        abort();
    }

    if (fmpz_poly_is_one(rop->den))
        return;

    fmpz_poly_init(gcd);
    fmpz_poly_gcd(gcd, rop->num, rop->den);
    if (!fmpz_poly_is_unit(gcd))
    {
        fmpz_poly_div(rop->num, rop->num, gcd);
        fmpz_poly_div(rop->den, rop->den, gcd);
    }
    fmpz_poly_clear(gcd);

    if (fmpz_sgn(fmpz_poly_lead(rop->den)) < 0)
    {
        fmpz_poly_neg(rop->num, rop->num);
        fmpz_poly_neg(rop->den, rop->den);
    }
}

