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

int fmpz_poly_q_is_canonical(const fmpz_poly_q_t op)
{
    int ans;
    fmpz_poly_t t;

    if (fmpz_poly_is_zero(op->den))
        return 0;

    if (fmpz_sgn(fmpz_poly_lead(op->den)) < 0)
        return 0;

    fmpz_poly_init(t);
    fmpz_poly_gcd(t, op->num, op->den);
    ans = fmpz_poly_is_one(t);
    fmpz_poly_clear(t);
    return ans;
}

