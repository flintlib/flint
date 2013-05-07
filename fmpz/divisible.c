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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int fmpz_divisible(const fmpz_t x, const fmpz_t p)
{
    fmpz y = *x;
    fmpz q = *p;

    if (y == 0L)
    {
        return 1;
    }

    if (!COEFF_IS_MPZ(y))
    {
        if (!COEFF_IS_MPZ(q))
        {
            return !(y % q);
        }
        else
        {
            return 0;
        }
    }
    else
    {
        if (!COEFF_IS_MPZ(q))
        {
            return mpz_divisible_ui_p(COEFF_TO_PTR(y), q);
        }
        else
        {
            return mpz_divisible_p(COEFF_TO_PTR(y), COEFF_TO_PTR(q));
        }
    }
}

