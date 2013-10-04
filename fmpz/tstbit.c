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

   Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include "fmpz.h"

int fmpz_tstbit(const fmpz_t f, ulong i)
{
    if (!COEFF_IS_MPZ(*f))
    {
        if (i < FLINT_BITS)
        {
            return ((WORD(1) << i) & *f) != 0;
        }
        else  /* i >= FLINT_BITS */
        {
            return *f < 0;
        }
    }
    else
    {
        return mpz_tstbit(COEFF_TO_PTR(*f), i);
    }
}

