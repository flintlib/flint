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

    Copyright (C) 2009 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"

ulong
fmpz_tdiv_ui(const fmpz_t g, ulong h)
{
    fmpz c1 = *g;

    if (h == 0UL)
    {
        printf("Exception (fmpz_tdiv_ui). Division by 0.\n");
        abort();
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
		/* We need the absolut value of the remainder and
		   C 90 guarantees truncation towards zero. */
		if (c1 < 0L)
			return -c1 % h;
		else
			return c1 % h;
    }
    else                        /* g is large */
    {
        return mpz_tdiv_ui(COEFF_TO_PTR(c1), h);
    }
}
